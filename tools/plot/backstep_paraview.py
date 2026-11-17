"""Backward-facing step of gallery 74 in ParaView: the streamwise velocity
with streamlines and the turbulent kinetic energy of several turbulence
models, one panel per model stacked into one picture per field, with the
bubble edge U = 0, the computed reattachment point (red) and the
experiment's 6.26 H (black). Run it with ParaView's Python (pvpython, or
the Python shell of the ParaView GUI) from the repository root:

    pvpython tools/plot/backstep_paraview.py docs/figures/gallery/backstep \\
        "k-ε (scalable)=output/backstep_keps" \\
        "k-ε-ζ-f₀=output/backstep_zeta_f0" \\
        "k-ω SST=output/backstep_sst"

(on macOS: /Applications/ParaView-6.2.0.app/Contents/bin/pvpython). DIR
is the output directory of a run of examples/studies/backstep/*.toml; the
script reads DIR/fields.vtkhdf and DIR/wall_ymin.csv and writes
OUT_PREFIX_velocity.png and OUT_PREFIX_k.png. It prints x_r/H per model.
tools/plot/backstep_chart.py draws the same figures with matplotlib.
"""
import csv
import os
import sys

from paraview.simple import (
    Box, CellDatatoPointData, Calculator, ColorBy, Contour, CreateView, Delete,
    ExtractBlock, GetColorTransferFunction, GetScalarBar, OpenDataFile,
    Render, SaveScreenshot, Show, Slice, Sphere, StreamTracer, Text,
)
from vtkmodules.vtkIOImage import vtkPNGReader, vtkPNGWriter
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
import numpy as np
from vtkmodules.vtkCommonDataModel import vtkImageData

X_EXP = 6.26                 # Driver & Seegmiller (1985), x_r/H
X0, X1, Y0, Y1 = -3.0, 14.0, 0.0, 3.0
PANEL = (2000, 480)          # pixels of one model panel
# the camera window: the field plus room for the axes, the title and the bar
VIEW = (X0 - 0.9, X1 + 1.9, None, None)
_H = (VIEW[1] - VIEW[0]) * PANEL[1] / PANEL[0]
VIEW = (VIEW[0], VIEW[1], Y0 - 0.5 * (_H - (Y1 - Y0)) - 0.2, Y1 + 0.5 * (_H - (Y1 - Y0)) - 0.2)
Z_CUT, Z_MASK, Z_AXES, Z_MARK = 0.05, 0.06, 0.08, 0.1   # the layers, bottom to top
U_RANGE = (-0.3, 1.1)


def opt(obj, name, value):
    """Set a display property that may not exist in every ParaView version."""
    try:
        setattr(obj, name, value)
    except Exception:  # noqa: BLE001  (also paraview.NotSupportedException)
        pass


def opt_fails(obj, name, value):
    """Like opt, but say whether the value was refused."""
    try:
        setattr(obj, name, value)
        return False
    except Exception:  # noqa: BLE001  (also paraview.NotSupportedException)
        return True


def solid(rep, colour):
    """One colour for a representation, whatever array it was coloured by."""
    try:
        ColorBy(rep, None)
    except Exception:  # noqa: BLE001  (ParaView refuses it when no array is set)
        opt(rep, "ColorArrayName", [None, ""])
    opt(rep, "AmbientColor", colour)
    opt(rep, "DiffuseColor", colour)
    opt(rep, "Ambient", 1.0)    # the colour as given, without shading
    opt(rep, "Diffuse", 0.0)


def safely(what, draw):
    """Draw one decoration; a ParaView version that refuses it costs only that."""
    try:
        draw()
    except Exception as e:  # noqa: BLE001
        print(f"  skipped {what}: {type(e).__name__}: {e}")


# the inferno map of matplotlib as (t, r, g, b): set as a table, because the
# preset's name differs between ParaView versions
INFERNO = (
    (0.0000, 0.001, 0.000, 0.014), (0.0625, 0.042, 0.028, 0.141), (0.1250, 0.129, 0.047, 0.291),
    (0.1875, 0.238, 0.037, 0.396), (0.2500, 0.342, 0.062, 0.429), (0.3125, 0.441, 0.099, 0.432),
    (0.3750, 0.541, 0.135, 0.415), (0.4375, 0.640, 0.171, 0.381), (0.5000, 0.736, 0.216, 0.330),
    (0.5625, 0.822, 0.275, 0.266), (0.6250, 0.894, 0.353, 0.194), (0.6875, 0.947, 0.449, 0.115),
    (0.7500, 0.978, 0.558, 0.035), (0.8125, 0.988, 0.675, 0.065), (0.8750, 0.975, 0.798, 0.206),
    (0.9375, 0.948, 0.917, 0.411), (1.0000, 0.988, 0.998, 0.645),
)


def colour_map(lut, preset, rng):
    """A preset by name, or (preset None) the inferno table above."""
    if preset:
        try:
            lut.ApplyPreset(preset, True)
        except Exception as e:  # noqa: BLE001
            print(f"  colour map {preset} not available ({e}); ParaView's default is used")
    else:
        lo, hi = rng
        lut.RGBPoints = [x for t, r, g, b in INFERNO for x in (lo + t * (hi - lo), r, g, b)]
        opt(lut, "ColorSpace", "RGB")
    lut.RescaleTransferFunction(*rng)


def reattachment(d):
    with open(os.path.join(d, "wall_ymin.csv")) as f:
        rows = [r for r in f if not r.startswith("#")]
    w = sorted((float(r["x"]), float(r["tau_x"])) for r in csv.DictReader(rows))
    xr = None
    for (x0, t0), (x1, t1) in zip(w, w[1:]):
        if x0 > 0 and t0 < 0 <= t1:
            xr = x0 + (x1 - x0) * (-t0) / (t1 - t0)
    return xr


def load(d):
    """The interior cells cut at mid-depth, with point data and U = (u, v)."""
    reader = OpenDataFile(os.path.join(d, "fields.vtkhdf"))
    interior = ExtractBlock(Input=reader)
    interior.Selectors = ["/Root/internal"]
    interior.UpdatePipeline()
    if "velocity_x" not in interior.CellData.keys():
        # an older reader without the assembly names: take everything
        interior = reader
    cut = Slice(Input=interior)
    cut.SliceType = "Plane"
    cut.SliceType.Origin = [0.0, 0.0, Z_CUT]
    cut.SliceType.Normal = [0.0, 0.0, 1.0]
    points = CellDatatoPointData(Input=cut)
    points.PassCellData = 1
    vel = Calculator(Input=points)
    vel.ResultArrayName = "U"
    vel.Function = "velocity_x*iHat + velocity_y*jHat"
    vel.UpdatePipeline()
    return reader, vel


def seeded(src, p1, p2, n, length):
    st = StreamTracer(Input=src, SeedType="Line")
    st.Vectors = ["POINTS", "U"]
    opt(st, "SurfaceStreamlines", 1)
    st.IntegrationDirection = "BOTH"
    st.MaximumStreamlineLength = length
    st.SeedType.Point1 = p1
    st.SeedType.Point2 = p2
    st.SeedType.Resolution = n
    return st


def unicode_font():
    """A font file with Greek letters and subscripts, or None."""
    files = []
    try:
        import matplotlib
        files.append(os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans.ttf"))
    except Exception:  # noqa: BLE001
        pass
    files += ["/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
              "/Library/Fonts/Arial Unicode.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
              "C:/Windows/Fonts/arial.ttf"]
    return next((f for f in files if os.path.exists(f)), None)


FONT = unicode_font()
PLAIN = {"ε": "epsilon", "ζ": "zeta", "ω": "omega", "₀": "0", "²": "^2"}


def text(s):
    """The string as it is with a Unicode font, in plain letters without one."""
    if FONT:
        return s
    for a, b in PLAIN.items():
        s = s.replace(a, b)
    return s


def font(obj, prefix=""):
    if FONT:
        opt(obj, prefix + "FontFamily", "File")
        opt(obj, prefix + "FontFile", FONT)


def view():
    v = CreateView("RenderView")
    v.ViewSize = list(PANEL)
    v.OrientationAxesVisibility = 0
    v.Background = [1.0, 1.0, 1.0]
    opt(v, "UseColorPaletteForBackground", 0)
    opt(v, "BackgroundColorMode", "Single Color")
    opt(v, "InteractionMode", "2D")
    g = getattr(v, "AxesGrid", None)
    if g is not None:
        for name, value in (("Visibility", 1), ("XTitle", "x/H"), ("YTitle", "y/H"),
                            ("ZTitle", ""), ("ShowEdges", 1), ("ShowTicks", 1),
                            ("ShowGrid", 0), ("AxesToLabel", 3), ("UseCustomBounds", 1),
                            ("CustomBounds", [X0, X1, Y0, Y1, Z_AXES, Z_AXES]),
                            ("XAxisUseCustomLabels", 1),
                            ("XAxisLabels", [float(x) for x in range(-2, 15, 2)]),
                            ("YAxisUseCustomLabels", 1),
                            ("YAxisLabels", [0.0, 1.0, 2.0, 3.0]),
                            ("XTitleFontSize", 14), ("YTitleFontSize", 14),
                            ("XLabelFontSize", 13), ("YLabelFontSize", 13)):
            opt(g, name, value)
        for a in ("XTitleColor", "YTitleColor", "XLabelColor", "YLabelColor",
                  "GridColor"):
            opt(g, a, [0.0, 0.0, 0.0])
    return v


def frame(v):
    """The camera on the window VIEW (parallel, looking down z)."""
    xc, yc = 0.5 * (VIEW[0] + VIEW[1]), 0.5 * (VIEW[2] + VIEW[3])
    v.CameraParallelProjection = 1
    v.CameraPosition = [xc, yc, 50.0]
    v.CameraFocalPoint = [xc, yc, 0.0]
    v.CameraViewUp = [0.0, 1.0, 0.0]
    v.CameraParallelScale = 0.5 * (VIEW[3] - VIEW[2])


def block(v, x0, x1, y0, y1, colour):
    """A flat rectangle above the field: the step, and masks outside the window."""
    b = Box(XLength=x1 - x0, YLength=y1 - y0, ZLength=0.004,
            Center=[0.5 * (x0 + x1), 0.5 * (y0 + y1), Z_MASK])
    solid(Show(b, v), colour)


def marker(v, x, colour):
    s = Sphere(Center=[x, 0.0, Z_MARK], Radius=0.09, ThetaResolution=24, PhiResolution=12)
    solid(Show(s, v), colour)


def panel(v, data, label, xr, field, preset, rng, streamlines):
    d = Show(data, v)
    ColorBy(d, ("POINTS", field))
    opt(d, "Ambient", 1.0)      # the colours of the map, without shading
    opt(d, "Diffuse", 0.0)
    lut = GetColorTransferFunction(field)
    colour_map(lut, preset, rng)
    opt(lut, "AutomaticRescaleRangeMode", "Never")
    bar = GetScalarBar(lut, v)
    for name, value in (("Title", text("U/U_ref" if field == "velocity_x" else "k/U_ref²")),
                        ("ComponentTitle", ""), ("TitleColor", [0.0, 0.0, 0.0]),
                        ("LabelColor", [0.0, 0.0, 0.0]), ("TitleFontSize", 14),
                        ("LabelFontSize", 13), ("AddRangeLabels", 0),
                        ("WindowLocation", "Any Location"), ("Orientation", "Vertical"),
                        ("Position", [0.925, 0.2]), ("ScalarBarLength", 0.6)):
        opt(bar, name, value)
    font(bar, "Title")
    d.SetScalarBarVisibility(v, True)

    white, far = [1.0, 1.0, 1.0], 200.0
    block(v, X0, 0.0, 0.0, 1.0, [0.55, 0.55, 0.55])        # the step
    block(v, X0 - far, X0, -far, far, white)               # left of the window
    block(v, X1, X1 + far, -far, far, white)               # right of it
    block(v, X0, X1, Y1, far, white)                       # above it

    def lines():
        xb = 0.5 * (xr or X_EXP)
        once = 1.1 * (xr or X_EXP)      # half the way around the bubble's edge
        seeds = [([X0 + 0.05, 1.08, Z_CUT], [X0 + 0.05, 2.9, Z_CUT], 8, 60.0),
                 ([0.25, 0.08, Z_CUT], [0.15, 0.25, Z_CUT], 1, 1.0)]
        # closed lines in the bubble, each about once around
        for y in (0.06, 0.15, 0.25, 0.35):
            seeds.append(([xb, y, Z_CUT], [xb, y, Z_CUT], 1, once * (1.0 - y / 0.5)))
        for p1, p2, n, length in seeds:
            sl = Show(seeded(data, p1, p2, n, length), v)
            solid(sl, [0.0, 0.0, 0.0])
            opt(sl, "LineWidth", 1.2)

    def bubble_edge():
        edge = Contour(Input=data)
        edge.ContourBy = ["POINTS", "velocity_x"]
        edge.Isosurfaces = [0.0]
        ed = Show(edge, v)
        solid(ed, [1.0, 1.0, 1.0] if streamlines else [0.0, 0.8, 0.9])
        opt(ed, "LineWidth", 2.0)

    def points():
        marker(v, X_EXP, [0.0, 0.0, 0.0])
        if xr:
            marker(v, xr, [0.84, 0.15, 0.16])

    def heading():
        title = Text()
        xt = f"x_r/H = {xr:.2f}" if xr else "no reattachment"
        title.Text = text(f"{label}:  {xt}  (experiment {X_EXP})" if streamlines else label)
        td = Show(title, v)
        opt(td, "Color", [0.0, 0.0, 0.0])
        opt(td, "FontSize", 16)
        font(td)
        if opt_fails(td, "WindowLocation", "Upper Left Corner"):
            opt(td, "WindowLocation", "UpperLeftCorner")

    if streamlines:
        safely("streamlines", lines)
    safely("the U = 0 line", bubble_edge)
    safely("the reattachment markers", points)
    safely("the title", heading)
    # ParaView resets the camera at the first render of a session: render,
    # then place the camera
    Render(v)
    frame(v)
    Render(v)


def stack(files, out):
    """One picture of the panels, the first on top (numpy, no VTK filters)."""
    arrays = []
    for f in files:
        r = vtkPNGReader()
        r.SetFileName(f)
        r.Update()
        im = r.GetOutput()
        nx, ny, _ = im.GetDimensions()
        nc = im.GetNumberOfScalarComponents()
        # VTK images start at the lower left: the first row is the bottom
        arrays.append(vtk_to_numpy(im.GetPointData().GetScalars()).reshape(ny, nx, nc))
    # the panel to be on top goes last in VTK's bottom-up order
    big = np.ascontiguousarray(np.concatenate(arrays[::-1], axis=0))
    img = vtkImageData()
    img.SetDimensions(big.shape[1], big.shape[0], 1)
    sc = numpy_to_vtk(big.reshape(-1, big.shape[2]), deep=True)
    img.GetPointData().SetScalars(sc)
    w = vtkPNGWriter()
    w.SetFileName(out)
    w.SetInputData(img)
    w.Write()
    for f in files:
        os.remove(f)


def main():
    prefix = sys.argv[1]
    runs = [a.split("=", 1) for a in sys.argv[2:]]
    models = []
    for label, d in runs:
        reader, data = load(d)
        xr = reattachment(d)
        models.append((label, data, xr))
        print(f"{label:22s} " + (f"x_r/H {xr:.2f} ({100 * (xr / X_EXP - 1):+.1f} % from {X_EXP})"
                                 if xr else "no reattachment"))
    kmax = max(m[1].PointData["k"].GetRange()[1] for m in models)
    for field, preset, rng, lines, name in (
        ("velocity_x", "Cool to Warm", U_RANGE, True, "velocity"),
        ("k", None, (0.0, kmax), False, "k"),
    ):
        parts = []
        for i, (label, data, xr) in enumerate(models):
            v = view()
            panel(v, data, label, xr, field, preset, rng, lines)
            f = f"{prefix}_{name}_{i}.png"
            SaveScreenshot(f, v, ImageResolution=list(PANEL))
            parts.append(f)
            Delete(v)
        stack(parts, f"{prefix}_{name}.png")
        print("wrote", os.path.abspath(f"{prefix}_{name}.png"))


if __name__ == "__main__":
    main()
