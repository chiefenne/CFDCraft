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
PANEL = (2000, 420)          # pixels of one model panel
U_RANGE = (-0.3, 1.1)


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
    cut.SliceType.Origin = [0.0, 0.0, 0.05]
    cut.SliceType.Normal = [0.0, 0.0, 1.0]
    points = CellDatatoPointData(Input=cut)
    points.PassCellData = 1
    vel = Calculator(Input=points)
    vel.ResultArrayName = "U"
    vel.Function = "velocity_x*iHat + velocity_y*jHat"
    vel.UpdatePipeline()
    return reader, vel


def seeded(src, p1, p2, n):
    st = StreamTracer(Input=src, SeedType="Line")
    st.Vectors = ["POINTS", "U"]
    st.SurfaceStreamlines = 1
    st.IntegrationDirection = "BOTH"
    st.MaximumStreamlineLength = 60.0
    st.SeedType.Point1 = p1
    st.SeedType.Point2 = p2
    st.SeedType.Resolution = n
    return st


def view():
    v = CreateView("RenderView")
    v.ViewSize = list(PANEL)
    v.OrientationAxesVisibility = 0
    v.Background = [1.0, 1.0, 1.0]
    try:
        v.UseColorPaletteForBackground = 0
    except AttributeError:
        pass
    v.CameraParallelProjection = 1
    # x from X0 to X1 plus room for the colour bar on the right
    xc, yc = 0.5 * (X0 + X1) + 1.2, 0.5 * (Y0 + Y1)
    v.CameraPosition = [xc, yc, 50.0]
    v.CameraFocalPoint = [xc, yc, 0.0]
    v.CameraViewUp = [0.0, 1.0, 0.0]
    v.CameraParallelScale = 0.5 * (Y1 - Y0) * 1.18
    v.AxesGrid.Visibility = 1
    v.AxesGrid.XTitle = "x/H"
    v.AxesGrid.YTitle = "y/H"
    v.AxesGrid.ZAxisVisibility = 0
    for a in ("XTitleColor", "YTitleColor", "XLabelColor", "YLabelColor", "GridColor"):
        setattr(v.AxesGrid, a, [0.0, 0.0, 0.0])
    v.AxesGrid.DataBoundsInflateFactor = 0.0
    v.AxesGrid.UseCustomBounds = 1
    v.AxesGrid.CustomBounds = [X0, X1, Y0, Y1, 0.0, 0.0]
    return v


def marker(v, x, colour):
    s = Sphere(Center=[x, 0.0, 0.1], Radius=0.09, ThetaResolution=24, PhiResolution=12)
    d = Show(s, v)
    d.AmbientColor = colour
    d.DiffuseColor = colour
    ColorBy(d, None)


def panel(v, data, label, xr, field, lut_name, rng, streamlines):
    d = Show(data, v)
    ColorBy(d, ("POINTS", field))
    lut = GetColorTransferFunction(field)
    lut.ApplyPreset(lut_name, True)
    lut.RescaleTransferFunction(*rng)
    try:
        lut.AutomaticRescaleRangeMode = "Never"
    except AttributeError:
        pass
    bar = GetScalarBar(lut, v)
    bar.Title = "U/U_ref" if field == "velocity_x" else "k/U_ref²"
    bar.ComponentTitle = ""
    bar.TitleColor = [0.0, 0.0, 0.0]
    bar.LabelColor = [0.0, 0.0, 0.0]
    bar.WindowLocation = "Any Location"
    bar.Position = [0.93, 0.12]
    bar.ScalarBarLength = 0.76
    d.SetScalarBarVisibility(v, True)

    # the step itself
    step = Box(XLength=4.0, YLength=1.0, ZLength=0.01, Center=[X0 - 1.0 + 2.0, 0.5, 0.06])
    sd = Show(step, v)
    ColorBy(sd, None)
    sd.AmbientColor = sd.DiffuseColor = [0.55, 0.55, 0.55]

    if streamlines:
        xb = 0.5 * (xr or X_EXP)
        for p1, p2, n in (([X0 + 0.05, 1.08, 0.05], [X0 + 0.05, 2.9, 0.05], 8),
                          ([xb, 0.06, 0.05], [xb, 0.75, 0.05], 5),
                          ([0.25, 0.08, 0.05], [0.15, 0.25, 0.05], 1)):
            st = seeded(data, p1, p2, n)
            sl = Show(st, v)
            ColorBy(sl, None)
            sl.AmbientColor = sl.DiffuseColor = [0.0, 0.0, 0.0]
            sl.LineWidth = 1.2

    edge = Contour(Input=data)
    edge.ContourBy = ["POINTS", "velocity_x"]
    edge.Isosurfaces = [0.0]
    ed = Show(edge, v)
    ColorBy(ed, None)
    ed.AmbientColor = ed.DiffuseColor = [1.0, 1.0, 1.0] if streamlines else [0.0, 0.8, 0.9]
    ed.LineWidth = 2.0

    marker(v, X_EXP, [0.0, 0.0, 0.0])
    if xr:
        marker(v, xr, [0.84, 0.15, 0.16])
    title = Text()
    xt = f"x_r/H = {xr:.2f}" if xr else "no reattachment"
    title.Text = f"{label}:  {xt}  (experiment {X_EXP})" if streamlines else label
    td = Show(title, v)
    td.Color = [0.0, 0.0, 0.0]
    td.FontSize = 18
    td.WindowLocation = "Upper Left Corner"
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
    for field, lut_name, rng, lines, name in (
        ("velocity_x", "Cool to Warm", U_RANGE, True, "velocity"),
        ("k", "Inferno (matplotlib)", (0.0, kmax), False, "k"),
    ):
        parts = []
        for i, (label, data, xr) in enumerate(models):
            v = view()
            panel(v, data, label, xr, field, lut_name, rng, lines)
            f = f"{prefix}_{name}_{i}.png"
            SaveScreenshot(f, v, ImageResolution=list(PANEL))
            parts.append(f)
            Delete(v)
        stack(parts, f"{prefix}_{name}.png")
        print("wrote", f"{prefix}_{name}.png")


if __name__ == "__main__":
    main()
