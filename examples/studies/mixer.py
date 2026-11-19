"""The baffled stirred vessel of gallery 96 and 98 to 102: the moving reference
frame (steady) and the sliding mesh by imprinting (three revolutions from
rest), laminar and with water and the k-ω SST model, the transient frozen
rotor and the vessel with two stirrers, with the field figures of the
gallery and the README.

    pip install gmsh matplotlib h5py numpy pillow
    cargo build --release -p foehncraft-cli
    python3 examples/studies/mixer.py [OUT] [--figures DIR] [--size 2.5e-3] [--force]
                                      [--turbulent | --rotors] [--blades 5]

The figures show the vessel with five blades and three baffles
(`--blades 5`, the default; meshes mixer5_*.msh and twin_mixer5.msh), never
a four-fold rotor pattern (CONTRIBUTING.md, "Images"). `--blades 4` runs the
four-blade vessel of the case files (the torque numbers of gallery 101 and
102) into OUT/four-blades and draws no figures.

Builds the two meshes with examples/meshes/mixer.py (when they are missing
or `--size` is given), runs examples/studies/mixer/mrf.toml into OUT/mrf
(16 300 cells, about 20 s) and sliding.toml into OUT/sliding (1 508 steps,
about seven minutes on four cores; a finished run is not run again unless
`--force`), and draws with tools/plot/rotor_fields.py into DIR (default OUT):

- mixer_mrf.png: |u| with streamlines, |u − Ω×r| in the rotor and ω_z of
  the MRF run;
- mixer_sliding.png: the same after three revolutions of the sliding mesh;
- mixer_sliding_startup.png: ω_z at t = 0.15, 0.3, 0.6 s and after three
  revolutions;
- mixer_sliding_interface.png: the imprinted interface at the tip of a blade,
  with the cell edges;
- mixer_relative.png: the fields of rotating machinery of the MRF run (D183,
  gallery 99): |w| with the streamlines of w (relative in the rotor, absolute
  in the stator), the relative flow angle β and the rothalpy in the rotor.

With `--turbulent` it runs mrf_turbulent.toml into OUT/mrf_turbulent (water,
Re = 20 000, k-ω SST, the coupled solver, about 30 s) and
sliding_turbulent.toml into OUT/sliding_turbulent (1 508 steps, about ten
minutes on four cores) instead, and draws (gallery 100):

- mixer_turbulent_mrf.png: |u| with streamlines, k and ν_t/ν of the frame run;
- mixer_turbulent_sliding.png: the same after three revolutions of the
  sliding mesh;
- stirred_vessel_mrf.jpg: the README tile of the frame run.

With `--rotors` it runs mrf_transient.toml into OUT/mrf_transient (the
transient frozen rotor, 480 steps) and examples/studies/twin_mixer/sliding.toml
on twin_mixer5.msh (three and five blades) into OUT/twin, after the
laminar sliding run, and draws (gallery 102 and the README):

- mixer_frozen_vs_sliding.png: ω_z of the frozen rotor and the sliding mesh
  at t = 0.3 and 0.6 s;
- twin_mixer_vorticity.png: ω_z of the twin vessel at t = 0.2 to 0.8 s;
- twin_mixer_speed.png: |u| with streamlines and ω_z at t = 0.8 s;
- twin_rotors_vorticity.jpg: a README tile (600 × 400; with `--readme DIR`
  written there; `--turbulent` writes the other one, stirred_vessel_mrf.jpg).

OUT defaults to output/mixer. Run from anywhere; paths are relative to the
repository.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT / "examples" / "studies" / "mixer"
MESHES = ROOT / "examples" / "meshes"
OMEGA, R_IFACE, R_VESSEL, NU_WATER = 10.0, 0.06, 0.10, 1e-6
GREY = (59, 59, 59)    # the colour of the solid parts in rotor_fields.py


def mesh_set(blades, split):
    """The `--set` that points a case of this folder at the mesh of `blades`."""
    stem = "mixer" if blades == 4 else f"mixer{blades}"
    kind = "split" if split else "conformal"
    return ["--set", f"mesh.file=../../meshes/{stem}_{kind}.msh"]


def run(cmd):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], cwd=ROOT, check=True)


def plot(out, fields, *opts):
    run([sys.executable, "-m", "tools.plot.rotor_fields", out, *fields,
         "--rotor-radius", R_IFACE, "--palette", *opts])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", nargs="?", default=str(ROOT / "output" / "mixer"))
    ap.add_argument("--figures", default=None, help="folder of the figures (default OUT)")
    ap.add_argument("--size", type=float, default=None, help="largest cell size (m), remeshes")
    ap.add_argument("--force", action="store_true", help="run the cases again")
    ap.add_argument("--turbulent", action="store_true",
                    help="the turbulent cases (water, k-ω SST) of gallery 100")
    ap.add_argument("--rotors", action="store_true",
                    help="the frozen rotor and the twin vessel of gallery 102")
    ap.add_argument("--blades", type=int, default=5, choices=(4, 5),
                    help="blades of the vessel (5: three baffles, the figures)")
    ap.add_argument("--readme", default=None, help="folder of the README tiles")
    a = ap.parse_args()
    out = Path(a.out).resolve()
    if a.blades == 4:
        out = out / "four-blades"
    fig = Path(a.figures).resolve() if a.figures else out
    fig.mkdir(parents=True, exist_ok=True)
    foehn = ROOT / "target" / "release" / ("foehn.exe" if os.name == "nt" else "foehn")
    if not foehn.is_file():
        raise SystemExit(f"{foehn} is missing: cargo build --release -p foehncraft-cli")

    stem = "mixer" if a.blades == 4 else f"mixer{a.blades}"
    meshes = [MESHES / f"{stem}_conformal.msh", MESHES / f"{stem}_split.msh"]
    if a.size is not None or not all(m.is_file() for m in meshes):
        shape = [] if a.blades == 4 else ["--blades", a.blades, "--baffles", 3]
        run([sys.executable, MESHES / "mixer.py"] + shape
            + (["--size", a.size] if a.size else []))
    if a.turbulent:
        turbulent(a, out, fig, foehn)
        return
    for name in ("mrf", "sliding"):
        if a.force or not (out / name / "result.json").is_file():
            run([foehn, "flow", CASES / f"{name}.toml", "--out", out / name]
                + mesh_set(a.blades, name == "sliding"))
    if a.blades == 4:
        print(f"four-blade runs in {out}; no figures (CONTRIBUTING.md, Images)")
        return
    if a.rotors:
        rotors(a, out, fig, foehn)
        return

    mrf, sld = out / "mrf", out / "sliding"
    three = [sld / "fields.vtkhdf"]
    common = ["--omega", OMEGA, "--outer-radius", R_VESSEL, "--dpi", 110,
              "--quantities", "speed,relative,vorticity"]
    plot(fig / "mixer_mrf.png", [mrf / "fields.vtkhdf"], *common,
         "--title", "Stirred vessel, moving reference frame (steady)")
    plot(fig / "mixer_sliding.png", three, *common,
         "--title", "Stirred vessel, sliding mesh (imprinting), after 3 revolutions")
    series = [sld / f"fields_{n:06d}.vtkhdf" for n in (120, 240, 480)] + three
    plot(fig / "mixer_sliding_startup.png", series, "--quantities", "vorticity",
         "--limits=-40:40", "--dpi", 100, "--title", "Start-up from rest: vorticity",
         "--labels", "t = 0.15 s,t = 0.3 s,t = 0.6 s,t = 1.89 s (3 rev)")
    plot(fig / "mixer_relative.png", [mrf / "fields.vtkhdf"], "--outer-radius", R_VESSEL,
         "--dpi", 110, "--quantities",
         "relative_velocity_magnitude,relative_flow_angle@rotor,rothalpy@rotor",
         "--title", "Stirred vessel (MRF): the relative flow")
    plot(fig / "mixer_sliding_interface.png", three, "--quantities", "vorticity",
         "--limits=-30:30", "--mesh", "--window", "0.04,0.07,-0.015,0.015", "--dpi", 150,
         "--labels", "Imprinted interface: 64 rotor / 96 stator intervals")


def turbulent(a, out, fig, foehn):
    """The turbulent runs and figures of gallery 100."""
    for name in ("mrf_turbulent", "sliding_turbulent"):
        if a.force or not (out / name / "result.json").is_file():
            run([foehn, "flow", CASES / f"{name}.toml", "--out", out / name]
                + mesh_set(a.blades, name.startswith("sliding")))
    if a.blades == 4:
        return
    common = ["--omega", OMEGA, "--outer-radius", R_VESSEL, "--nu", NU_WATER,
              "--dpi", 110, "--quantities", "speed,k,viscosity_ratio",
              "--limits", "0:0.45,0:0.005,0:150"]
    plot(fig / "mixer_turbulent_mrf.png", [out / "mrf_turbulent" / "fields.vtkhdf"], *common,
         "--title", "Stirred vessel, water (Re 20 000), k-ω SST: moving reference frame (steady)")
    plot(fig / "mixer_turbulent_sliding.png", [out / "sliding_turbulent" / "fields.vtkhdf"],
         *common, "--title",
         "Stirred vessel, water (Re 20 000), k-ω SST: sliding mesh, after 3 revolutions")
    readme = Path(a.readme).resolve() if a.readme else fig
    tile(out / "mrf_turbulent" / "fields.vtkhdf", readme / "stirred_vessel_mrf.jpg",
         "--quantities", "speed", "--outer-radius", R_VESSEL, "--limits", "0:0.45")



def tile(fields, dst, *opts):
    """A 600 × 400 README tile: one panel of `fields` drawn at high
    resolution by rotor_fields.py, cut to the axes (no title, no colour
    bar), cut (wide axes) or padded with grey (square axes) to 3:2."""
    from PIL import Image
    tmp = Path(dst).with_suffix(".tile.png")
    run([sys.executable, "-m", "tools.plot.rotor_fields", tmp, fields, "--dpi", 300, *opts])
    im = Image.open(tmp).convert("RGB")
    tmp.unlink()
    px = np.asarray(im).astype(int)
    # the axes: the box of the grey outside the vessel (left of the colour bar)
    grey = (np.abs(px - np.array(GREY)).sum(axis=2) < 12)
    grey[:, int(0.85 * px.shape[1]):] = False
    ys = np.where(grey.sum(axis=1) > 20)[0]      # not the grey of the text
    xs = np.where(grey.sum(axis=0) > 20)[0]
    box = im.crop((xs.min() + 10, ys.min() + 10, xs.max() - 10, ys.max() - 10))
    w, h = box.size
    if w / h > 1.5:
        cut = int((w - 1.5 * h) / 2)
        box = box.crop((cut, 0, w - cut, h))
    else:
        canvas = Image.new("RGB", (int(1.5 * h), h), GREY)
        canvas.paste(box, ((canvas.width - w) // 2, 0))
        box = canvas
    box.resize((600, 400), Image.LANCZOS).save(dst, quality=88)


def rotors(a, out, fig, foehn):
    """The frozen rotor against the sliding mesh, and the twin vessel with
    three and five blades (gallery 102, README tiles)."""
    if a.force or not (out / "mrf_transient" / "result.json").is_file():
        run([foehn, "flow", CASES / "mrf_transient.toml", "--out", out / "mrf_transient"]
            + mesh_set(a.blades, False))
    twin = MESHES / "twin_mixer5.msh"
    if a.size is not None or not twin.is_file():
        run([sys.executable, MESHES / "twin_mixer.py", "--right-blades", 5]
            + (["--size", a.size] if a.size else []))
    if a.force or not (out / "twin" / "result.json").is_file():
        run([foehn, "flow", ROOT / "examples" / "studies" / "twin_mixer" / "sliding.toml",
             "--out", out / "twin", "--set", "mesh.file=../../meshes/twin_mixer5.msh"])
    frz, sld, tw = out / "mrf_transient", out / "sliding", out / "twin"
    plot(fig / "mixer_frozen_vs_sliding.png",
         [frz / "fields_000240.vtkhdf", sld / "fields_000240.vtkhdf",
          frz / "fields_000480.vtkhdf", sld / "fields_000480.vtkhdf"],
         "--quantities", "vorticity", "--limits=-60:60", "--dpi", 100,
         "--labels", "frozen rotor t = 0.3 s,sliding mesh t = 0.3 s,"
         "frozen rotor t = 0.6 s,sliding mesh t = 0.6 s")
    series = [tw / f"fields_{n:06d}.vtkhdf" for n in (200, 400, 600)] + [tw / "fields.vtkhdf"]
    run([sys.executable, "-m", "tools.plot.rotor_fields", fig / "twin_mixer_vorticity.png",
         *series, "--quantities", "vorticity", "--limits=-60:60", "--dpi", 100, "--palette",
         "--labels", "t = 0.2 s,t = 0.4 s,t = 0.6 s,t = 0.8 s"])
    run([sys.executable, "-m", "tools.plot.rotor_fields", fig / "twin_mixer_speed.png",
         tw / "fields.vtkhdf", "--quantities", "speed,vorticity", "--limits", "0:0.7,-60:60",
         "--dpi", 100, "--palette"])
    readme = Path(a.readme).resolve() if a.readme else fig
    tile(tw / "fields.vtkhdf", readme / "twin_rotors_vorticity.jpg",
         "--quantities", "vorticity", "--limits=-60:60")


if __name__ == "__main__":
    main()
