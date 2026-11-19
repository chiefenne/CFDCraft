"""The baffled stirred vessel of gallery 98 and 99: the moving reference frame
(steady) and the sliding mesh by imprinting (three revolutions from rest),
with the field figures of the gallery.

    pip install gmsh matplotlib h5py numpy pillow
    cargo build --release -p foehncraft-cli
    python3 examples/studies/mixer.py [OUT] [--figures DIR] [--size 2.5e-3] [--force]

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

OUT defaults to output/mixer. Run from anywhere; paths are relative to the
repository.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT / "examples" / "studies" / "mixer"
MESHES = ROOT / "examples" / "meshes"
OMEGA, R_IFACE, R_VESSEL = 10.0, 0.06, 0.10


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
    a = ap.parse_args()
    out = Path(a.out).resolve()
    fig = Path(a.figures).resolve() if a.figures else out
    fig.mkdir(parents=True, exist_ok=True)
    foehn = ROOT / "target" / "release" / ("foehn.exe" if os.name == "nt" else "foehn")
    if not foehn.is_file():
        raise SystemExit(f"{foehn} is missing: cargo build --release -p foehncraft-cli")

    meshes = [MESHES / "mixer_conformal.msh", MESHES / "mixer_split.msh"]
    if a.size is not None or not all(m.is_file() for m in meshes):
        run([sys.executable, MESHES / "mixer.py"] + (["--size", a.size] if a.size else []))
    for name in ("mrf", "sliding"):
        if a.force or not (out / name / "result.json").is_file():
            run([foehn, "flow", CASES / f"{name}.toml", "--out", out / name])

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


if __name__ == "__main__":
    main()
