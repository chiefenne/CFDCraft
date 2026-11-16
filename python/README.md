# foehncraft (Python)

The Python API of [Foehncraft](https://github.com/chiefenne/Foehncraft)
(D61): read, change and write case files, build expressions with Python
operators, run cases in process or through the `foehn` command line with
callbacks at every iteration or time step, read results and fields, query
the material database, and run parametric studies.

```
pip install ./python            # builds the native module (Rust, maturin)
python -m pytest python/tests
```

See the [gallery](../docs/gallery/python.md) and the
[recipe](../docs/recipes/python.md).
