- **Kernels: density per cell and per face (stage 7.0b).** The flux,
  momentum, scalar and turbulence kernels read ρ from cell and face fields
  instead of one number, on the CPU and the GPU. Results are bit for bit
  unchanged; the fields carry the variable density of stage 7.
