- **The flow a blade sees: relative and frame velocities of rotating
  machinery (D183, gallery 99, recipe 99).** With an `[[mrf]]` or a
  `[[sliding]]` table, `output.derived` takes five new fields in the
  notation of turbomachinery (absolute **c**, blade speed **u** = Ω × r,
  relative **w** = **c** − **u**): `relative-velocity` (**w** in each turning
  zone, **c** elsewhere, so a pump shows the relative flow in the impeller
  and the absolute flow in the volute), `frame-velocity`, `turbomachinery`
  (c_r, c_u, c_a, c_m, w_u about the axis and the flow angles α and β),
  `relative-total-pressure` and `rothalpy`. `output.frame = "<zone>"` shows
  every cell in the frame of one zone. Sample lines get the columns of **w**
  and of the cylindrical components, point probes of a transient run those
  of **w**; on the sliding mesh at the current position of the cell.
  `tools/plot/rotor_fields.py` draws any dataset of the file, also in the
  rotor only (`NAME@rotor`).
- **Fix:** line samples and the reference comparison of the sliding mesh by
  imprinting located their points on the mesh of the start of the run, whose
  rotor cells had moved on; they use the mesh of the end of the run now.
