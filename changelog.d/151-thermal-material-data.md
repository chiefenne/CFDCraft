- **Materials: thermal data (stage 7.0c).** Every material may carry a
  specific heat c_p and a thermal conductivity k (any property form, or
  orthotropic with two axes). Solids need density, c_p and k (MAT006).
  User material database files (`[materials] include`, `FOEHN_MATERIALS`)
  extend the built-in database; `foehn materials --include FILE` lists them.
  Water and air have c_p and k tables from the IAPWS-95 and IAPWS 2011 (water)
  and Helmholtz-equation (air) formulations, evaluated with CoolProp 8.0 (MIT).
  NASA-7 and NASA-9 specific-heat forms and a CHEMKIN reader are in.
