- **Torque, shaft power and machine performance of rotating machinery
  (D184, gallery 101, recipe 101).** For `[[mrf]]` and `[[sliding]]`: the
  torque of the fluid on every wall about the axis of its rotor (the walls of
  a rotor are its `walls`, by default the walls inside its zone), in the
  boundary table, `boundaries.csv`, `forces.csv` and the report; the shaft
  power −Ω·T of each rotor, booked in the step and the final energy balance
  (`balance.csv` column `shaft_power`, the rows of the ledger tables) and in
  the compact balance line; the angular momentum balance about the axis of
  each rotor. `[[output.machines]]` names a rotor, an inlet and an outlet and
  gives the volume flow, the total pressure rise, the head, the hydraulic and
  total-to-static efficiency and the Euler work ṁ Ω Δ(r c_u). New log
  columns `rotor.<zone>.torque|power|angle` (the angle of every sliding rotor,
  in degrees, each with its own speed) and `machine.<name>.flow|pressure_rise|head|efficiency`
  (all appear by default). The convected kinetic energy of the ledger uses the
  flux relative to a turning zone. The torque of a moving wall includes the
  transposed viscous stress: the Laplacian form of the momentum rows misses
  the −u_θ/r of the shear stress (a turning cylinder: 35 % too small before,
  0.5 % against the exact torque now). `tools/plot/rotor_movie.py` writes the
  frames of an animation of a field series in ParaView.
