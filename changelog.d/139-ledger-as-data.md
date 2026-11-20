- **Ledger as data (D119, stage 7 part 7.0a).** The conservation ledger
  books every quantity as an account of named terms (quantity, unit,
  region, term kind, name, exchange-pair id) and closes it with one
  routine, so the energy module and later modules add their quantities
  without new structs. Mass and kinetic energy are booked this way; every
  report (`balance.csv`, log tables, `flow.md`, `result.json`) is bitwise
  as before.
