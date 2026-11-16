#!/usr/bin/env bash
# D87: the flow core without default features depends on no C library.
# Fails if the normal dependency tree of foehncraft-flow (no HDF5) holds a
# `-sys` crate or a crate that builds C code (`cc`, `cmake`).
set -euo pipefail
tree=$(cargo tree -p foehncraft-flow --no-default-features -e normal,build --prefix none --format '{p}' | sort -u)
bad=$(printf '%s\n' "$tree" | grep -E '^([A-Za-z0-9_-]+-sys|cc|cmake|pkg-config) ' || true)
if [ -n "$bad" ]; then
  echo "foehncraft-flow (no default features) depends on C code:"
  printf '%s\n' "$bad"
  exit 1
fi
echo "foehncraft-flow (no default features) is C-free: $(printf '%s\n' "$tree" | wc -l) crates"
