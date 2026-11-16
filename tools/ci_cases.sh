#!/usr/bin/env bash
# Runs one set of the case catalogue the way CI does (R5: each case checks
# its own results and exits non-zero on failure).
#
#   tools/ci_cases.sh SET [FOEHN]     run the set; FOEHN defaults to target/release/foehn
#   tools/ci_cases.sh --check         every case file is in exactly one Linux set
#   tools/ci_cases.sh --list SET      print the case files of the set
#
# Sets (times are the Linux CI runner, main run of 2026-11-15):
#   linux-1  mesh, scalar, poisson, the 3D cavities           about 2.2 min
#   linux-2  cylinder Re 20, cavity Re 1000, Poiseuille       about 2.7 min
#   linux-3  Taylor-Green, transient cavity, the rest         about 2.2 min
#   long     cylinder Re 100 (Schaefer-Turek 2D-2, ~4.6 min)  full tier only
#   quick    macOS and Windows: mesh, scalar, poisson and the quick flows
set -euo pipefail

mesh_cases() { ls examples/cases/mesh/*_box.toml examples/cases/mesh/*_2d.toml examples/cases/mesh/poly_partition.toml; }
scalar_cases() { ls examples/cases/scalar/*.toml | grep -v '_gpu\.toml$'; }
poisson_cases() { ls examples/cases/poisson/*.toml | grep -v '_gpu\.toml$'; }

F=examples/cases/flow
flow_linux_1="$F/cavity_re100.toml $F/cavity_re100_energy.toml $F/cavity_re1000_poly3d.toml $F/couette_hex.toml $F/couette_mixed.toml $F/couette_poly.toml $F/couette_tet.toml"
flow_linux_2="$F/cylinder_re20.toml $F/cavity_re1000.toml $F/taylor_green_cfl10.toml $F/poiseuille_tet.toml $F/poiseuille_hex.toml $F/poiseuille_poly.toml"
flow_linux_3="$F/cavity_re100_budget.toml $F/taylor_green_bdf2.toml $F/taylor_green_cn.toml $F/poiseuille_mixed.toml $F/cavity_re100_transient.toml $F/taylor_green_space.toml $F/cavity_re100_poly3d.toml $F/cavity_re100_tet3d.toml"
flow_long="$F/cylinder_re100.toml"
flow_quick="$F/cavity_re100_budget.toml $F/couette_hex.toml $F/couette_mixed.toml $F/couette_poly.toml $F/couette_tet.toml $F/cavity_re100.toml $F/cavity_re100_poly3d.toml $F/poiseuille_hex.toml $F/poiseuille_poly.toml $F/taylor_green_space.toml"

# prints "COMMAND FILE" lines
list() {
  case "$1" in
    linux-1)
      mesh_cases | sed 's/^/mesh /'
      scalar_cases | sed 's/^/scalar /'
      poisson_cases | sed 's/^/poisson /'
      for f in $flow_linux_1; do echo "flow $f"; done ;;
    linux-2) for f in $flow_linux_2; do echo "flow $f"; done ;;
    linux-3) for f in $flow_linux_3; do echo "flow $f"; done ;;
    long)    for f in $flow_long; do echo "flow $f"; done ;;
    quick)
      mesh_cases | sed 's/^/mesh /'
      scalar_cases | sed 's/^/scalar /'
      poisson_cases | sed 's/^/poisson /'
      for f in $flow_quick; do echo "flow $f"; done ;;
    *) echo "unknown set '$1'" >&2; exit 2 ;;
  esac
}

if [ "${1:-}" = "--check" ]; then
  # A new flow case that is in no set would silently never run: refuse it.
  assigned=$(for s in linux-1 linux-2 linux-3 long; do list "$s"; done | awk '$1=="flow"{print $2}' | sort)
  dupes=$(printf '%s\n' "$assigned" | uniq -d)
  missing=$(comm -23 <(ls $F/*.toml | grep -v '_gpu\.toml$' | sort) <(printf '%s\n' "$assigned" | sort -u))
  [ -z "$dupes" ] || { echo "in more than one set: $dupes" >&2; exit 1; }
  [ -z "$missing" ] || { echo "in no set (add to tools/ci_cases.sh): $missing" >&2; exit 1; }
  echo "every flow case is in exactly one Linux set"; exit 0
fi
if [ "${1:-}" = "--list" ]; then list "$2"; exit 0; fi

set_name="${1:?usage: tools/ci_cases.sh SET [FOEHN]}"
foehn="${2:-target/release/foehn}"
status=0
while read -r cmd file; do
  echo "::group::$cmd $file"
  start=$SECONDS
  if ! "$foehn" "$cmd" "$file"; then echo "::error::$cmd $file failed"; status=1; fi
  echo "::endgroup::"
  echo "$cmd $file: $((SECONDS - start)) s"
done < <(list "$set_name")
exit $status
