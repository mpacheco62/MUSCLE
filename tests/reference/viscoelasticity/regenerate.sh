#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the golden values of the viscoelasticity family
# (tests/reference/viscoelasticity/viscoelasticity_dp.txt) by running the golden test
# tests/viscoelasticity/test_visco_golden.f90 with `--write FILE`.
#
# ONLY run it on code whose behaviour IS the reference: this file was captured from the original
# (pre-FYPP, OO) modules src/viscoelasticity/muscle_visco_{state,frozen,base,maxwell,kelvin}.f90
# (with the defaults of src/elasticity/muscle_elasticity_base.f90 that Elasticity_frozen inherits),
# and the test exists to detect any change against them. Regenerating after a change would hide it.
# The test run by ctest never writes this file.
#
# Approved changes since the capture (D-017.13a/b, recorded in the "Updated" lines of the file header):
# the current reference is what the original code gives with the gamma series threshold 1D-4 (instead
# of 1D-5) in Generalized_maxwell%coefficients and without the dt < 0 cases (load history 8 and the
# dt = -0.01 entry of the Maxwell coefficient tables), and it is also exactly what the current code
# (FYPP templates, kernel API with status) writes. Regenerating rewrites the header: the "Updated"
# lines are lost and must be copied back by hand.
#
# Usage: tests/reference/viscoelasticity/regenerate.sh BUILD_DIR
#   BUILD_DIR: a Debug build (no -ffast-math) with `cmake --build BUILD_DIR --target build_all_tests` done.
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 BUILD_DIR" >&2
    exit 2
fi
build_dir=$(cd "$1" && pwd)
ref_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$ref_dir/../../.." && pwd)
commit=$(git -C "$repo" rev-parse --short HEAD)
dirty=""
if ! git -C "$repo" diff --quiet HEAD -- src/viscoelasticity src/elasticity; then
    dirty=" (with uncommitted changes in src/viscoelasticity or src/elasticity)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)

exe="$build_dir/tests/viscoelasticity/test_visco_golden"
out="$ref_dir/viscoelasticity_dp.txt"
tmp=$(mktemp)
"$exe" --write "$tmp"
{
    echo "# Captured from commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type,"
    echo "# by tests/reference/viscoelasticity/regenerate.sh (test_visco_golden --write)."
    cat "$tmp"
} > "$out"
rm -f "$tmp"
echo "wrote $out ($(grep -vc '^#' "$out") records)"
