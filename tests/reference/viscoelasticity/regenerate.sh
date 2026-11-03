#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the golden values of the viscoelasticity family
# (tests/reference/viscoelasticity/viscoelasticity_dp.txt) by running the golden test
# tests/viscoelasticity/test_visco_golden.fypp (kernel API, D-018.18f) with `--write FILE`.
# The comment header of the existing file (its capture, approved updates and removals) is kept, and so
# is every record that the new one matches for the golden comparison (NaN = NaN;
# tests/reference/merge_golden.py); the other records are rewritten.
#
# ONLY run it on code whose behaviour IS the reference: this file was captured from the original
# (pre-FYPP, OO) modules src/viscoelasticity/muscle_visco_{state,frozen,base,maxwell,kelvin}.f90
# (with the defaults of src/elasticity/muscle_elasticity_base.f90 that Elasticity_frozen inherited),
# with the approved changes recorded in the "Updated" lines of the header (D-017.13a/b); the kernel
# API reproduces it bit for bit, and the test exists to detect any change against them. Regenerating
# after a change would hide it. The test run by ctest never writes this file.
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
"$ref_dir/../write_golden.sh" "$exe" "$out" "commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type" \
    "tests/reference/viscoelasticity/regenerate.sh (test_visco_golden --write)"
