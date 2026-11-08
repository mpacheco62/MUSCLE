#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the golden values of the elasticity family
# (tests/reference/elasticity/elasticity_dp.txt) by running the golden test
# tests/elasticity/test_elas_golden.f90 (kernel API, D-018.18f) with `--write FILE`
# (tests/reference/write_golden.sh: the header of the existing file and its matching records are kept).
#
# ONLY run it on code whose behaviour IS the reference: this file was captured from the original
# (pre-FYPP, OO) modules src/elasticity/muscle_elasticity_{base,linear}.f90, the kernel API reproduces
# it bit for bit, and the test exists to detect a change beyond the tolerance of its category (D-021). Regenerating after a change
# would hide it. The test run by ctest never writes this file.
#
# Usage: tests/reference/elasticity/regenerate.sh BUILD_DIR
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
if ! git -C "$repo" diff --quiet HEAD -- src/elasticity; then
    dirty=" (with uncommitted changes in src/elasticity)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)

exe="$build_dir/tests/elasticity/test_elas_golden"
out="$ref_dir/elasticity_dp.txt"
"$ref_dir/../write_golden.sh" "$exe" "$out" "commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type" \
    "tests/reference/elasticity/regenerate.sh (test_elas_golden --write)"
