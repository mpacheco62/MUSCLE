#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the golden values of the yield criteria
# (tests/reference/plasticity/yield_criteria/yield_criteria_dp.txt) by running the golden test
# tests/plasticity/yield_criteria/test_yield_golden.f90 with `--write FILE`.
#
# ONLY run it on code whose behaviour IS the reference: this file was captured from the original
# (pre-FYPP, OO) modules src/plasticity/yield_criteria/muscle_yield_*.f90 at commit 3d9a230
# (whose base numeric derivatives already use the D-016.4 finite-difference step at a zero
# argument), and the test exists to detect any change against them. Regenerating after a change
# would hide it. The test run by ctest never writes this file.
#
# Usage: tests/reference/plasticity/yield_criteria/regenerate.sh BUILD_DIR
#   BUILD_DIR: a Debug build (no -ffast-math) with `cmake --build BUILD_DIR --target build_all_tests` done.
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 BUILD_DIR" >&2
    exit 2
fi
build_dir=$(cd "$1" && pwd)
ref_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$ref_dir/../../../.." && pwd)
commit=$(git -C "$repo" rev-parse --short HEAD)
dirty=""
if ! git -C "$repo" diff --quiet HEAD -- src/plasticity/yield_criteria; then
    dirty=" (with uncommitted changes in src/plasticity/yield_criteria)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)

exe="$build_dir/tests/plasticity/yield_criteria/test_yield_golden"
out="$ref_dir/yield_criteria_dp.txt"
tmp=$(mktemp)
"$exe" --write "$tmp"
{
    echo "# Captured from commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type,"
    echo "# by tests/reference/plasticity/yield_criteria/regenerate.sh (test_yield_golden --write)."
    cat "$tmp"
} > "$out"
rm -f "$tmp"
echo "wrote $out ($(grep -vc '^#' "$out") records)"
