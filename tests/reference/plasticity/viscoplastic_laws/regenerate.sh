#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the golden values of the viscoplastic laws
# (tests/reference/plasticity/viscoplastic_laws/vp_laws_dp.txt) by running the golden test
# tests/plasticity/viscoplastic_laws/test_vp_golden.f90 with `--write FILE`.
#
# ONLY run it on code whose behaviour IS the reference: this file was captured from the original
# (pre-FYPP, OO) modules src/plasticity/viscoplastic_laws/muscle_vp_*.f90 at commit 5ab2802
# (which already includes the D-016.4 finite-difference step at a zero argument), composed with
# the hardening laws of src/plasticity/hardening_laws (golden-tested themselves), and the test
# exists to detect any change against them. Regenerating after a change would hide it. The test run by ctest never writes this file.
#
# Usage: tests/reference/plasticity/viscoplastic_laws/regenerate.sh BUILD_DIR
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
if ! git -C "$repo" diff --quiet HEAD -- src/plasticity/viscoplastic_laws; then
    dirty=" (with uncommitted changes in src/plasticity/viscoplastic_laws)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)

exe="$build_dir/tests/plasticity/viscoplastic_laws/test_vp_golden"
out="$ref_dir/vp_laws_dp.txt"
tmp=$(mktemp)
"$exe" --write "$tmp"
{
    echo "# Captured from commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type,"
    echo "# by tests/reference/plasticity/viscoplastic_laws/regenerate.sh (test_vp_golden --write)."
    cat "$tmp"
} > "$out"
rm -f "$tmp"
echo "wrote $out ($(grep -vc '^#' "$out") records)"
