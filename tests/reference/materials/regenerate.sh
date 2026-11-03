#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the dp regression golden of the finite strain materials
# (tests/reference/materials/mat_kinematics_dp.txt, D-018.15) by running the test
# tests/materials/test_mat_kinematics.fypp with `--write FILE`.
#
# There is no historical composition of a return mapping with a kinematics: the file is a regression
# reference of the kernel itself (checked against physical oracles by the same test). ONLY run it on
# code whose behaviour IS the reference (a reviewed change of these materials); regenerating after an
# unreviewed change would hide it. The test run by ctest never writes this file.
#
# Usage: tests/reference/materials/regenerate.sh BUILD_DIR
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
if ! git -C "$repo" diff --quiet HEAD -- src tools/fypp; then
    dirty=" (with uncommitted changes in src or tools/fypp)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)

exe="$build_dir/tests/materials/test_mat_kinematics"
out="$ref_dir/mat_kinematics_dp.txt"
tmp=$(mktemp)
"$exe" --write "$tmp"
{
    echo "# Captured from commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type,"
    echo "# by tests/reference/materials/regenerate.sh (test_mat_kinematics --write)."
    cat "$tmp"
} > "$out"
rm -f "$tmp"
echo "wrote $out ($(grep -vc '^#' "$out") records)"
