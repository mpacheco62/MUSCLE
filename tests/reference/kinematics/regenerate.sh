#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the golden values of the kinematics family
# (tests/reference/kinematics/kinematics_dp.txt) by running the golden test
# tests/kinematics/test_kin_golden.fypp (kernel API, D-018.18f) with `--write FILE`
# (tests/reference/write_golden.sh: the header of the existing file and its matching records are kept).
#
# ONLY run it on code whose behaviour IS the reference: this file was captured from the original
# (pre-FYPP, OO) modules src/kinematics/*.f90 (HEAD 62b7c18), with the approved updates of the header;
# the kernel API reproduces it bit for bit, and the test exists to detect a change beyond the tolerance of its category (D-021).
# Regenerating after a change would hide it. The test run by ctest never writes this file.
#
# Usage: tests/reference/kinematics/regenerate.sh BUILD_DIR
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
if ! git -C "$repo" diff --quiet HEAD -- src/kinematics; then
    dirty=" (with uncommitted changes in src/kinematics)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)

exe="$build_dir/tests/kinematics/test_kin_golden"
out="$ref_dir/kinematics_dp.txt"
"$ref_dir/../write_golden.sh" "$exe" "$out" "commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type" \
    "tests/reference/kinematics/regenerate.sh (test_kin_golden --write)"
