#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the minimal reference of the math modules (tests/reference/math/math_reference_dp.txt, plan R)
# by running tests/math/test_muscle_math_reference.f90 with `--write FILE` (tests/reference/write_reference.sh: the
# header of the existing file and its matching records are kept).
#
# ONLY run it on code whose behaviour IS the reference: the records were captured from the original
# (pre-FYPP) modules src/math/muscle_math_{operations,spectral_derivs}.f90, and the test exists to detect a
# change beyond the tolerance of its category (D-021). Regenerating after a change would hide it. The test
# run by ctest never writes this file.
#
# Usage: tests/reference/math/regenerate.sh BUILD_DIR
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
if ! git -C "$repo" diff --quiet HEAD -- src/math; then
    dirty=" (with uncommitted changes in src/math)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
lapack=$(grep -m1 '^LAPACK_lapack_LIBRARY:' "$build_dir/CMakeCache.txt" | cut -d= -f2 || true)
if [[ -n "$lapack" ]]; then
    lapack=", LAPACK $(readlink -f "$lapack")"
fi

"$ref_dir/../write_reference.sh" "$build_dir/tests/math/test_math_reference" \
    "$ref_dir/math_reference_dp.txt" \
    "commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type$lapack" \
    "tests/reference/math/regenerate.sh (test_math_reference --write)"
