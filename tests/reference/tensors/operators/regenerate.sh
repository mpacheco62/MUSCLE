#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the golden values of the operators between tensor types
# (tests/reference/tensors/operators/ops_<module>_dp.txt) by running the golden tests of
# tests/tensors/operators/test_ops_golden_<module>.fypp with `--write FILE`.
#
# ONLY run it on code whose behaviour IS the reference: these files were captured from the
# original (pre-FYPP) operator modules, and the tests exist to detect any change against them.
# Regenerating after a change would hide it. The test run by ctest never writes these files.
#
# Usage: tests/reference/tensors/operators/regenerate.sh BUILD_DIR
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
if ! git -C "$repo" diff --quiet HEAD -- src/tensors; then
    dirty=" (with uncommitted changes in src/tensors)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)

for module in addition_subtraction assignment contraction_double contraction_single dyadic transform; do
    exe="$build_dir/tests/tensors/operators/test_ops_golden_$module"
    out="$ref_dir/ops_${module}_dp.txt"
    tmp=$(mktemp)
    "$exe" --write "$tmp"
    {
        echo "# Captured from commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type,"
        echo "# by tests/reference/tensors/operators/regenerate.sh (test_ops_golden_$module --write)."
        cat "$tmp"
    } > "$out"
    rm -f "$tmp"
    echo "wrote $out ($(grep -vc '^#' "$out") records)"
done
