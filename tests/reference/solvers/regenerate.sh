#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the golden values of the five return mapping solvers
# (tests/reference/solvers/{closest_point,closest_point_perzyna,closest_point_time,closest_point_2,
# cutting_plane}_dp.txt) by running the kernel golden tests tests/materials/test_mat_<solver>.fypp of the
# declared materials (test registry) with `--write FILE`. The comment header of an existing file (its
# capture and the approved updates) is kept, and so is every record that the new one matches for the
# golden comparison (NaN = NaN; tests/reference/merge_golden.py); the other records are rewritten.
#
# ONLY run it on code whose behaviour IS the reference: these files were captured from the original
# (pre-FYPP, OO) solvers (HEAD c27dba7, Phase 7 step 1, D-018; removed in Phase 7 step 9, D-018.18),
# composed with the model families of src/ (golden-tested themselves), and the kernels of the
# materials reproduce them bit for bit; the tests exist to detect any change against them.
# Regenerating after a change would hide it. The tests run by ctest never write these files. Only the
# records changed by an approved decision (e.g. D-018.3) may be updated afterwards, with a line in the
# header.
#
# Usage: tests/reference/solvers/regenerate.sh BUILD_DIR
#   BUILD_DIR: a Debug build (no -ffast-math) of the test registry of materials (no
#   -DMUSCLE_MATERIALS) with `cmake --build BUILD_DIR --target build_all_tests` done.
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
if ! git -C "$repo" diff --quiet HEAD -- src; then
    dirty=" (with uncommitted changes in src)"
fi
compiler=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$build_dir/CMakeCache.txt" | cut -d= -f2)
version=$("$compiler" --version | head -n1)
build_type=$(grep -m1 '^CMAKE_BUILD_TYPE:' "$build_dir/CMakeCache.txt" | cut -d= -f2)

for file in closest_point closest_point_perzyna closest_point_time closest_point_2 cutting_plane; do
    exe="$build_dir/tests/materials/test_mat_${file}"
    out="$ref_dir/${file}_dp.txt"
    "$ref_dir/../write_golden.sh" "$exe" "$out" "commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type" \
        "tests/reference/solvers/regenerate.sh (test_mat_${file} --write)"
done
