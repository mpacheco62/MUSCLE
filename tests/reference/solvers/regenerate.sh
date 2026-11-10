#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Regenerates the reference values of the return mapping solvers that keep a reference
# (tests/reference/solvers/{closest_point_perzyna,closest_point_time}_dp.txt; reduced to a few cases in
# plan R, R6: closest_point, closest_point_2 and cutting_plane have none) by running the kernel reference
# tests tests/materials/test_muscle_material_<solver>.fypp of the declared materials (test registry) with
# `--write FILE` (only the cases with a case number write records). The comment header of an existing file (its
# capture and the approved updates) is kept, and so is every record that the new one matches for the
# reference comparison (NaN = NaN; tests/reference/merge_reference.py); the other records are rewritten.
#
# ONLY run it on code whose behaviour IS the reference: these files were captured from the original
# (pre-FYPP, OO) solvers (HEAD c27dba7, Phase 7 step 1, D-018; removed in Phase 7 step 9, D-018.18),
# composed with the model families of src/ (reference-tested themselves), and the kernels of the
# materials reproduce them bit for bit; the tests exist to detect a change beyond the tolerance of their category (D-021).
# Regenerating after a change would hide it. The tests run by ctest never write these files. Only the
# records changed by an approved decision (e.g. D-018.3) may be updated afterwards, with a line in the
# header.
#
# A NEW material of the test registry needs a case in its test (case number 0: no records); the reference
# records are not regenerated for it.
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

for file in closest_point_perzyna closest_point_time; do
    exe="$build_dir/tests/materials/test_mat_${file}"
    out="$ref_dir/${file}_dp.txt"
    "$ref_dir/../write_reference.sh" "$exe" "$out" "commit $commit$dirty, $version, CMAKE_BUILD_TYPE=$build_type" \
        "tests/reference/solvers/regenerate.sh (test_mat_${file} --write)"
done
