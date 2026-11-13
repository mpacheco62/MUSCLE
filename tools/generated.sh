#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Prints the Fortran file(s) that fypp generates from a template, to read or open them in an editor
# (language servers such as fortls understand the generated code, not the templates).
#
#   tools/generated.sh TEMPLATE [BUILD_DIR]
#
# BUILD_DIR defaults to build-check. A template of src/ or tests/ gives BUILD_DIR/fypp/<same path>.f90;
# the material template and its includes (src/materials/muscle_material*.fypp) give the module of every
# declared material, BUILD_DIR/fypp/src/materials/muscle_mat_<name>.f90.

set -u
cd "$(dirname "$0")/.."
[[ $# -ge 1 ]] || { echo "usage: $0 TEMPLATE [BUILD_DIR]" >&2; exit 2; }
tpl="${1#./}"
build_dir="${2:-build-check}"
case "$tpl" in
    src/materials/muscle_material*.fypp|src/solvers/*.fypp)
        ls "$build_dir"/fypp/src/materials/muscle_mat_*.f90 2> /dev/null || {
            echo "no generated materials in $build_dir (build it first)" >&2; exit 1; } ;;
    *.fypp)
        out="$build_dir/fypp/${tpl%.fypp}.f90"
        if [[ ! -f "$out" ]]; then
            if ! grep -qiE '^[[:space:]]*(module|program)[[:space:]]' "$tpl"; then
                echo "$tpl has no output of its own: it is included by other templates" >&2
            else
                echo "$out does not exist (build $build_dir first)" >&2
            fi
            exit 1
        fi
        echo "$out" ;;
    *) echo "$tpl is not a .fypp template" >&2; exit 2 ;;
esac
