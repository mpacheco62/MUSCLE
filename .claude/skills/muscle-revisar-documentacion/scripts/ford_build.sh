#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Finds FORD (and installs it with --install), builds the documentation of MUSCLE and summarizes the warnings.
#
#   ford_build.sh <repo> [--install] [--warn] [OUTPUT_DIR]       (OUTPUT_DIR default <repo>/build-check/doc)
#
# FORD is looked up, in order: in the PATH; in a conda environment named `ford` (the one environment.yml creates);
# with --install, by creating that environment from environment.yml (conda/mamba), or else a Python venv in
# <repo>/build-check/.ford-venv with `pip install ford fypp` (needs the network). Installing is never automatic.
# The build runs from the repository root (my_project.md preprocesses the templates with fypp). It must end without
# errors; with --warn FORD also lists undocumented entities: the summary groups them and leaves out the known noise
# (local variables, the `module procedure` entries of a generic, the `public` of re-exported names).
# tools/check.sh looks FORD up in the same places.
set -u
repo=$(readlink -f "$1"); shift
install=0; warn=0; out=""
for a in "$@"; do
    case "$a" in
        --install) install=1 ;;
        --warn) warn=1 ;;
        *) out=$a ;;
    esac
done
out=${out:-$repo/build-check/doc}
bin=""
if command -v ford > /dev/null; then
    bin=$(dirname "$(command -v ford)")
else
    for c in conda mamba "$HOME/miniforge3/condabin/conda" "$HOME/miniconda3/condabin/conda" "$HOME/anaconda3/condabin/conda"; do
        if command -v "$c" > /dev/null 2>&1; then
            prefix=$("$c" env list 2>/dev/null | awk '$1 == "ford" {print $NF}')
            [[ -n "$prefix" && -x "$prefix/bin/ford" ]] && { bin="$prefix/bin"; break; }
            if [[ -z "$prefix" && $install -eq 1 ]]; then
                echo "creating the conda environment ford from environment.yml (with $c)"
                "$c" env create -f "$repo/environment.yml" || exit 1
                prefix=$("$c" env list | awk '$1 == "ford" {print $NF}'); bin="$prefix/bin"; break
            fi
        fi
    done
    if [[ -z "$bin" && -x "$repo/build-check/.ford-venv/bin/ford" ]]; then
        bin="$repo/build-check/.ford-venv/bin"
    elif [[ -z "$bin" && $install -eq 1 ]]; then
        echo "creating a venv with FORD and fypp in build-check/.ford-venv"
        python3 -m venv "$repo/build-check/.ford-venv" && "$repo/build-check/.ford-venv/bin/pip" install -q ford 'fypp>=3.1' \
            || exit 1
        bin="$repo/build-check/.ford-venv/bin"
    fi
fi
if [[ -z "$bin" ]]; then
    echo "FORD not found (PATH, conda environment ford, build-check/.ford-venv): run again with --install"
    exit 2
fi
echo "FORD: $("$bin/ford" --version 2>&1 | head -1) ($bin)"
log="$out.log"; mkdir -p "$(dirname "$out")"
flags=(--output_dir "$out"); (( warn )) && flags+=(--warn)
start=$SECONDS
(cd "$repo" && PATH="$bin:$PATH" PYTHONDONTWRITEBYTECODE=1 "$bin/ford" "${flags[@]}" my_project.md) > "$log" 2>&1
rc=$?
echo "build: exit $rc in $((SECONDS - start)) s, log $log, pages $out/index.html"
grep -iE '^(error|critical)|Traceback|Error:' "$log" | head -20
if (( warn )); then
    echo "warnings (without the known noise):"
    grep -E '^Warning' "$log" | grep -vE "Undocumented (variable|moduleprocedure)|Unknown entity '[^']*' with attribute 'public'" \
        | sed -E "s/'[^']*'/'...'/; s/ in file .*//; s/ in module .*//" | sort | uniq -c | sort -rn | head -20
fi
exit $rc
