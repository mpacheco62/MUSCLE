#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
# Lists the ids of the mutations of one source file in a mutation area of MUSCLE, ready for --only.
#
# Usage: mutation_ids.sh REPO AREA FILE
#   REPO: root of the MUSCLE repository; AREA: a tests/mutations/<area>.toml; FILE: path relative to REPO.
# The hand-written [[manual]] mutations of an area (ids manualN) are not listed: add them by hand.
# Prints one comma-separated line; run the check with
#   /usr/bin/python3 tools/mutation_check.py AREA --samples all --only "<that line>"
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo "usage: $0 REPO AREA FILE" >&2
    exit 2
fi
repo=$1 area=$2 file=$3
cd "$repo"
# tools/mutation_check.py needs tomllib (Python >= 3.11): prefer the system interpreter
py=/usr/bin/python3
[[ -x $py ]] || py=python3
# (grep finds nothing for a file outside the area: `|| true`, or set -o pipefail ends the script here)
ids=$("$py" tools/mutation_check.py "$area" --samples all --list | { grep -E "^[a-z]+:${file}:" || true; } | paste -sd, -)
if [[ -z $ids ]]; then
    echo "no mutation of $file in area $area (is the file in the sources of tests/mutations/$area.toml?)" >&2
    exit 1
fi
echo "$ids"
