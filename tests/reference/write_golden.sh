#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Writes one golden reference file for the regenerate.sh scripts: runs the golden test EXE with
# `--write`, and keeps the comment header of an existing OUT (its capture, approved updates and
# removals) and every record that the new one matches for the golden comparison (NaN = NaN,
# merge_golden.py); a new OUT gets the capture header (CAPTURE: commit, compiler, build type; BY: the
# script and test that wrote it).
#
# Usage: write_golden.sh EXE OUT CAPTURE BY
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "usage: $0 EXE OUT CAPTURE BY" >&2
    exit 2
fi
exe=$1
out=$2
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
tmp=$(mktemp)
"$exe" --write "$tmp" > /dev/null
if [[ -f "$out" ]]; then
    python3 -I "$here/merge_golden.py" "$out" "$tmp" > "$out.new"
else
    {
        echo "# Captured from $3,"
        echo "# by $4."
        cat "$tmp"
    } > "$out.new"
fi
mv "$out.new" "$out"
rm -f "$tmp"
echo "wrote $out ($(grep -vc '^#' "$out") records)"
