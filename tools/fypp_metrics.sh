#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Size of the FYPP layer (Phase 9, D-020): shows the progress of the readability refactor.
#   tools/fypp_metrics.sh [--list]
# --list also prints, per library template, its macro calls and lines.

set -u
cd "$(dirname "$0")/.."

# Macros that build declarations from metadata or logic (forbidden by D-020 in new code).
forbidden='MODEL_FUNCTION|MODEL_PARAMS_TYPE|MODEL_INIT|MODEL_PUBLIC|SOLVE_DUMMIES|FAMILY_USES|FAMILY_EXPORTS|YIELD_ZERO_AT_NULL'

lines() { cat "$@" 2> /dev/null | wc -l; }
count() { grep -hoE "$1" "${@:2}" 2> /dev/null | wc -l; }

lib=$(find src -name '*.fypp' | sort)
tst=$(find tests -name '*.fypp' | sort)
inc=$(find include -name '*.fypp' | sort)
printf '%-44s %7s\n' 'library templates (src/*.fypp)' "$(echo "$lib" | wc -l)"
printf '%-44s %7s\n' '  lines' "$(lines $lib)"
printf '%-44s %7s\n' '  with model/solver macros (forbidden list)' \
    "$(grep -lE "[\$@]:($forbidden)\b|[\$@]:[A-Z_]+_BODY" $lib | wc -l)"
printf '%-44s %7s\n' '  calls to the forbidden-list macros' "$(count "[\$@]:($forbidden)\b" $lib)"
printf '%-44s %7s\n' '  calls to @:DEVICE (kept)' "$(count '@:DEVICE\b' $lib)"
printf '%-44s %7s\n' 'include/*.fypp lines' "$(lines $inc)"
printf '%-44s %7s\n' '  #:def macros' "$(count '^\s*#:def ' $inc)"
printf '%-44s %7s\n' 'tools/fypp/muscle_meta.py lines' "$(lines tools/fypp/muscle_meta.py)"
printf '%-44s %7s\n' '  functions' "$(count '^def ' tools/fypp/muscle_meta.py)"
printf '%-44s %7s\n' 'test templates (tests/*.fypp)' "$(echo "$tst" | wc -l)"
printf '%-44s %7s\n' '  lines' "$(lines $tst)"

if [[ "${1:-}" == --list ]]; then
    echo; printf '%-70s %6s %6s\n' template macros lines
    for f in $lib; do
        printf '%-70s %6s %6s\n' "$f" "$(count '[$@]:[A-Z_]+' "$f")" "$(wc -l < "$f")"
    done
fi
