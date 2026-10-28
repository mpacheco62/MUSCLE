#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Verification gate of the FYPP migration (docs/fypp/PLAN.md).
#
#   tools/check.sh [--lint-only] [BUILD_DIR]
#
# 1. Lint: mechanical rules from .claude/rules/ on the migrated sources (*.fypp) and tab
#    checks on every Fortran source.
# 2. Configure BUILD_DIR (default: build-check) as Debug with every real kind, build
#    `build_all_tests` and run ctest. Skipped with --lint-only.
#
# Opt-outs:
# - `! lint: allow(<rule>)` at the end of a line skips that rule on that line.
# - `! lint: host-only` anywhere in a template exempts the whole file from the GPU kernel
#   rules (oop, alloc, procarg, io): use it for host code such as dispatchers or drivers.
# Exit status != 0 on any failure.

set -u
cd "$(dirname "$0")/.."

lint_only=0
if [[ "${1:-}" == "--lint-only" ]]; then
    lint_only=1
    shift
fi
build_dir="${1:-build-check}"

failures=0

# report <rule> <description> <hits>
report() {
    if [[ -n "$3" ]]; then
        echo "LINT [$1] $2:"
        echo "$3" | sed 's/^/    /'
        failures=$((failures + 1))
    fi
}

# code_lines <files...>: "file:line:code" with comment lines, fypp directives, inline fypp
# expressions, strings and trailing comments removed. Lines carrying
# `! lint: allow(<rule>)` end with an `#allow(<rule>)` marker.
code_lines() {
    perl -ne '
        next if /^\s*(!|#)/;
        chomp;
        my $allow = /!\s*lint:\s*allow\((\w+)\)/ ? " #allow($1)" : "";
        s/\$\{.*?\}\$/FYPP/g;
        s/\x27[^\x27]*\x27|"[^"]*"/""/g;
        s/!.*//;
        print "$ARGV:$.:$_$allow\n";
    } continue { close ARGV if eof' "$@"
}

# rule <name> <description> <perl regex> <files...>
rule() {
    local name="$1" desc="$2" regex="$3"
    shift 3
    (( $# == 0 )) && return
    local hits
    hits=$(code_lines "$@" | RULE="$name" RE="$regex" perl -ne '
        my ($loc, $code) = /^([^:]+:\d+):(.*)$/ or next;
        next if $code =~ /#allow\($ENV{RULE}\)/;
        $code =~ s/\s*#allow\(\w+\)$//;
        print "$loc: $code\n" if $code =~ /$ENV{RE}/;
    ')
    report "$name" "$desc" "$hits"
}

listed() { git ls-files --cached --others --exclude-standard -- "$@"; }

mapfile -t fypp_files < <(listed '*.fypp')
mapfile -t fortran_files < <(listed '*.f90' '*.F90' '*.inc' '*.fypp')
mapfile -t py_files < <(listed 'tools/fypp/*.py')
# Library templates subject to the GPU kernel rules
mapfile -t kernel_files < <(printf '%s\n' "${fypp_files[@]}" | grep -E '^(src|include)/' \
                            | xargs -r grep -L '! lint: host-only')

echo "== lint: ${#fypp_files[@]} templates (${#kernel_files[@]} kernel), ${#fortran_files[@]} Fortran sources"

# Formatting (every Fortran source)
report tabs "tab characters" "$(grep -nP '\t' "${fortran_files[@]}" 2>/dev/null)"

# Headers of migrated sources
for f in "${fypp_files[@]}" "${py_files[@]}"; do
    head -n 3 "$f" | grep -q 'SPDX-License-Identifier' || report spdx "missing SPDX header" "$f"
done

# One preprocessor per template (fypp-templates.md)
report cpp "cpp directive in a template, use #:if defined(...)" \
    "$(grep -nP '^\s*#\s*(if|ifdef|ifndef|elif|else|endif|define|undef|include)\b' \
       "${fypp_files[@]}" 2>/dev/null)"

# Kinds and literals (precision-kinds.md)
rule kind "hard-coded or default real kind, use real(\${k}\$) or real(dp)" \
    '(?i)\breal\s*\(\s*(kind\s*=\s*)?(4|8|real32|real64)\s*\)|\breal\s*\*\s*[48]\b|\bdouble\s+precision\b|^\s*real\s*(::|,)' \
    "${fypp_files[@]}"
rule kind "dble()/sngl(), use real(x, \${k}\$)" '(?i)\b(dble|sngl)\s*\(' "${fypp_files[@]}"
rule literal "D-exponent literal, use 1.0e-5_\${k}\$" '(?i)(?<![\w.])\d+\.?\d*d[+-]?\d+\b' "${fypp_files[@]}"
rule literal "real literal without kind (default real is single precision)" \
    '(?i)(?<![\w.])(\d+\.\d*|\.\d+)(e[+-]?\d+)?(?![\w.])|(?<![\w.])\d+e[+-]?\d+(?![\w.])' \
    "${fypp_files[@]}"

# GPU kernel rules (gpu-static-dispatch.md)
rule oop "dynamic dispatch (class/abstract/deferred/select type)" \
    '(?i)\bclass\s*\(|\babstract\b(?!\s+interface)|\bdeferred\b|\bselect\s+type\b' "${kernel_files[@]}"
rule alloc "allocatable/pointer" '(?i)\b(allocatable|pointer)\b' "${kernel_files[@]}"
rule procarg "procedure arguments/pointers" '(?i)\bprocedure\s*\(' "${kernel_files[@]}"
rule io "I/O or stop in kernel code" \
    '(?i)^\s*(\d+\s+)?(if\s*\(.*\)\s*)?(print\b|write\s*\(|read\s*\(|stop\b|error\s+stop\b)' \
    "${kernel_files[@]}"

if (( failures > 0 )); then
    echo "== lint: ${failures} rule(s) failed"
    exit 1
fi
echo "== lint: OK"
(( lint_only )) && exit 0

echo "== build: ${build_dir} (Debug, MUSCLE_REAL_KINDS=sp;dp)"
cmake -S . -B "$build_dir" -G Ninja -DCMAKE_BUILD_TYPE=Debug "-DMUSCLE_REAL_KINDS=sp;dp" \
    > /dev/null || exit 1
cmake --build "$build_dir" --target build_all_tests -j "$(nproc)" || exit 1
ctest --test-dir "$build_dir" --output-on-failure -j "$(nproc)"
