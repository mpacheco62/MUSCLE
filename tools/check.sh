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
#    `build_all_tests`, check the operator signatures (tools/fypp/ops_signatures.py) and run
#    ctest; build the FORD documentation if FORD is installed. Skipped with --lint-only.
#
# Opt-outs:
# - `! lint: allow(<rule>)` at the end of a line skips that rule on that line.
# - `! lint: host-only` anywhere in a template exempts the whole file from the GPU kernel
#   rules device, alloc and io: use it for host code such as drivers (LAPACK, UDTIO). The rules
#   oop and procarg (no dynamic dispatch, no procedure arguments, D-018.18b) apply to every library
#   source (src/, include/), host-only or not; `class(...)` only in the type-bound wrappers of
#   the tensor types (`*_tbp`, historical API kept by D-018.18e).
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
# Every library source (templates and plain Fortran, host-only included): no dynamic dispatch
mapfile -t library_files < <(printf '%s\n' "${fortran_files[@]}" | grep -E '^(src|include)/')

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

# The library is compiled with -cpp, which reads every line, comments included: an odd number of
# apostrophes (or double quotes) in a comment is an unterminated character literal for cpp and
# gives a "missing terminating ' character" warning in every module that includes the line. Rewrite
# the comment without them (no possessive 's, no lone primes).
if (( ${#library_files[@]} > 0 )); then
    report quote "odd number of ' or \" in a comment of a library source (cpp warning)" "$(perl -ne '
        # fypp directives and their continuation lines are Python, not Fortran
        my $directive = /^\s*(#|[\$\@]:)/ || ($cont && /^\s*&/);
        $cont = $directive && /&\s*$/;
        next if $directive || /lint:\s*allow\(quote\)/;
        chomp(my $l = $_);
        $l =~ s/\$\{.*?\}\$/FYPP/g;
        # comment: the text after the first ! outside a character literal
        $l =~ /^((?:\x27[^\x27]*\x27|"[^"]*"|[^\x27"!])*)!(.*)$/ or next;
        my $comment = $2;
        my $sq = () = $comment =~ /\x27/g;
        my $dq = () = $comment =~ /"/g;
        print "$ARGV:$.: $_" if $sq % 2 || $dq % 2;
    } continue { if (eof) { close ARGV; $cont = 0 } }' "${library_files[@]}")"
fi

# Kinds and literals (precision-kinds.md)
rule kind "hard-coded or default real kind, use real(\${k}\$) or real(dp)" \
    '(?i)\breal\s*\(\s*(kind\s*=\s*)?(4|8|real32|real64)\s*\)|\breal\s*\*\s*[48]\b|\bdouble\s+precision\b|^\s*real\s*(::|,)' \
    "${fypp_files[@]}"
rule kind "dble()/sngl(), use real(x, \${k}\$)" '(?i)\b(dble|sngl)\s*\(' "${fypp_files[@]}"
# Sources not yet migrated: at least no literal real kinds (Phase 1 removed them). They still
# use real(real64), which is replaced when each module is migrated.
mapfile -t legacy_files < <(listed '*.f90' '*.F90' '*.inc')
rule legacy-kind "literal kind in a non-template source, use real(dp)" \
    '(?i)\breal\s*\(\s*(kind\s*=\s*)?[48]\s*\)|\breal\s*\*\s*[48]\b' \
    "${legacy_files[@]}"
rule literal "D-exponent literal, use 1.0e-5_\${k}\$" '(?i)(?<![\w.])\d+\.?\d*d[+-]?\d+\b' "${fypp_files[@]}"
rule literal "real literal without kind (default real is single precision)" \
    '(?i)(?<![\w.])(\d+\.\d*|\.\d+)(e[+-]?\d+)?(?![\w.])|(?<![\w.])\d+e[+-]?\d+(?![\w.])' \
    "${fypp_files[@]}"

# Static dispatch in the whole library (gpu-static-dispatch.md, D-018.18b): no abstract types,
# deferred bindings, select type or procedure arguments/pointers, host-only code included
rule oop "dynamic dispatch (abstract/deferred/select type)" \
    '(?i)\babstract\b(?!\s+interface)|\bdeferred\b|\bselect\s+type\b' "${library_files[@]}"
rule procarg "procedure arguments/pointers" '(?i)\bprocedure\s*\(' "${library_files[@]}"
# A `class(...)` dummy only in the type-bound wrappers of the tensor types (`*_tbp`, D-011,
# D-018.18e): src/tensors/ and their common macros (include/muscle_tensor_common.fypp)
if (( ${#library_files[@]} > 0 )); then
    report oop "class(...) outside a type-bound wrapper of a tensor type (*_tbp)" "$(perl -ne '
        if (/^\s*(?:pure\s+|elemental\s+|recursive\s+)*(?:function|subroutine)\s+([\w\$\{\}.]+)/i) { $p = $1 }
        $p = "" if /^\s*end\s+(function|subroutine)\b/i;
        $tbp = $p =~ /_tbp$/ && $ARGV =~ m{^(src/tensors/|include/muscle_tensor_)};
        print "$ARGV:$.: $_" if /\bclass\s*\(/i && !/^\s*!/ && !$tbp && !/lint:\s*allow\(oop\)/;
    } continue { close ARGV if eof' "${library_files[@]}")"
fi

# Per-procedure kernel rule (D-011): every pure/elemental procedure marked @:DEVICE().
if (( ${#kernel_files[@]} > 0 )); then
    report device "pure/elemental kernel procedure without @:DEVICE()" "$(perl -ne '
        if (/^\s*((?:pure|elemental)\s+(?:\w+\s+)*?(?:function|subroutine))\s+([\w\$\{\}.]+)/i) {
            ($name, $line, $dev) = ($2, $., 0);
        }
        $dev = 1 if /^\s*\@:DEVICE\(\)/;
        if (/^\s*end\s+(function|subroutine)\b/i && $name) {
            print "$ARGV:$line: $name\n" unless $dev || $name =~ /_tbp$/;
            $name = "";
        }
    } continue { close ARGV if eof' "${kernel_files[@]}")"
fi
rule alloc "allocatable/pointer" '(?i)\b(allocatable|pointer)\b' "${kernel_files[@]}"
rule io "I/O or stop in kernel code" \
    '(?i)^\s*(\d+\s+)?(if\s*\(.*\)\s*)?(print\b|write\s*\(|read\s*\(|stop\b|error\s+stop\b)' \
    "${kernel_files[@]}"

if (( failures > 0 )); then
    echo "== lint: ${failures} rule(s) failed"
    exit 1
fi
echo "== lint: OK"
(( lint_only )) && exit 0

# Module files of deleted or renamed modules stay in an incremental build directory and would let
# a stale `use` compile: remove every library .mod whose module no source defines any more.
prune_stale_modules() {
    local b="$1" dir="$1/modules/muscle" mod name f
    [[ -d "$dir" ]] || return 0
    local sources=()
    mapfile -t sources < <(listed 'src/*.f90' 'src/*.F90')
    for f in $(listed 'src/*.fypp'); do
        [[ -f "$b/fypp/${f%.fypp}.f90" ]] && sources+=("$b/fypp/${f%.fypp}.f90")
    done
    for mod in "$dir"/*.mod; do
        [[ -e "$mod" ]] || continue
        name=$(basename "$mod" .mod)
        if ! grep -qiE "^[[:space:]]*module[[:space:]]+${name}[[:space:]]*(!.*)?$" "${sources[@]}"; then
            echo "   removing stale module file $mod"
            rm -f "$mod" "$mod.stamp"
        fi
    done
}
prune_stale_modules "$build_dir"
prune_stale_modules "${build_dir}-udtio"

echo "== build: ${build_dir} (Debug, MUSCLE_REAL_KINDS=sp;dp)"
cmake -S . -B "$build_dir" -G Ninja -DCMAKE_BUILD_TYPE=Debug "-DMUSCLE_REAL_KINDS=sp;dp" \
    > /dev/null || exit 1
cmake --build "$build_dir" --target build_all_tests -j "$(nproc)" || exit 1

# Generated code must fit the free-form limit (comments may run longer)
long=$(find "$build_dir/fypp" -name '*.f90' -print0 2>/dev/null | xargs -0 -r perl -ne '
    next if /^#/;                      # #line markers (paths can be long in worktrees)
    my $code = $_; $code =~ s/!.*//;
    print "$ARGV:$.: ", length($code), " columns\n" if length($code) > 133;
    close ARGV if eof')
if [[ -n "$long" ]]; then
    echo "LINT [line-length] generated code longer than 132 columns:"; echo "$long" | sed 's/^/    /'
    exit 1
fi

# Literals written by the Python helpers (tools/fypp/muscle_meta.py) only appear in the
# generated code: apply the literal-kind rule to it too
mapfile -t generated < <(for f in "${fypp_files[@]}"; do
    [[ "$f" == src/* && -f "$build_dir/fypp/${f%.fypp}.f90" ]] && echo "$build_dir/fypp/${f%.fypp}.f90"
done)
rule literal "real literal without kind in generated code" \
    '(?i)(?<![\w.])(\d+\.\d*|\.\d+)(e[+-]?\d+)?(?![\w.])|(?<![\w.])\d+e[+-]?\d+(?![\w.])' \
    "${generated[@]}"
(( failures > 0 )) && exit 1

# The generic operators keep every historical double precision signature (D-015). One source
# per operator module (the template output if it exists): stale outputs of renamed templates
# left in the build directory must not hide a missing signature.
ops=()
for m in $(PYTHONDONTWRITEBYTECODE=1 python3 -c 'import sys; sys.path.insert(0, "tools/fypp")
import muscle_meta; print(*muscle_meta.OP_MODULES)'); do
    src="src/tensors/operators/muscle_tensor_ops_$m"
    if [[ -f "$src.fypp" ]]; then ops+=("$build_dir/fypp/$src.f90"); else ops+=("$src.f90"); fi
done
python3 tools/fypp/ops_signatures.py --check tools/fypp/ref/ops_signatures_dp.txt "${ops[@]}" || {
    echo "LINT [ops-signatures] operator signatures missing (tools/fypp/ref/ops_signatures_dp.txt)"
    exit 1
}

# Formatted output (UDTIO) is disabled for gfortran < 16.3: compile that path anyway
# (library and tests; only the tensor tests run there: the OpenMP tests hit the gfortran bug)
echo "== build: ${build_dir}-udtio (MUSCLE_UDTIO=ON, tensor tests)"
cmake -S . -B "${build_dir}-udtio" -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMUSCLE_UDTIO=ON > /dev/null || exit 1
cmake --build "${build_dir}-udtio" --target build_all_tests -j "$(nproc)" > /dev/null || {
    echo "== build with MUSCLE_UDTIO=ON failed"; exit 1; }
ctest --test-dir "${build_dir}-udtio" -R 'tensors-types-.*_kinds' --output-on-failure -j "$(nproc)" \
    > "${build_dir}-udtio/ctest-udtio.log" || {
    echo "== UDTIO tensor tests failed (see ${build_dir}-udtio/ctest-udtio.log)"; exit 1; }

# Documentation: FORD preprocesses the templates with fypp (my_project.md). Optional: skipped if
# FORD is not installed; it must run without errors (warnings are listed in the log).
if command -v ford > /dev/null; then
    echo "== doc: FORD (${build_dir}/doc)"
    PYTHONDONTWRITEBYTECODE=1 ford --quiet --output_dir "${build_dir}/doc" my_project.md \
        > "${build_dir}/ford.log" 2>&1 || {
        echo "== FORD failed (see ${build_dir}/ford.log)"; tail -20 "${build_dir}/ford.log"; exit 1; }
else
    echo "== doc: FORD not installed, skipped"
fi

ctest --test-dir "$build_dir" --output-on-failure -j "$(nproc)"
