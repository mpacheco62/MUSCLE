#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Verification gate of MUSCLE (pages/contributing.md, Development cycle).
#
#   tools/check.sh [--lint-only | --doc | --quick] [--gpu-full] [BUILD_DIR]
#
# Modes: --lint-only (stage 1, seconds, while editing); --doc (stage 1, the FORD documentation and the drift between
# the documentation and the code, ~2 min, before a pull request that touches documentation or a public API);
# --quick (stages 1-2 without UDTIO and FORD, ~1-2 min, before every commit); none (every stage, ~25 min, to close a
# stage or before a large pull request).
#
# 1. Lint: mechanical rules from .claude/rules/ on the migrated sources (*.fypp) and tab
#    checks on every Fortran source.
# 2. Configure BUILD_DIR (default: build-check) as Debug with every real kind, build
#    `build_all_tests`, check the operator signatures (tools/fypp/ops_signatures.py) and run
#    ctest; build the FORD documentation if FORD is found (PATH, conda environment ford or
#    build-check/.ford-venv). Skipped with --lint-only and --doc.
#    With --quick, stop after ctest: no UDTIO build, no FORD and none of the GPU stages below
#    (the everyday check; run the full gate before closing a stage).
# 3. GPU builds (D-019.3, D-019.4b), library and target tests (GPU_TARGET_TESTS): BUILD_DIR-nvptx and
#    BUILD_DIR-amdgcn (gfortran, compiled and linked with offload, both at once) and BUILD_DIR-spir64
#    (ifx, the target tests run on the Intel GPU with memory and time limits). With --gpu-full,
#    NVPTX and AMDGCN compile and link every test. ifx comes from PATH or $ONEAPI_ROOT/setvars.sh.
# 4. Audit of the device images (tools/audit_kernels.sh): the NVPTX, AMDGCN and SPIR64 images above and
#    the NVPTX image of a Release dp library (BUILD_DIR-nvptx-release) call only allowed runtime symbols.
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
doc_only=0
quick=0
gpu_full=0
while [[ "${1:-}" == --* ]]; do
    case "$1" in
        --lint-only) lint_only=1 ;;
        --doc) doc_only=1 ;;
        --quick) quick=1 ;;
        --gpu-full) gpu_full=1 ;;
        *) echo "usage: $0 [--lint-only | --doc | --quick] [--gpu-full] [BUILD_DIR]" >&2; exit 2 ;;
    esac
    shift
done
build_dir="${1:-build-check}"
if (( lint_only + doc_only + quick > 1 )); then
    echo "usage: $0 [--lint-only | --doc | --quick] [--gpu-full] [BUILD_DIR]: one mode at most" >&2
    exit 2
fi

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
# Other sources with a license header: tests, benchmarks (Fortran and Python), the shell tools and the
# scripts of the project skills
mapfile -t header_files < <(listed 'tests/*.f90' 'tests/*.F90' 'tests/*.py' 'benchmarks/*.f90' \
                            'benchmarks/*.F90' 'benchmarks/*.py' 'tools/*.sh' '.claude/skills/*.sh' \
                            '.claude/skills/*.py')
# Library templates subject to the GPU kernel rules
mapfile -t kernel_files < <(printf '%s\n' "${fypp_files[@]}" | grep -E '^(src|include)/' \
                            | xargs -r grep -L '! lint: host-only')
# Every library source (templates and plain Fortran, host-only included): no dynamic dispatch
mapfile -t library_files < <(printf '%s\n' "${fortran_files[@]}" | grep -E '^(src|include)/')

echo "== lint: ${#fypp_files[@]} templates (${#kernel_files[@]} kernel), ${#fortran_files[@]} Fortran sources"

# Formatting (every Fortran source)
report tabs "tab characters" "$(grep -nP '\t' "${fortran_files[@]}" 2>/dev/null)"
# Every source line, comments included, fits the free-form limit (the generated code is checked below)
report line-length "source line longer than 132 columns" "$(perl -CSD -ne '
    chomp; print "$ARGV:$.: ", length($_), " columns\n" if length($_) > 132 && !/!\s*lint:\s*allow\(line-length\)/;
    close ARGV if eof' "${fortran_files[@]}")"

# License headers (sources of the library, tests, benchmarks and tools)
for f in $(printf '%s\n' "${fypp_files[@]}" "${py_files[@]}" "${header_files[@]}" | sort -u); do
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
    } continue { if (eof) { close ARGV; $cont = 0 }' "${library_files[@]}")"
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
# ifx implements exponent() with a runtime routine that does not load with the emulated FP64 of Intel
# GPUs (__isinfd/__isnand): fp_exponent of muscle_math_fp
rule exponent "intrinsic exponent() in kernel code, use fp_exponent" '(?i)(?<![\w%])exponent\s*\(' "${kernel_files[@]}"
# ifx compiles scale() on real32 with double precision operations (for_scale4_v): a sp kernel needs FP64,
# emulated on Intel GPUs without it, where the calls into the emulation are miscompiled on fused EUs
# (PROBLEMAS.md): fp_scale of muscle_math_fp
rule scale "intrinsic scale() in kernel code, use fp_scale" '(?i)(?<![\w%])scale\s*\(' "${kernel_files[@]}"
# ifx 2026.2 keeps only the last element of a leading array item of an array constructor in the host
# code of a declare target procedure ([x, y] with x(3) gives [0, 0, x(3), y]): assign the sections
# (detects whole components and sections; a leading whole array variable needs review)
# ifx implements leadz() with a runtime routine (for_<j|k>ileadz) in device code (tools/audit_kernels.sh)
rule fortran-rt "intrinsic leadz() in kernel code" '(?i)(?<![\w%])leadz\s*\(' "${kernel_files[@]}"
rule arrayctor "array constructor led by an array component or section (ifx host code)" \
    '\[\s*[A-Za-z_]\w*(\([^()]*\))?(%\w+)*%vals\s*[,\]]|\[\s*[A-Za-z_][\w%]*\([^()]*:[^()]*\)\s*,' \
    "${kernel_files[@]}"
rule io "I/O or stop in kernel code" \
    '(?i)^\s*(\d+\s+)?(if\s*\(.*\)\s*)?(print\b|write\s*\(|read\s*\(|stop\b|error\s+stop\b)' \
    "${kernel_files[@]}"

if (( failures > 0 )); then
    echo "== lint: ${failures} rule(s) failed"
    exit 1
fi
echo "== lint: OK"
(( lint_only )) && exit 0

# build_doc <required>: the FORD documentation in ${build_dir}/doc. FORD preprocesses the templates with fypp
# (my_project.md) and must run without errors (warnings are listed in ${build_dir}/ford.log). FORD is looked
# up in the PATH, then in the conda environment `ford` that environment.yml creates, then in the venv
# build-check/.ford-venv (.claude/skills/muscle-revisar-documentacion/scripts/ford_build.sh --install creates
# either); its directory goes first in the PATH, so that FORD runs the fypp of the same environment. Without
# FORD: skipped, or a failure if <required> is 1.
build_doc() {
    local required="$1" ford_bin="" conda_cmd conda_path ford_env seen=" "
    if command -v ford > /dev/null; then
        ford_bin=$(dirname "$(command -v ford)")
    else
        for conda_cmd in conda mamba "$HOME/miniforge3/condabin/conda" "$HOME/miniconda3/condabin/conda" \
                         "$HOME/anaconda3/condabin/conda"; do
            command -v "$conda_cmd" > /dev/null 2>&1 || continue
            # the conda of the PATH and the same one by its absolute path: ask once
            conda_path=$(readlink -f "$(command -v "$conda_cmd")")
            [[ "$seen" == *" $conda_path "* ]] && continue
            seen+="$conda_path "
            ford_env=$("$conda_cmd" env list 2>/dev/null | awk '$1 == "ford" {print $NF}')
            if [[ -n "$ford_env" && -x "$ford_env/bin/ford" ]]; then ford_bin="$ford_env/bin"; break; fi
        done
        [[ -z "$ford_bin" && -x build-check/.ford-venv/bin/ford ]] && ford_bin="$PWD/build-check/.ford-venv/bin"
    fi
    if [[ -z "$ford_bin" ]]; then
        echo "== doc: FORD not found (PATH, conda environment ford, build-check/.ford-venv)$( (( required )) \
            && echo '; install it with .claude/skills/muscle-revisar-documentacion/scripts/ford_build.sh . --install' \
            || echo ', skipped')"
        return "$required"
    fi
    mkdir -p "$build_dir"
    echo "== doc: FORD from $ford_bin (${build_dir}/doc)"
    PATH="$ford_bin:$PATH" PYTHONDONTWRITEBYTECODE=1 ford --quiet --output_dir "${build_dir}/doc" my_project.md \
        > "${build_dir}/ford.log" 2>&1 || {
        echo "== FORD failed (see ${build_dir}/ford.log)"; tail -20 "${build_dir}/ford.log"; return 1; }
}

if (( doc_only )); then
    build_doc 1 || exit 1
    echo "== doc: drift between the documentation and the code (doc_drift.py)"
    drift=$(python3 .claude/skills/muscle-revisar-documentacion/scripts/doc_drift.py .)
    drift_status=$?
    echo "$drift" | grep -v ': 0 finding(s)$'
    (( drift_status == 0 )) || { echo "== check --doc: documentation drift found"; exit 1; }
    echo "== check --doc: OK"
    exit 0
fi

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
for gpu in nvptx nvptx-release amdgcn spir64; do prune_stale_modules "${build_dir}-${gpu}"; done

echo "== build: ${build_dir} (Debug, MUSCLE_REAL_KINDS=sp;dp)"
cmake -S . -B "$build_dir" -G Ninja -DCMAKE_BUILD_TYPE=Debug "-DMUSCLE_REAL_KINDS=sp;dp" \
    > /dev/null || exit 1
# The full check rebuilds from scratch, so that a warning of a file not recompiled since is seen too
clean=()
(( quick )) || clean=(--clean-first)
cmake --build "$build_dir" --target build_all_tests -j "$(nproc)" "${clean[@]}" 2>&1 | tee "$build_dir/build.log"
(( PIPESTATUS[0] == 0 )) || exit 1

# No compiler warnings in the library or the tests (the flags that silence the intentional ones are in
# src/CMakeLists.txt and tests/CMakeLists.txt). gfortran prints the location (the template line, through
# the #line markers) a few lines before each "Warning:"; with --quick only the recompiled files are seen.
warned=$(perl -ne '
    $loc = "$1:$2" if m{^(\S+\.(?:fypp|f90|F90|inc)):(\d+):\d+:\s*$};
    print "$loc: $1\n" if /^Warning: (.*)/' "$build_dir/build.log" | sed "s#^$PWD/##" | sort -u)
if [[ -n "$warned" ]]; then
    echo "LINT [warnings] the build has compiler warnings:"; echo "$warned" | sed 's/^/    /'
    exit 1
fi

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

if (( ! quick )); then
# Formatted output (UDTIO) is disabled for gfortran < 16.3: compile that path anyway
# (library and tests; only the tensor tests run there: the OpenMP tests hit the gfortran bug)
echo "== build: ${build_dir}-udtio (MUSCLE_UDTIO=ON, tensor tests)"
cmake -S . -B "${build_dir}-udtio" -G Ninja -DCMAKE_BUILD_TYPE=Debug -DMUSCLE_UDTIO=ON > /dev/null || exit 1
cmake --build "${build_dir}-udtio" --target build_all_tests -j "$(nproc)" > /dev/null || {
    echo "== build with MUSCLE_UDTIO=ON failed"; exit 1; }
ctest --test-dir "${build_dir}-udtio" -R 'tensors-types-.*_kinds' --output-on-failure -j "$(nproc)" \
    > "${build_dir}-udtio/ctest-udtio.log" || {
    echo "== UDTIO tensor tests failed (see ${build_dir}-udtio/ctest-udtio.log)"; exit 1; }

build_doc 0 || exit 1
fi

ctest --test-dir "$build_dir" --output-on-failure -j "$(nproc)" || failures=$((failures + 1))
if (( quick )); then
    (( failures > 0 )) && { echo "== check --quick: ctest failed"; exit 1; }
    echo "== check --quick: OK (UDTIO, FORD and GPU stages skipped)"
    exit 0
fi

# GPU builds (D-019.3, D-019.4b), each in its own directory, every real kind, of the library and the
# target tests (ctest label `target`, built by the target build_target_tests: the tests that run kernels
# in an OpenMP target region, muscle_target_test_options in tests/CMakeLists.txt):
# - gfortran NVPTX and AMDGCN, at the same time: compiled and linked with offload (no device to run
#   them here); every test with --gpu-full;
# - ifx SPIR64: run on the Intel GPU under memory and time limits (testing.md). Each test binary holds
#   one device image with the kernels of its kinds: with dp the dp emulation of the GPU must be on for
#   the image to load. GPU_SKIP_TESTS: ctest name patterns (regex) not run on the GPU, with the reason.
GPU_TARGETS=(muscle build_target_tests)
# GPU_SKIP_TESTS: tools/gpu_skip_tests.sh
source tools/gpu_skip_tests.sh
gpu_build() { # gpu_build <MUSCLE_GPU> <dir> <jobs> <targets...> [-- cmake args]
    local gpu="$1" dir="$2" jobs="$3"
    shift 3
    local targets=() args=()
    while (( $# )) && [[ "$1" != "--" ]]; do targets+=("$1"); shift; done
    (( $# )) && shift && args=("$@")
    local t0=$SECONDS
    # The offload compilers write large temporary files (device assembly) at link time: keep them
    # in the build directory, not in /tmp (often a small tmpfs)
    mkdir -p "$dir/tmp"
    local -x TMPDIR="$PWD/$dir/tmp"
    [[ "$dir" == /* ]] && TMPDIR="$dir/tmp"
    cmake -S . -B "$dir" -G Ninja -DCMAKE_BUILD_TYPE=Debug "-DMUSCLE_REAL_KINDS=sp;dp" \
        "-DMUSCLE_GPU=$gpu" "${args[@]}" > "$dir/configure.log" 2>&1 || {
        echo "== GPU $gpu: configure failed (see $dir/configure.log)"; tail -5 "$dir/configure.log"; return 1; }
    # GPU_MEM: memory cap of the build (systemd-run scope, no swap)
    local run=()
    [[ -n "${GPU_MEM:-}" ]] && run=(systemd-run --user --scope -p "MemoryMax=$GPU_MEM" -p MemorySwapMax=0 --quiet)
    "${run[@]}" cmake --build "$dir" --target "${targets[@]}" -j "$jobs" > "$dir/build.log" 2>&1 || {
        echo "== GPU $gpu: build failed (see $dir/build.log)"; grep -E '^FAILED' "$dir/build.log" | head; return 1; }
    echo "== GPU $gpu: built ${targets[*]} in $((SECONDS - t0)) s"
}
mkdir -p "${build_dir}-nvptx" "${build_dir}-amdgcn" "${build_dir}-spir64"
gnu_targets=("${GPU_TARGETS[@]}")
(( gpu_full )) && gnu_targets=(build_all_tests)
half=$(( ($(nproc) + 1) / 2 ))
GPU_MEM=6G gpu_build NVPTX "${build_dir}-nvptx" "$half" "${gnu_targets[@]}" > "${build_dir}-nvptx/check.log" 2>&1 &
nvptx_pid=$!
GPU_MEM=6G gpu_build AMDGCN "${build_dir}-amdgcn" "$half" "${gnu_targets[@]}" > "${build_dir}-amdgcn/check.log" 2>&1 &
amdgcn_pid=$!
nvptx_ok=1
amdgcn_ok=1
wait "$nvptx_pid" || { failures=$((failures + 1)); nvptx_ok=0; }
cat "${build_dir}-nvptx/check.log"
wait "$amdgcn_pid" || { failures=$((failures + 1)); amdgcn_ok=0; }
cat "${build_dir}-amdgcn/check.log"

# Audit of the device images (tools/audit_kernels.sh: runtime symbols called from device code, per object),
# in the background while SPIR64 builds and runs: the Debug images of NVPTX and AMDGCN, and the NVPTX
# image of a Release dp library (optimized code: the runtime calls that -O0 leaves are inlined there,
# and the heap temporaries of array expressions only remain there). SPIR64 after its build.
audit_dirs=()
(( nvptx_ok )) && audit_dirs+=("${build_dir}-nvptx")
(( amdgcn_ok )) && audit_dirs+=("${build_dir}-amdgcn")
mkdir -p "${build_dir}-nvptx-release"
(
    rc=0
    GPU_MEM=6G gpu_build NVPTX "${build_dir}-nvptx-release" "$half" muscle \
        -- -DCMAKE_BUILD_TYPE=Release -DMUSCLE_REAL_KINDS=dp || rc=1
    (( rc == 0 )) && audit_dirs+=("${build_dir}-nvptx-release")
    (( ${#audit_dirs[@]} > 0 )) && { systemd-run --user --scope -p MemoryMax=6G -p MemorySwapMax=0 --quiet \
        tools/audit_kernels.sh "${audit_dirs[@]}" || rc=1; }
    exit $rc
) > "${build_dir}-nvptx/audit.log" 2>&1 &
audit_pid=$!

# ifx from the oneAPI environment (ONEAPI_ROOT, default /opt/intel/oneapi), in a subshell
oneapi="${ONEAPI_ROOT:-/opt/intel/oneapi}"
if ! command -v ifx > /dev/null && [[ ! -f "$oneapi/setvars.sh" ]]; then
    echo "== GPU SPIR64: ifx not found (PATH, $oneapi/setvars.sh)"; failures=$((failures + 1))
else
    (
        # setvars.sh reads unset variables
        set +u
        command -v ifx > /dev/null || source "$oneapi/setvars.sh" > /dev/null 2>&1
        set -u
        command -v ifx > /dev/null || { echo "== GPU SPIR64: ifx not found after $oneapi/setvars.sh"; exit 1; }
        # The NVPTX Release build (half of the jobs) may still run: the other half, under the same cap
        GPU_MEM=6G gpu_build SPIR64 "${build_dir}-spir64" "$(( $(nproc) - half > 0 ? $(nproc) - half : 1 ))" \
            "${GPU_TARGETS[@]}" \
            -- -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icx || exit 1
        skip=()
        for t in "${!GPU_SKIP_TESTS[@]}"; do
            echo "   skipped on the GPU: $t (${GPU_SKIP_TESTS[$t]})"
            skip+=("$t")
        done
        exclude="^($(IFS='|'; echo "${skip[*]}"))\$"
        echo "== GPU SPIR64: ctest -L target on the Intel GPU"
        t0=$SECONDS
        systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0 --quiet \
            env OMP_TARGET_OFFLOAD=MANDATORY OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 \
            timeout -k 30 3600 ctest --test-dir "${build_dir}-spir64" -L target -E "$exclude" --timeout 300 --output-on-failure \
            > "${build_dir}-spir64/ctest-gpu.log" 2>&1
        rc=$?
        grep -E 'tests passed|\*\*\*' "${build_dir}-spir64/ctest-gpu.log"
        echo "   ($((SECONDS - t0)) s, log: ${build_dir}-spir64/ctest-gpu.log)"
        systemd-run --user --scope -p MemoryMax=6G -p MemorySwapMax=0 --quiet \
            tools/audit_kernels.sh "${build_dir}-spir64" || rc=1
        exit $rc
    ) || failures=$((failures + 1))
fi
wait "$audit_pid" || failures=$((failures + 1))
cat "${build_dir}-nvptx/audit.log"

if (( failures > 0 )); then
    echo "== check: ${failures} stage(s) failed"
    exit 1
fi
echo "== check: OK"
