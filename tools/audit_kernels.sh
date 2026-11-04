#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Audit of the kernel path in the device images of the library (Phase 8, D-019).
#
#   tools/audit_kernels.sh [--lint] BUILD_DIR...
#
# Each BUILD_DIR is a GPU build (MUSCLE_GPU = NVPTX, AMDGCN or SPIR64) whose library is already built.
# The audit lists, for each object of the library, the runtime symbols its device code calls that are
# not on the allowed list below, and fails if there is any:
# - NVPTX, AMDGCN (gfortran): builds the target muscle_audit_probe (src/CMakeLists.txt: the whole
#   library linked into a program with a target region, -save-temps) and reads the device assembly
#   (PTX, GCN) that the offload compiler wrote: calls to functions the image does not define. Each
#   device function is attributed to the library object that defines the same host symbol.
# - SPIR64 (ifx): extracts the device bitcode of each library object (clang-offload-bundler) and
#   lists its undefined symbols (llvm-nm) that no library object defines.
# Allowed: device math library (libm names), OpenMP device runtime, memcpy/memset, compiler support
# routines (__powidf2, __divsi3, ...), SPIR-V/OpenCL builtins; in a Debug build also the -fcheck
# diagnostics (_gfortran_runtime_error_at/_warning_at) and the matmul routines (not inlined at -O0);
# for SPIR64 the routine of scale() (for_scale<4|8>_v, see ALLOWED_SPIR64).
# Everything else (Fortran runtime: _gfortran_*, for_*; malloc/free; I/O; stop; ieee_* entry/exit)
# is reported. --lint also runs the source-level rules of tools/check.sh --lint-only (device,
# alloc, io, exponent, ...), which tools/check.sh already runs itself.
# Exit status != 0 on any finding or error.

set -u
cd "$(dirname "$0")/.."

lint=0
if [[ "${1:-}" == "--lint" ]]; then lint=1; shift; fi
(( $# > 0 )) || { echo "usage: $0 [--lint] BUILD_DIR..." >&2; exit 2; }

failures=0
if (( lint )); then
    tools/check.sh --lint-only || failures=$((failures + 1))
fi

cache_value() { # cache_value <build dir> <variable>
    sed -n "s/^$2:[A-Z]*=//p" "$1/CMakeCache.txt" | head -n 1
}

for dir in "$@"; do
    if [[ ! -f "$dir/CMakeCache.txt" ]]; then
        echo "== audit: $dir is not a build directory"; failures=$((failures + 1)); continue
    fi
    gpu=$(cache_value "$dir" MUSCLE_GPU | tr '[:lower:]' '[:upper:]')
    config=$(cache_value "$dir" CMAKE_BUILD_TYPE)
    lib="$dir/src/libmuscle.a"
    work="$dir/audit"
    mkdir -p "$work"
    case "$gpu" in
        NVPTX|AMDGCN)
            # Stale assembly of an earlier probe must not be read
            rm -f "$work"/muscle_audit_probe*
            cmake --build "$dir" --target muscle_audit_probe > "$work/build.log" 2>&1 || {
                echo "== audit $gpu ($dir): probe link failed (see $work/build.log)"; tail -5 "$work/build.log"
                failures=$((failures + 1)); continue; }
            if [[ "$gpu" == NVPTX ]]; then
                mapfile -t asm < <(ls "$work"/muscle_audit_probe.xnvptx-none.mkoffload*.s 2>/dev/null | head -n 1)
            else
                mapfile -t asm < <(ls -v "$work"/muscle_audit_probe.xamdgcn-amdhsa.mkoffload.*.s 2>/dev/null | tail -n 1)
            fi
            [[ -f "${asm[0]:-}" ]] || {
                echo "== audit $gpu ($dir): no device assembly in $work"; failures=$((failures + 1)); continue; }
            nm --defined-only -A "$lib" 2>/dev/null | awk '$(NF-1) ~ /^[TtWw]$/ {split($1, a, ":"); print a[2], $NF}' \
                > "$work/host_symbols.txt"
            mode=gnu
            inputs=("${asm[0]}" "$work/host_symbols.txt")
            ;;
        SPIR64)
            tools_dir=""
            for c in "$(dirname "$(command -v ifx 2>/dev/null || echo /nonexistent)")/compiler" \
                     "${ONEAPI_ROOT:-/opt/intel/oneapi}/compiler/latest/bin/compiler"; do
                [[ -x "$c/clang-offload-bundler" && -x "$c/llvm-nm" ]] && { tools_dir="$c"; break; }
            done
            [[ -n "$tools_dir" ]] || {
                echo "== audit SPIR64 ($dir): clang-offload-bundler/llvm-nm of oneAPI not found"
                failures=$((failures + 1)); continue; }
            rm -rf "$work/objects"
            mkdir -p "$work/objects"
            (cd "$work/objects" && ar x "$(realpath "$OLDPWD/$lib" 2> /dev/null || echo "$lib")") || {
                echo "== audit SPIR64 ($dir): cannot extract $lib"; failures=$((failures + 1)); continue; }
            : > "$work/device_symbols.txt"
            for o in "$work"/objects/*.o; do
                b="${o%.o}.bc"
                "$tools_dir/clang-offload-bundler" -type=o -targets=openmp-spir64 -input="$o" -output="$b" \
                    -unbundle 2> /dev/null || continue
                "$tools_dir/llvm-nm" --defined-only "$b" 2> /dev/null | awk -v o="$(basename "$o")" '{print o, "D", $NF}'
                "$tools_dir/llvm-nm" --undefined-only "$b" 2> /dev/null | awk -v o="$(basename "$o")" '{print o, "U", $NF}'
            done >> "$work/device_symbols.txt"
            mode=spir64
            inputs=("$work/device_symbols.txt")
            ;;
        *)
            echo "== audit: $dir has MUSCLE_GPU='$gpu' (not a GPU build)"; failures=$((failures + 1)); continue ;;
    esac

    MODE="$mode" CONFIG="$config" GPU="$gpu" python3 -I - "${inputs[@]}" <<'EOF' || failures=$((failures + 1))
import collections, os, re, sys

mode, config, gpu = os.environ["MODE"], os.environ["CONFIG"], os.environ["GPU"]
# Runtime symbols allowed in device code
ALLOWED = re.compile(r"""^(?:
    # device math library
    (?:a?(?:sin|cos|tan)h?|atan2|sincos|exp|exp2|expm1|log|log10|log2|log1p|pow|sqrt|cbrt|hypot|fabs|fmod
       |fmin|fmax|floor|ceil|trunc|round|lround|llround|rint|lrint|nearbyint|copysign|scalbn|ldexp|frexp
       |erf|erfc|tgamma|lgamma|fma)f?
    # OpenMP device runtime
    |GOMP_\w+|gomp_\w+|omp_\w+|__kmpc_\w+|ompx_\w+|__omp_\w+
    |memcpy|memset|memmove
    # compiler support routines (integer division, powers, wide multiplication, ...)
    |__(?:pow|mul|div|mod|udiv|umod|neg|ashl|ashr|lshr|cmp|ucmp|clz|ctz|popcount|ffs|parity|bswap)\w*\d
    # SPIR-V / OpenCL builtins and ITT instrumentation (ifx)
    |_Z\d+__spirv_\w+|__spirv_\w+|_Z\d+get_\w+|__itt_offload_\w+
)$""", re.X)
# Debug builds: -fcheck diagnostics and the matmul routine (gfortran inlines matmul only when optimizing)
ALLOWED_DEBUG = re.compile(r"^_gfortran_(?:runtime_error_at|runtime_warning_at|matmul_r[48])$")

# SPIR64: the device routine of scale() (yield criteria, spectral_decomposition). Replacing scale() by
# products with powers of 2 (same bits on the CPU) makes IGC lose the device in three sp material kernels
# that pass with it (2026-11-03, PROBLEMAS.md): kept until a driver/compiler without that bug.
ALLOWED_SPIR64 = re.compile(r"^for_scale[48]_v$")

def allowed(sym):
    return (bool(ALLOWED.match(sym)) or (config == "Debug" and bool(ALLOWED_DEBUG.match(sym)))
            or (gpu == "SPIR64" and bool(ALLOWED_SPIR64.match(sym))))

def base(name):
    """Host name of a device function: without the suffixes of outlined regions and clones."""
    return re.split(r"[.$]_omp_fn|\.(?:constprop|isra|part|cold|localalias)\b|\$", name)[0]

findings = collections.defaultdict(lambda: collections.defaultdict(set))  # object -> symbol -> callers
if mode == "gnu":
    asm, host = sys.argv[1], sys.argv[2]
    owner = {}
    for line in open(host):
        obj, sym = line.split()
        owner.setdefault(sym, obj)
    defined, calls, cur = set(), collections.defaultdict(set), None
    if gpu == "NVPTX":
        fdef = re.compile(r"^// BEGIN (?:GLOBAL )?FUNCTION DEF: (\S+)")
        fcall = re.compile(r"\bcall(?:\.uni)?\s+(?:\([^)]*\)\s*,\s*)?([A-Za-z_$][\w$.]*)")
    else:
        fdef = re.compile(r"^\s*\.type\s+([^,\s]+),\s*@function")
        ddef = re.compile(r"^\s*\.type\s+([^,\s]+),\s*@object")
        fcall = re.compile(r"(?<![\w$.])([A-Za-z_$][\w$.]*)@(?:rel32|gotpcrel32|rel64|abs32)@lo")
    ddef = locals().get("ddef")
    for line in open(asm, errors="replace"):
        if ddef:
            m = ddef.match(line)
            if m:
                defined.add(m.group(1))
                continue
        m = fdef.match(line)
        if m:
            cur = m.group(1)
            defined.add(cur)
            continue
        if cur:
            for m in fcall.finditer(line):
                calls[cur].add(m.group(1))
    for caller, callees in calls.items():
        for sym in callees:
            if sym in defined or allowed(sym) or sym.startswith("."):
                continue
            obj = owner.get(base(caller), "(no library object: %s)" % base(caller))
            findings[obj][sym].add(base(caller))
else:
    defined, undefined = set(), collections.defaultdict(set)
    for line in open(sys.argv[1]):
        obj, kind, sym = line.split(None, 2)
        sym = sym.strip()
        (defined.add(sym) if kind == "D" else undefined[obj].add(sym))
    for obj, syms in undefined.items():
        for sym in syms:
            if sym not in defined and not allowed(sym) and not sym.startswith("llvm."):
                findings[obj][sym]

if findings:
    n = sum(len(s) for s in findings.values())
    print("== audit %s (%s): %d runtime symbol(s) in device code not allowed:" % (gpu, config, n))
    for obj in sorted(findings):
        for sym in sorted(findings[obj]):
            callers = sorted(findings[obj][sym])
            more = " (+%d)" % (len(callers) - 4) if len(callers) > 4 else ""
            print("    %s: %s%s" % (obj, sym, (" <- " + ", ".join(callers[:4]) + more) if callers else ""))
    sys.exit(1)
print("== audit %s (%s): device code calls only allowed runtime symbols" % (gpu, config))
EOF
done

exit $(( failures > 0 ))
