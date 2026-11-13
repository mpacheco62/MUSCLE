#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Registers and spills of the batch kernels of the declared materials (Phase 8, D-019.7c).
#
#   tools/kernel_resources.sh [--amdgcn BUILD_DIR] [--spir64 BUILD_DIR] [--kernel NAME]
#
# Prints one table per back end, one row per material and kind, for the target region of the batch loop
# NAME (default batch_loop_tangent: every point, with the tangent; e.g. batch_loop, batch_loop_dev):
# - AMDGCN (gfortran, BUILD_DIR a MUSCLE_GPU=AMDGCN build, Release recommended: the Debug device code is
#   not optimized): builds muscle_audit_probe (the whole library, -save-temps, as tools/audit_kernels.sh)
#   and reads the GCN assembly of every material: the kernel descriptor (`.vgpr_count`, `.agpr_count`,
#   `.sgpr_count`, `.private_segment_fixed_size`) and, since GCC does not inline the point procedures
#   into the kernel, the call tree of the kernel: the largest VGPR used by any function of the tree and
#   the deepest stack (frame = local variables, spills included, + callee saves + outgoing arguments,
#   summed along the deepest call path), the per work item scratch the device stack must hold.
# - SPIR64 (ifx, BUILD_DIR a MUSCLE_GPU=SPIR64 build with the target tests built): runs each material
#   target test (materials-mat_kernels_target_*, materials-mat_batch_dev_*; not those of
#   tools/gpu_skip_tests.sh) on the Intel GPU under memory and time limits with the IGC shader dumps on,
#   and reads the GRF count and the spill size of each kernel from the IGC assembly (IGC inlines the
#   whole kernel). Each test image holds the kernels of its material in every kind (dp with the emulated
#   FP64). Only the materials with a target test that runs on the GPU have a row.
# NVPTX: needs ptxas of CUDA (not available here), not measured.
# Not part of tools/check.sh. Exit status != 0 on any error.

set -u
cd "$(dirname "$0")/.."

amdgcn=""
spir64=""
kernel="batch_loop_tangent"
while (( $# )); do
    case "$1" in
        --amdgcn) amdgcn="$2"; shift 2 ;;
        --spir64) spir64="$2"; shift 2 ;;
        --kernel) kernel="$2"; shift 2 ;;
        *) echo "usage: $0 [--amdgcn BUILD_DIR] [--spir64 BUILD_DIR] [--kernel NAME]" >&2; exit 2 ;;
    esac
done
[[ -n "$amdgcn$spir64" ]] || { echo "usage: $0 [--amdgcn BUILD_DIR] [--spir64 BUILD_DIR] [--kernel NAME]" >&2; exit 2; }

failures=0
cache_value() { # cache_value <build dir> <variable>
    sed -n "s/^$2:[A-Z]*=//p" "$1/CMakeCache.txt" | head -n 1
}

if [[ -n "$amdgcn" ]]; then
    [[ "$(cache_value "$amdgcn" MUSCLE_GPU | tr '[:lower:]' '[:upper:]')" == AMDGCN ]] || {
        echo "== $amdgcn is not a MUSCLE_GPU=AMDGCN build"; exit 1; }
    work="$amdgcn/audit"
    mkdir -p "$work" "$amdgcn/tmp"
    rm -f "$work"/muscle_audit_probe*
    TMPDIR="$(realpath "$amdgcn")/tmp" systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=0 --quiet \
        cmake --build "$amdgcn" --target muscle_audit_probe > "$work/build.log" 2>&1 || {
        echo "== AMDGCN ($amdgcn): probe build failed (see $work/build.log)"; exit 1; }
    asm=$(ls -v "$work"/muscle_audit_probe.xamdgcn-amdhsa.mkoffload.*.s 2>/dev/null | tail -n 1)
    [[ -f "$asm" ]] || { echo "== AMDGCN ($amdgcn): no device assembly in $work"; exit 1; }
    echo "== AMDGCN ($(cache_value "$amdgcn" CMAKE_BUILD_TYPE), $(cache_value "$amdgcn" MUSCLE_GPU_ARCH)): $kernel"
    KERNEL="$kernel" python3 -I - "$asm" <<'EOF' || failures=$((failures + 1))
import os, re, sys

kernel = os.environ["KERNEL"]
funcs = {}       # name -> dict(frame, vgpr, calls)
meta = {}        # kernel name -> descriptor fields
cur = None
FIELDS = ('vgpr_count', 'agpr_count', 'sgpr_count', 'private_segment_fixed_size')
name_re = re.compile(r'^([A-Za-z_.$][\w.$]*):\s*$')
call_re = re.compile(r'([A-Za-z_.$][\w.$]*)@rel32@lo')
vreg_re = re.compile(r'\bv\[?(\d+)(?::(\d+))?')
frame_re = re.compile(r';\s*(local vars|callee save|outgoing args) size:\s*(\d+)')
md_name = None
types = set()
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('\t.type\t') and '@function' in line:
            types.add(line.split()[1].split(',')[0])
            continue
        m = name_re.match(line)
        if m and m.group(1) in types:
            cur = funcs.setdefault(m.group(1), {'frame': 0, 'vgpr': -1, 'calls': set()})
            continue
        s = line.strip()
        if s.startswith('- .name:'):
            md_name = s.split()[-1]
            meta[md_name] = {}
            continue
        if md_name and s.startswith('.') and ':' in s:
            key, _, val = s[1:].partition(':')
            if key in FIELDS:
                meta[md_name][key] = int(val)
            continue
        if cur is None or s.startswith('.'):
            continue
        m = frame_re.search(s)
        if m:
            cur['frame'] += int(m.group(2))
            continue
        if s.startswith(';'):
            continue
        for c in call_re.findall(s):
            cur['calls'].add(c)
        for a, b in vreg_re.findall(s):
            cur['vgpr'] = max(cur['vgpr'], int(b or a))

def tree(name, seen=()):
    """(largest VGPR index, deepest stack) of the call tree of a function."""
    f = funcs.get(name)
    if f is None or name in seen:
        return -1, 0
    vg, depth = f['vgpr'], 0
    for c in f['calls']:
        cv, cd = tree(c, seen + (name,))
        vg, depth = max(vg, cv), max(depth, cd)
    return vg, f['frame'] + depth

pat = re.compile(r'^__muscle_mat_(\w+)_MOD_' + re.escape(kernel) + r'_(sp|dp)\._omp_fn\.\d+$')
rows = sorted((m.group(1), m.group(2), k) for k in meta for m in [pat.match(k)] if m)
if not rows:
    sys.exit(f'no kernel {kernel} in the assembly')
w = max(len(r[0]) for r in rows)
print(f'{"material":<{w}}  kind  VGPR  AGPR  SGPR  private  tree VGPR  tree stack (B)')
for mat, kind, k in rows:
    d = meta[k]
    vg, depth = tree(k)
    print(f'{mat:<{w}}  {kind:<4}  {d.get("vgpr_count", -1):>4}  {d.get("agpr_count", 0):>4}  '
          f'{d.get("sgpr_count", -1):>4}  {d.get("private_segment_fixed_size", -1):>7}  {vg + 1:>9}  {depth:>14}')
EOF
fi

if [[ -n "$spir64" ]]; then
    [[ "$(cache_value "$spir64" MUSCLE_GPU | tr '[:lower:]' '[:upper:]')" == SPIR64 ]] || {
        echo "== $spir64 is not a MUSCLE_GPU=SPIR64 build"; exit 1; }
    oneapi="${ONEAPI_ROOT:-/opt/intel/oneapi}"
    set +u
    command -v ifx > /dev/null || source "$oneapi/setvars.sh" > /dev/null 2>&1
    set -u
    source tools/gpu_skip_tests.sh
    skip="^($(IFS='|'; echo "${!GPU_SKIP_TESTS[*]}"))\$"
    # name and executable of every material target test that runs on the GPU
    mapfile -t tests < <(ctest --test-dir "$spir64" -L target -R '^materials-mat_(kernels_target|batch_dev)_' \
        -E "$skip" --show-only=json-v1 | python3 -I -c '
import json, sys
for t in json.load(sys.stdin)["tests"]:
    print(t["name"], t["command"][0])')
    (( ${#tests[@]} > 0 )) || { echo "== SPIR64 ($spir64): no material target test"; exit 1; }
    dump="$spir64/igc_dump"
    out="$spir64/kernel_resources_spir64.txt"
    : > "$out"
    echo "== SPIR64 ($(cache_value "$spir64" CMAKE_BUILD_TYPE), Intel GPU): $kernel"
    for t in "${tests[@]}"; do
        name="${t%% *}"
        exe="${t#* }"
        rm -rf "$dump"
        mkdir -p "$dump"
        systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0 --quiet \
            env OMP_TARGET_OFFLOAD=MANDATORY OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 \
            IGC_ShaderDumpEnable=1 IGC_DumpToCustomDir="$(realpath "$dump")" \
            timeout -k 30 300 "$exe" > "$dump/run.log" 2>&1 || {
            echo "   $name: run failed (see $dump/run.log)"; failures=$((failures + 1)); continue; }
        for a in "$dump"/*_entry_*.asm; do
            k=$(grep -m1 -oE '^//\.kernel [^ ]+' "$a" | cut -d' ' -f2)
            grf=$(grep -m1 -oE 'numGRF=[0-9]+' "$a" | cut -d= -f2)
            spill=$(grep -m1 -oE '^//\.spill size [0-9]+' "$a" | awk '{print $3}')
            simd=$(basename "$a" | grep -oE 'simd[0-9]+')
            echo "$k ${grf:-?} ${spill:-0} $simd"
        done >> "$out"
    done
    rm -rf "$dump"
    KERNEL="$kernel" python3 -I - "$out" <<'EOF' || failures=$((failures + 1))
import os, re, sys

kernel = os.environ["KERNEL"]
pat = re.compile(r'_muscle_mat_(\w+)_mp_' + re.escape(kernel) + r'_(sp|dp)__l\d+')
rows = {}
for line in open(sys.argv[1]):
    k, grf, spill, simd = line.split()
    m = pat.search(k)
    if m:
        rows[(m.group(1), m.group(2))] = (grf, spill, simd)
if not rows:
    sys.exit(f'no kernel {kernel} in the IGC dumps')
w = max(len(m) for m, _ in rows)
print(f'{"material":<{w}}  kind  SIMD    GRF  spill (B)')
for (mat, kind), (grf, spill, simd) in sorted(rows.items()):
    print(f'{mat:<{w}}  {kind:<4}  {simd:<6}  {grf:>3}  {spill:>9}')
EOF
fi

exit $(( failures > 0 ))
