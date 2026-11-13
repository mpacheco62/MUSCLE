#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Size and risk report of the device kernels of one executable on the Intel GPU (IGC, SPIR64 build of ifx).
#
#   igc_kernel_report.sh <executable> [kernel-regex] [-- args of the executable]
#
# Runs the executable once on the GPU (memory and time limits, emulated FP64 on, no kernel cache) with the IGC
# shader dumps in a directory next to it, and prints one row per kernel whose name matches kernel-regex
# (default: every kernel):
# - private: per-thread private memory reserved by IGC (.zeinfo private_size). 32768 for a SIMD8 kernel is the
#   4 KB per work item default that IGC uses when it cannot bound the stack and the kernel has no stack calls
#   (D-030): a sign of a run-time-size temporary;
# - spill: register spill size (bytes) of the kernel;
# - stackcalls: has_stack_calls of the .zeinfo (calls with a stack frame);
# - calls: call instructions in the kernel ISA;
# - dyn_alloca: allocas with a run-time size in the optimized LLVM IR of the whole image (temporaries that IGC
#   cannot bound, D-030), with the functions that hold them;
# - fp64emu: calls of the FP64 emulation (__igcbuiltin_dp_*) in the kernel ISA: a sp kernel should have none
#   (scale/exponent intrinsics, a dp literal, ...; D-031).
# Also prints whether the run passed (the PASS/FAIL lines of the executable) and IGC stack overflow messages are
# checked in a second run when the first one loses the device. Environment: ifx (setvars.sh) in the PATH.
set -u
if [[ $# -lt 1 || ! -f "$1" || ! -x "$1" ]]; then
    echo "usage: $0 <executable> [kernel-regex] [-- args of the executable]   (an existing executable)" >&2
    exit 2
fi
exe=$(readlink -f "$1"); shift
regex='.'
if [[ $# -gt 0 && "$1" != "--" ]]; then regex=$1; shift; fi
[[ $# -gt 0 && "$1" == "--" ]] && shift
dump="$(dirname "$exe")/igc_dump_$(basename "$exe")"
rm -rf "$dump"; mkdir -p "$dump"
run() { systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0 --quiet env OMP_TARGET_OFFLOAD=MANDATORY \
    OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 NEO_CACHE_PERSISTENT=0 "$@" timeout -k 30 1200 "$exe" "${ARGS[@]}"; }
ARGS=("$@")
out=$(run IGC_ShaderDumpEnable=1 IGC_DumpToCustomDir="$dump" 2>&1)
echo "run: $(grep -m1 -oE 'PASS[^ ]*|FAIL[^ ]*|DEVICE_LOST|not found in the program|CHECK [^ ]+ .*|failures [0-9]+' <<< "$out" \
    || echo 'no status line')"
if grep -q DEVICE_LOST <<< "$out"; then
    so=$(run IGC_StackOverflowDetection=1 2>&1 | grep -c 'Stack overflow detected')
    echo "device lost; with IGC_StackOverflowDetection=1: $so stack overflow message(s)"
fi
zeinfo=$(ls "$dump"/*.zeinfo 2>/dev/null | head -1)
[[ -n "$zeinfo" ]] || { echo "no IGC dump in $dump (did the image compile?)"; exit 1; }
python3 - "$dump" "$regex" <<'PY'
import glob, re, sys
dump, regex = sys.argv[1], re.compile(sys.argv[2])
kernels = {}
for z in glob.glob(f'{dump}/*.zeinfo'):
    cur = None
    for line in open(z):
        m = re.match(r'\s*- name:\s+(\S+)', line)
        if m and line.startswith('  - name'):
            cur = kernels.setdefault(m.group(1), {})
        elif cur is not None:
            for key in ('private_size', 'spill_size', 'has_stack_calls'):
                m = re.match(rf'\s+{key}:\s+(\S+)', line)
                if m and key not in cur:
                    cur[key] = m.group(1)
asm = {}
for a in glob.glob(f'{dump}/*_entry_*.asm'):
    first = open(a).readline()
    m = re.match(r'//\.kernel\s+(\S+)', first)
    if m:
        text = open(a).read()
        asm[m.group(1)] = (len(re.findall(r'\bcall \(', text)), len(re.findall(r'__igcbuiltin_dp_', text)))
dyn = {}
for ll in glob.glob(f'{dump}/*_optimized.ll'):
    fn = None
    for line in open(ll):
        if line.startswith('define'):
            fn = re.search(r'@([\w$.]+)\(', line).group(1)
        elif re.search(r' alloca [^,]+, i(32|64) %', line):
            dyn[fn] = dyn.get(fn, 0) + 1
print(f"{'kernel':70s} {'private':>8s} {'spill':>6s} {'stackcalls':>10s} {'calls':>5s} {'fp64emu':>7s}")
for k, v in sorted(kernels.items()):
    if k.startswith('Intel_Symbol_Table') or not regex.search(k):
        continue
    c, e = asm.get(k, ('?', '?'))
    flag = '  <- 4 KB/work item default (D-030?)' if v.get('private_size') == '32768' else ''
    print(f"{k[-70:]:70s} {v.get('private_size', '0'):>8s} {v.get('spill_size', '0'):>6s} "
          f"{v.get('has_stack_calls', 'false'):>10s} {c:>5} {e:>7}{flag}")
print('dynamic allocas (function: count):', ', '.join(f'{f}: {n}' for f, n in dyn.items()) or 'none')
PY
echo "dumps: $dump"
