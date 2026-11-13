# Wrong results in an OpenMP target region on Intel GPU with FP64 emulation (ifx + IGC)

## Summary

A Fortran OpenMP offload kernel built with `ifx -fiopenmp -fopenmp-targets=spir64` returns wrong values for some
work-items on an Intel UHD Graphics 770 when FP64 emulation is enabled, at `-O1` and above. The same binary with
`OMP_TARGET_OFFLOAD=DISABLED`, the `-O0` build, and gfortran on the host all give the correct result.

Two issues are involved:

1. **ifx:** the intrinsic `SCALE` on `real(real32)` arguments is compiled with double-precision operations in
   device code. A kernel with nothing but `y = scale(x, -3)` (single precision) fails to load on a GPU without
   native FP64 ("Double type is not supported on this platform") unless FP64 emulation is enabled.
2. **IGC (FP64 emulation):** with `OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1`, the attached kernel,
   which uses `SCALE` and also indexes small local arrays with indices read from `PARAMETER` arrays, produces
   wrong results in about 19 % of the work-items.

## Reproducer

`mre_igc_scale_fp64emu.F90` (and the same with designed inputs, `mre_igc_scale_fp64emu_patterns.F90`): single file, single precision, no external library (459 lines, about 430 of code).
Each work-item runs an elastic-plastic return mapping on one point. With the default input 233 of the 256 points are
elastic (the loop exits at its first iteration, `iters` = 0) and 23 are plastic (Newton iterations through the
branch with the eigen-decompositions).

```
ifx -O2 -fiopenmp -fopenmp-targets=spir64 -fpp mre_igc_scale_fp64emu.F90 -o mre
OMP_TARGET_OFFLOAD=MANDATORY OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 ./mre
OMP_TARGET_OFFLOAD=DISABLED ./mre
```

| Build / run | Result |
|---|---|
| `-O0`, GPU | `bad points: 0 of 256` (correct) |
| `-O1`, `-O2`, `-O3`, GPU | `bad points: 48 of 256`, the same every run |
| any level, `OMP_TARGET_OFFLOAD=DISABLED` | `bad points: 0 of 256` |
| gfortran `-fopenmp` (host) | `bad points: 0 of 256` |

In the wrong work-items the stress is bit-identical to the host result, but `iters` is 50 instead of 0, as if the
loop counter were corrupted.

### Reduced version

`mre_igc_reduced.F90` (280 lines, about 245 of code) keeps the failure with the designed input (one plastic lane per
block of 16): 128 of 256 wrong at `-O2`/`-O3` with both driver stacks, 0 at `-O0`, 0 with SCALE replaced by a
multiplication, 0 with literal rotation indices, 0 with `IGC_EmulationFunctionControl=1`. It drops the closed-form
eigenvalues, the yield function and the 6x6 tangent of the original; SCALE remains only in the plastic branch, so
the work-items that come back wrong never call an emulation function. Further reductions are blocked by the
sensitivity of the failure to register allocation: changes that only replace a value loaded at run time by a
constant (the elastic constants, the initial plastic strain) or move it to another derived type make it disappear,
and so do small rewrites of the eigenvector update. Synthetic kernels written from scratch with the same
ingredients (a divergent branch per block of 16, Jacobi rotations indexed through `PARAMETER` tables, SCALE, up to
3.7 KB of spill) did not fail.

Side observation while building those kernels: with enough unrolled code (a synthetic kernel, not included here, with 7 or more
local 3x3 matrices updated through SCALE, without `IGC_EmulationFunctionControl=1`) IGC produces a module without the kernel and without a build error; the
runtime then reports `kernel '__omp_offloading_..._batch__l89' not found in the program` and the program crashes.

## Controls (each removes the failure; same arithmetic)

| Change | GPU result |
|---|---|
| `scale(x, -e)` → `x*2.0**(-e)`, three places | 0 bad points; the image no longer needs FP64 |
| Jacobi rotations indexed through `P_OF(i)`, `Q_OF(i)`, `R_OF(i)` → the same three rotations with literal indices | 0 bad points |
| removing the plastic branch (never executed for this input) or the 6×6 tangent arguments of the kernel | 0 bad points |

The wrong work-items are elastic: they never execute the plastic branch, but other work-items do. Results that
depend only on code a work-item does not execute point to code generation, not to the arithmetic.

## Exact trigger (designed inputs)

The build `-DPATTERN=<n>` of `mre_igc_scale_fp64emu_patterns.F90` chooses which work-items take the "plastic" branch
(the branch with the eigen-decompositions) inside every aligned block of 16 work-items. The kernel runs in SIMD8
(`IGC_ForceOCLSIMDWidth=8` gives the same result), so a block of 16 is two adjacent SIMD8 hardware threads.

| Plastic work-items in each block of 16 | Wrong work-items on the GPU |
|---|---|
| one in the first half (lane 3) | all 8 lanes of the second half, in every block |
| one in the second half (lane 11) | all 8 lanes of the first half |
| one in each half | none |
| none / all | none |
| the whole first half | all 8 lanes of the second half |
| one in the first half, odd blocks only | the second half of the odd blocks only |

Rule: **when exactly one of two adjacent SIMD8 threads executes the divergent branch, all 8 lanes of the other
thread return a wrong value** (here the iteration counter: 50 instead of 0; the stress is correct). The work-group
size does not matter (`thread_limit(8)` to `thread_limit(128)` give the same blocks), which points to the fused EU
pairs of Gen12LP (the dumps report platform TGLLP).

What the dumps (`IGC_ShaderDumpEnable=1`) show:

- `SCALE` on `real(real32)` is lowered to `__igcbuiltin_sp_to_dp`, `__igcbuiltin_dp_mul` and `__igcbuiltin_dp_to_sp`:
  emulated FP64 builtins, `noinline`, called from the kernel (36 `call` instructions; with `SCALE` replaced by a
  real32 multiplication there are no calls and no failure).
- The run-time indexing leaves five private arrays (`alloca`) in the optimized IR (the two 3x3 matrices of the
  Jacobi rotations, the two eigenvector matrices and the `[3 x i16]` index table), against one with literal indices;
  the register spill goes from 2112 to 3392 bytes.

| Condition | Effect on the designed input that fails (one plastic lane per block) |
|---|---|
| `SCALE` → `x*2.0**(-e)` (no emulated FP64) | 0 wrong |
| literal indices in the rotations | 0 wrong |
| `-O0` | 0 wrong |
| `IGC_EmulationFunctionControl=1` (emulation functions inlined) | 0 wrong |
| `IGC_EmulationFunctionControl=0, 2, 3` | 128 wrong (unchanged) |
| `IGC_EnableCallWA=0/1`, `IGC_SeparateSpillPvtScratchSpace=0/1`, `IGC_ForceInlineStackCallWithImplArg=1` | unchanged |

## Machine-code analysis: a fused-EU hazard

`IGC_DisableEuFusion=1` makes every reproducer correct (`mre_igc_reduced.F90` 128 → 0, the original 48 → 0 and
128 → 0). The dumps with and without the flag contain the same instructions (the `.asm` files are identical apart
from comments); the only difference is the attribute `require_disable_eufusion: true` in the kernel's `.zeinfo`, so
the runtime dispatches the same code with the EUs not fused. The miscompile is therefore not in the arithmetic: the
same code is correct on unfused EUs and wrong on the fused EU pairs of Gen12LP.

IGC knows the hardware issue (visa/FlowGraph.cpp, `findNestedDivergentBBs`): on fused EUs the "fused mask" can stay
wrong after nested divergent control flow (a `goto` out of a divergent loop, a `join`), and the EU whose channels are
all off keeps executing NoMask instructions. The workaround (`-noMaskWA`, `-fusedCallWA 1` in the vISA options) marks
the affected blocks and guards NoMask instructions with the per-EU "any channel on" flag (`(W&f0.0.any8h)`). In the
failing kernel all blocks of the return-mapping loop and of the called emulation functions are marked `[NoMaskWA]`
(307 blocks, 630 guarded instructions), but several kinds of NoMask instructions stay unguarded in those blocks and
are executed by the "off" EU: flag saves and restores (`(W) mov r60.0 f0.0`, `(W) mov f0.0 r60.0`), spill and fill
messages to scratch (`(W) send.dc0 ... r126`), the stack-call prologue and epilogue of the emulation functions
(`(W) ... (1|M16)` SP/FP updates, the save and restore of r127 with the return address), and the SP/FP copies of the
leaf builtins (`(W) mov (4|M0) r59.4 r127.0`).

The wrong value matches a wrong path through the loop exits of the elastic thread: the loop counter of the 8 wrong
lanes is set by `mov (8|M0) r78 = 0x33` (i = 51, `iters` = 50), which is only in block `_0_100`, the join reached
when the counter is exhausted; the elastic lanes leave the loop at block `_0_057` with r78 = i = 1. The stress of
those lanes is right, so they did not execute the plastic update. Why the three conditions matter, from the dumps:
SCALE adds the stack calls into the emulation functions (11 `call` in the kernel, each with NoMask prologue,
epilogue and caller-save spills) inside the divergent loop; the run-time indices and -O2 raise the register
pressure (2176 bytes of spill in this file), which puts the loop counter and other loop-carried values in scratch
behind NoMask fills and spills; `IGC_EmulationFunctionControl=1` inlines the emulation and removes the calls. Without
the hardware documentation of the fused mask it was not possible to name the single instruction that goes wrong
(`IGC_ShaderOverride` is not honoured by the release build of IGC, so the ISA could not be edited and re-run).

Workarounds that give correct results: `IGC_DisableEuFusion=1`, `IGC_EmulationFunctionControl=1`.

## Related reports

- intel/intel-graphics-compiler#397: FP64 emulation left in call form at a required SIMD32 returns an unrelated
  live value of the caller; `IGC_EmulationFunctionControl=1` restores correct results (Arc A770, IGC 2.30.1). Our
  case has no required SIMD size (SIMD8) and comes from ifx/OpenMP, with the same workaround.
- intel/intel-graphics-compiler#413: silent miscompile of software-emulated FP64 at -O1/-O2, correct at -O0
  (Alder Lake-P, DPC++); closed in October 2026 by a commit whose relation to the bug is not stated.

## Environment

- ifx (IFX) 2026.2.0 20261001, `-fiopenmp -fopenmp-targets=spir64`
- Intel UHD Graphics 770 (integrated, Raptor Lake, no native FP64; IGC platform TGLLP)
- Reproduced with two driver stacks, same results (48 wrong with the default input, 128 with `-DPATTERN=1`, 0 with
  `IGC_EmulationFunctionControl=1`):
  - intel-opencl-icd / libze-intel-gpu1 26.05.37020, libigc2 2.28.4 (distribution packages);
  - compute-runtime 26.40.40196.8 with IGC 2.42.6 and gmmlib 22.10.2, the latest releases (2026-11-03 and
    2026-10-26), loaded from the extracted packages (`LD_LIBRARY_PATH`, `ZE_ENABLE_ALT_DRIVERS`,
    `NEO_CACHE_PERSISTENT=0`); the IGC dump confirms a recompilation (spill size 3328 bytes against 3392).
- Linux 6.18 (WSL2), Intel Core i9-13900K

## How it was found

The kernel is a reduction of a material model of a constitutive library: its batch kernel gave wrong results on
the GPU for an eigen-decomposition written with indices read from `PARAMETER` arrays, while the same routine was
correct alone and in a point kernel. The reduction kept, at every step, a correct result on the host
(`OMP_TARGET_OFFLOAD=DISABLED`), builds without warnings about undefined or unset variables (`gfortran -Wall`:
uninitialized, INTENT(OUT) not set; `ifx -warn all`: #6178, #6843; the remaining warnings are unused arguments and
variables), and a host run with
`-finit-real=snan -ffpe-trap=invalid,zero,overflow`, so no undefined behaviour was introduced.
