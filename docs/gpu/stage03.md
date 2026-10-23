# Validated CPU/NVIDIA baseline before Fypp integration

The stage 3 feature probes passed on 2026-10-23 using GNU Fortran 12.3,
Fypp 3.3 and OpenMP Target NVPTX (`sm_80`) on an NVIDIA GeForce RTX 5060.
The original MUSCLE sources and CMake build are unchanged in this checkpoint.

The [terminal run summary](../../validation/evidence/stage03/summary.json)
records 2/2 NVIDIA tests, 2/2 tests in each of three CPU profiles, 51/51 MUSCLE
Release and 51/51 Debug tests, and 462 Fypp tests. These are historical
measurements from the original workspace, whose absolute paths remain in the
raw logs; use the commands below to reproduce feature tests from this checkout.

The [CPU/GPU comparison](../../validation/evidence/stage03/nvidia-gcc12-comparison.json)
contains 257 points, exact agreement of identifiers/statuses, and no differences
outside the preset real64 tolerances. It records one device, initial-device
false, zero calculation points on host and zero reduction iterations on host.
The maximum stress difference is 2.84e-14 and checksum difference 5.82e-11.
No tolerance was relaxed. A successful build or host fallback cannot pass.

## Reproduce

Provide Fypp 3.3 on PATH or in a sibling `../fypp` checkout (release commit
`0db20548ca930cd4f9a201bed5f5acd4fc81bc73`). No network downloads occur at build time.
From the MUSCLE repository root:

```sh
cmake -S validation/openmp-fypp --preset nvidia-gcc12
cmake --build build/profiles/nvidia-gcc12 --parallel 4
ctest --test-dir build/profiles/nvidia-gcc12 --output-on-failure -V
```

Use `cpu-serial`, `cpu-openmp` or `cpu-target-host` for CPU feature checks.
The NVIDIA compiler wrapper limits PATH inside the compiler process to avoid
CUDA ptxas incompatibilities; it does not change the user's shell.
A terminal with GPU access is required for the NVIDIA tests.

## Scope and next stage

Stage 4 will integrate real Fypp generation into the library build while
preserving the public CPU API and sharing numerical implementations with
concrete specializations. The complete constitutive GPU solve belongs to
stage 5. AMD hardware/toolchain validation remains pending. CPU-001
(Drucker–Prager tension nonconvergence) remains a separate numerical issue.

The [CPU reference program](../../validation/reference/cpu_reference.f90)
exports 768 Von Mises states with history, stresses, tangent and residuals.
Its stored Release CSV and exact Debug/Release comparison are the numerical
baseline for stage 4. The one-point Drucker–Prager reproducer is also retained.
