# MUSCLE: Mechanics User-defined Solids Constitutive Library & Extensions

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Fortran](https://img.shields.io/badge/Language-Fortran_2008-orange.svg)](https://gcc.gnu.org/fortran/)

**MUSCLE** is a high-performance Fortran 2008 library for the development, testing and implementation of
constitutive material models in Computational Solid Mechanics (elasticity, plasticity, viscoplasticity,
viscoelasticity, finite strain kinematics).

The sources are [FYPP](https://fypp.readthedocs.io) templates: every tensor type, model and material is
generated in single (`sp`) and double (`dp`) precision from one source, with static dispatch only (no
polymorphism, `pure` procedures), so that the same code runs on CPUs (OpenMP) and can be compiled for GPUs
(OpenMP offload). Constitutive models are combined into **declared materials**: you list the combinations
you need in a small Python registry and the build generates one Fortran module per material, with a point
API for the user subroutines of FE codes (Abaqus UMAT, ANSYS USERMAT, LS-DYNA umat41) and a batch API over
arrays of points.

---

## Key Features

### 1. Tensor Engine (`muscle_tensors`)
Voigt-notation tensor algebra with manually unrolled loops, in 3D and 2D, `sp` and `dp`:
*   **2nd order tensors**, symmetric and general (`ten_3D2Osym_dp`, `ten_2D2O_sp`, ...).
*   **4th order tensors** with minor, major and full symmetries (`ten_3D4O3sym_dp`, ...).
*   **Identities**, plain and scaled (`iden_2O`, `iden_4O4T`, `iden_2OS_dp`, ...).
*   **Overloaded operators:** `+ - * / =`, double contraction (`.ddot.`), dyadic products (`.tdot.`,
    `.tdotsym.`), `.transform.`, deviatoric part (`.dev.`) and inverse (`.inv.`).
*   Names without kind suffix (`ten_3D2Osym`, ...) are the `dp` types.

### 2. Model families (free `pure` functions over concrete parameter types)
*   **Elasticity:** linear isotropic; frozen viscoelastic elasticity driven by a viscoelastic model.
*   **Viscoelasticity:** generalized Maxwell and Kelvin-Voigt.
*   **Yield criteria:** von Mises, Hill48, Drucker-Prager, CPB06, CPB06ex2 (Cazacu-Plunkett-Barlat) and
    Yld2004-18p, with analytical first and second derivatives.
*   **Isotropic hardening:** bilinear, Ludwik, Swift, modified Voce.
*   **Overstress viscoplasticity:** Perzyna, Chaboche.
*   **Viscoplastic flow stress laws:** Johnson-Cook (and full Johnson-Cook), Rusinek-Klepaczko (standard and
    modified), Nemat-Nasser-Li, Voyiadjis-Almasri.
*   **Kinematics:** small strain, total and updated Lagrangian, material and spatial logarithmic,
    corotational (strain rotated by the FE code).

### 3. Return-mapping solvers (bodies included by every material)
*   **Closest point projection** (Newton-Raphson, consistent tangent), with overstress (`closest_point_perzyna`)
    and consistency (`closest_point_time`) viscoplasticity.
*   **Monolithic closest point** (`closest_point_2`, 7x7 Newton with an own LU solver).
*   **Cutting plane** (explicit, continuum elastoplastic tangent).

### 4. Declared materials (`muscle_mat_<name>`)
For each material of the registry, in `sp` and `dp`:
*   `mat_<name>_<k>`: the material type (one parameter component per model, `tol`, `iter_nw`).
*   `mat_init_<name>(mat, props, status)` and `MAT_<NAME>_NPROPS`: the material from a flat parameter array.
*   `muscle_point_<name>(...)`: point API over flat arrays (`strain(6)` or `F_n(9)`, `F_np1(9)`), for FE
    user subroutines.
*   `muscle_update_<name>(...)`: batch API over structure-of-arrays data (`eps(npts, 6)`), OpenMP on the
    CPU, optional index list for multiphase models.
*   `MAT_<NAME>_NSTATE`, `mat_pack_<name>`/`mat_unpack_<name>`, `mat_rotate_hsv_<name>`: the history layout.

---

## Requirements

*   **Compiler:** a Fortran 2008 compiler with OpenMP (the GPU build needs OpenMP 5.1 offload). Older
    versions than the tested ones below may work but are not verified.
*   **Build:** CMake 3.20+ (Ninja recommended) and **FYPP 3.1+** (`pip install fypp`).
*   **Python 3** (FYPP and the material registry).
*   **BLAS & LAPACK** (or Intel MKL), required by CMake.
*   Optional: [FORD](https://forddocs.readthedocs.io) for the documentation. `environment.yml` installs
    FYPP and FORD.

Tested toolchain (Linux):

| Tool | Version | Scope |
|------|---------|-------|
| GNU Fortran | 15.2 | CPU build and every test (`sp;dp`, Debug); NVPTX and AMDGCN offload builds compile and link (not run on NVIDIA/AMD hardware) |
| Intel `ifx` (oneAPI) | 2026.2 | SPIR64 offload build; target tests run on an Intel UHD 770 GPU (some cases excluded by IGC bugs, `tools/gpu_skip_tests.sh`) |
| CMake | 4.2 (3.20+ required) | |
| Ninja | 1.13 | |
| FYPP | 3.1 and 3.3 | |
| Python | 3.10 and 3.14 | |

Formatted output of the tensors (UDTIO) is disabled by default with GNU Fortran < 16.3, whose UDTIO is
not thread-safe under OpenMP (`MUSCLE_UDTIO=AUTO`).

---

## Building

```bash
git clone https://github.com/mpacheco62/muscle.git
cd muscle

# Configure and build the library
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Tests (not built by the default target)
cmake --build build --target build_all_tests
ctest --test-dir build --output-on-failure
```

Options:

| Option | Effect |
|--------|--------|
| `-DCMAKE_BUILD_TYPE=Debug\|Release\|Profile\|Profile_unittests` | Build type (`Profile_unittests`: coverage, target `coverage_report`). |
| `-DMUSCLE_MATERIALS=<path>` | Your material registry (default: the test registry `tools/fypp/muscle_materials.py`). |
| `-DMUSCLE_REAL_KINDS=dp` | Generate only double precision (default `sp;dp`). |
| `-DMUSCLE_BUILD_BENCHMARKS=ON` | Build the benchmarks (`benchmarks/`, target `build_all_benchmarks`). |
| `-DMUSCLE_GPU=OFF\|NVPTX\|AMDGCN\|SPIR64` | OpenMP offload target of the kernels (default `OFF`, see below). |
| `-DMUSCLE_GPU_ARCH=<arch>` | Device architecture (NVPTX: e.g. `sm_80`; AMDGCN: default `gfx90a`; not used with SPIR64). |
| `-DMUSCLE_UDTIO=AUTO\|ON\|OFF` | Formatted output of the tensors (default `AUTO`, see Requirements). |

### GPU build (OpenMP offload)

| `MUSCLE_GPU` | Compiler | What you need installed |
|--------------|----------|-------------------------|
| `NVPTX` | GNU Fortran | GCC built with the nvptx offload backend (e.g. the distribution package `gcc-<version>-offload-nvptx`) |
| `AMDGCN` | GNU Fortran | GCC built with the amdgcn offload backend (`gcc-<version>-offload-amdgcn`) |
| `SPIR64` | Intel `ifx` | oneAPI with its environment loaded (`source /opt/intel/oneapi/setvars.sh`) and the Intel GPU driver (Level Zero) |

```bash
# NVIDIA (gfortran with nvptx offload)
cmake -S . -B build-nvptx -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_GPU=NVPTX
# Intel GPU (ifx)
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-spir64 -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_GPU=SPIR64 \
      -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icx
```

CMake checks the combination at configure time by linking a small program with a target region. Run
GPU programs with `OMP_TARGET_OFFLOAD=MANDATORY` to fail instead of falling back to the host; on Intel GPUs
without native double precision, `OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1` enables the
emulated FP64. The batch API on the device is described in [FEA integration](pages/integration.md).

Documentation: `cmake --build build --target doc` (output in `build/doc`, when CMake finds `ford`), or
`ford my_project.md` from the repository root.

---

## Declaring materials

A registry is a Python file that defines `MATERIALS` (see `tools/fypp/muscle_materials.py`, which also
documents every model and solver it covers):

```python
# my_materials.py
from muscle_meta import Material

MATERIALS = (
    Material('steel', 'closest_point', elas='linear', yield_='vonmises', hard='swift',
             doc='von Mises with Swift hardening'),
    Material('alu_rate', 'closest_point_time', elas='linear', yield_='hill48', vp='jc', hard='swift',
             doc='Hill48, Johnson-Cook, Swift hardening'),
    Material('steel_ul', 'closest_point', elas='linear', yield_='vonmises', hard='bilinear',
             kin='updated_lagrangian', tangent_rate='jaumann', doc='Updated Lagrangian'),
)
```

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_MATERIALS=$PWD/my_materials.py
python3 tools/fypp/muscle_meta.py props steel my_materials.py   # order of the props of mat_init_steel
```

The registry is validated at configure time (known solver, one model of every family the solver needs,
valid names). Changing a combination of models requires rebuilding; the parameters do not.

---

## Quick Usage Example

With the material `cp_vm_swift` of the default registry (von Mises, Swift hardening, closest point):

```fortran
program main
    use muscle_kinds, only : dp
    use muscle_status, only : STATUS_OK
    use muscle_mat_cp_vm_swift, only : mat_cp_vm_swift_dp, mat_init_cp_vm_swift, muscle_point_cp_vm_swift, &
        muscle_update_cp_vm_swift, MAT_CP_VM_SWIFT_NSTATE
    implicit none

    integer, parameter :: NS = MAT_CP_VM_SWIFT_NSTATE, NPTS = 1000
    type(mat_cp_vm_swift_dp) :: mat
    real(dp) :: strain(6), stress(6), ddsdde(6, 6), hsv_n(NS), hsv_np1(NS)
    real(dp), allocatable :: eps(:, :), hsv_bn(:, :), hsv_bnp1(:, :), sig(:, :)
    integer, allocatable :: st(:), it(:)
    integer :: status, iters

    ! young, poisson, Swift k, e0, n, Newton tolerance, iteration limit (muscle_meta.py props cp_vm_swift)
    call mat_init_cp_vm_swift(mat, [210000.0_dp, 0.3_dp, 500.0_dp, 0.01_dp, 0.2_dp, 1.0e-10_dp, 200.0_dp], status)
    if (status /= STATUS_OK) error stop 'invalid parameters'

    ! One point: total strain at t_n+1, tensorial Voigt order (xx, yy, zz, xy, yz, xz), no factor 2
    strain = [0.01_dp, -0.003_dp, -0.003_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    hsv_n = 0.0_dp                       ! virgin state (plastic strain, equivalent plastic strain)
    call muscle_point_cp_vm_swift(mat, strain, 0.0_dp, hsv_n, hsv_np1, stress, status, iters, ddsdde)
    if (status == STATUS_OK) print *, 'sigma_xx =', stress(1), ' iterations:', iters

    ! A batch of points in structure-of-arrays layout (OpenMP parallel loop)
    allocate(eps(NPTS, 6), hsv_bn(NPTS, NS), hsv_bnp1(NPTS, NS), sig(NPTS, 6), st(NPTS), it(NPTS))
    eps = spread(strain, 1, NPTS)
    hsv_bn = 0.0_dp
    call muscle_update_cp_vm_swift(mat, eps, 0.0_dp, hsv_bn, hsv_bnp1, sig, st, it)
end program main
```

The same calls exist in single precision with `mat_cp_vm_swift_sp` and `real(sp)` arrays. The FE code
(or the caller) keeps the converged history `hsv_n`: the update writes `hsv_np1` and never changes
`hsv_n`, so a rejected increment is simply retried from `hsv_n`.

---

## FEA Integration

MUSCLE ships no user subroutine: each material provides the pieces an adapter needs (`mat_init_<name>`
from the material constants, `muscle_point_<name>` for the point update and the tangent `ddsdde` with the
DDSDDE shear convention, `mat_rotate_hsv_<name>` for codes that do not rotate their state variables). The
FORD page *FEA Integration* (`pages/integration.md`) has the component orders of Abaqus, ANSYS and LS-DYNA
and skeletons of UMAT, USERMAT and umat41. Coming from the former object-oriented API (`Elasticity_linear`,
`Closest_point`, `Plastic_material_history`, ...)? See `pages/migration.md`.

---

## Licensing & Commercial Use

This library is licensed under the **GNU General Public License v3.0 (GPL-3.0-or-later)**. You are free to use, modify, and redistribute this software under the terms of this license.

### Commercial Licensing Notice
> **Important Note:** If your organization wishes to integrate MUSCLE (or parts of its solver engine) into proprietary commercial software, or use it under licensing terms different from the GPLv3 (to avoid releasing your proprietary source code), please contact the author to discuss a commercial license agreement.

---

## Author & Contact

*   **Author:** Matias Pacheco-Alarcon
*   **Email:** [matias.pacheco.a@gmail.com](mailto:matias.pacheco.a@gmail.com)

*Feel free to reach out for collaborations, academic inquiries, or support with FEA integration.*
