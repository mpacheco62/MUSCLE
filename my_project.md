---
project: MUSCLE
summary: Mechanics User-defined Solids Constitutive Library & Extensions. A high-performance Fortran library of constitutive material models for FEA (Abaqus, ANSYS APDL, LS-DYNA), in single and double precision, CPU and GPU.
src_dir: ./src
output_dir: ./doc
author: Matias Pacheco-Alarcon
email: matias.pacheco.a@usach.cl
fpp_extensions: fpp
                fypp
preprocess: true
preprocessor: fypp -M tools/fypp -m muscle_meta
include: ./include
macro: ENABLE_UDTIO
display: public
         protected
source: false
graph: true
search: true
license: by-nc
max_frontpage_items: 4
page_dir: pages
md_extensions: markdown.extensions.toc
---

# Introduction

**MUSCLE** is a specialized Fortran library designed for the development and implementation of advanced constitutive laws in Computational Solid Mechanics. It provides a robust, modular, and high-performance framework to define material behaviors that can be integrated into Finite Element Analysis (FEA) codes such as **LS-DYNA (User Defined Materials)** and **ANSYS APDL (USERMAT)**.

The library is written in Fortran 2008 as [FYPP](https://fypp.readthedocs.io) templates: tensors, models and solvers are generated from a single source for single (`sp`) and double (`dp`) precision (and, for the tensors and elasticity, for 3D and 2D), with static dispatch only (no polymorphism, `pure` procedures), so that the same code runs on CPUs (OpenMP) and can be compiled for GPUs (OpenMP offload).

Models are combined into **declared materials**: a Python registry (`tools/fypp/muscle_materials.py`, or your own file with `-DMUSCLE_MATERIALS=<path>`) lists the combinations of elasticity, yield criterion, hardening or viscoplastic law, kinematics and solver, and the build generates one module `muscle_mat_<name>` per material. Users coming from the former object-oriented API: see the [migration guide](page/migration.html).

# Key Features

*   **Tensor Engine:** Optimized 2nd and 4th order tensor algebra using Voigt notation, in 3D and 2D, single and double precision, with `pure` operators usable in GPU kernels.
*   **Advanced Yield Criteria:** Support for Von Mises, Hill48, and the asymmetric CPB06 (Cazacu-Plunkett-Barlat) model.
*   **Plasticity Solvers:** Closest point projection (Newton-Raphson with consistent tangent; overstress and consistency viscoplasticity), monolithic closest point and cutting plane.
*   **Finite Strains:** total and updated Lagrangian, logarithmic (material and spatial) and corotational kinematics.
*   **Declared Materials:** per material a point API for FE user subroutines and a batch API over arrays of points (OpenMP, GPU-ready).
*   **Viscoplasticity:** A suite of high-strain-rate models including Johnson-Cook, Rusinek-Klepaczko, and Nemat-Nasser-Li.
*   **Performance:** Manual loop unrolling for tensor contractions and integration with Intel MKL/LAPACK for linear system solutions.

# Mathematical Foundation

The library solves the return mapping problem by minimizing the distance between a trial stress state and the yield surface. The flow rule is defined as:

\begin{align} 
    \pmb{\dot\varepsilon}^p &= \gamma \frac{\partial f\left(\pmb{\sigma}, \pmb{q}\right)}{\partial\pmb{\sigma}} \\
    \pmb{\dot \alpha} &= \gamma \frac{\partial f\left(\pmb{\sigma}, \pmb{q}\right)}{\partial\pmb{q}} 
\end{align}

Where \(f\) represents the yield function and \(\gamma\) the plastic multiplier.

# Architecture Overview

The project is organized into several functional layers:

### 1. Tensor Types (`muscle_tensors`)
The foundation of the library: symmetric and general 2nd and 4th order tensors and the identities, all available through `use muscle_tensors` (see the [Tensor Algebra Guide](page/tensors.html)).
- One type per dimension and kind, `ten_<dim><family>_<kind>`: e.g. `ten_3D2Osym_dp` (2nd order symmetric, 6 components), `ten_2D4O3sym_sp` (4th order fully symmetric, 10 components). The historical names without suffix (`ten_3D2Osym`, ...) are the double precision types.
- Identities: the unscaled ones (`iden_2O`, `iden_4O3T`, `iden_4O4T`) have no kind and serve both precisions; the scaled ones have one type per kind (`iden_2OS_sp`, `iden_2OS_dp`, ...).
- Operators generic over dimension and kind, also between different types: `+`, `-`, `*`, `/`, `=`, `.ddot.` (double contraction), `.tdot.` and `.tdotsym.` (dyadic products), `.transform.`, `.dev.`, `.inv.`. Both operands must have the same kind (no mixed `sp`/`dp` operators).
- Free `pure` functions for kernels (`norm`, `det`, `inverse`, ...) in `muscle_tensor_functions`; the type-bound procedures (`%init()`, `%xx()`, `%norm()`, ...) are kept for host code. Formatted output (`write`) is available only when built with `ENABLE_UDTIO`.

### 2. Model families
Each model is a concrete parameter type (`<family>_<model>_<k>`, e.g. `hard_swift_dp`) with free `pure` functions; the family module gives generics over every model and kind (`hard_stress`, `yield_stress_eq`, `elas_stress`, ...), resolved at compile time.
- **Elasticity** (`muscle_elas_laws`): linear isotropic, and the frozen viscoelastic elasticity of a step.
- **Viscoelasticity** (`muscle_visco_laws`): generalized Maxwell, Kelvin-Voigt.
- **Kinematics** (`muscle_kin_laws`): small strain, total/updated Lagrangian, material/spatial logarithmic, corotational.

### 3. Yield Criteria (`muscle_yield_laws`)
Implements the surface that delimits the elastic-plastic transition.
- [[muscle_yield_vonmises]]: Classic J2 plasticity.
- [[muscle_yield_hill48]], [[muscle_yield_druckerprager]], [[muscle_yield_yld2004]].
- [[muscle_yield_cpb06]], [[muscle_yield_cpb06ex2]]: criteria for materials with strength differential effects (asymmetry between tension and compression).

### 4. Hardening and Viscoplastic Laws (`muscle_hard_laws`, `muscle_vp_laws`, `muscle_os_laws`)
Defines the evolution of the flow stress:
- **Bilinear, Swift / Ludwik, modified Voce:** isotropic hardening.
- **Johnson-Cook, Rusinek-Klepaczko, Nemat-Nasser-Li, Voyiadjis-Almasri:** rate-dependent flow stress.
- **Perzyna, Chaboche:** overstress viscoplasticity.

### 5. Solvers and Declared Materials (`src/solvers`, `src/materials`)
The return mapping algorithms are bodies (`src/solvers/muscle_<solver>_body.fypp`) included by the material template, which calls only the family generics:
- **`closest_point`:** implicit closest point projection with the consistent tangent (quadratic convergence of the global FEA Newton iterations); **`closest_point_perzyna`** and **`closest_point_time`**: overstress and consistency viscoplasticity; **`closest_point_2`**: monolithic 7x7 Newton; **`cutting_plane`**: explicit, continuum tangent.
- Each generated module `muscle_mat_<name>` provides the material type `mat_<name>_<k>`, `mat_init_<name>` (flat parameter array), the point API `muscle_point_<name>` for FE user subroutines (see [FEA Integration](page/integration.html)) and the batch API `muscle_update_<name>` over SoA arrays (`eps(npts, 6)`, optional index list).
- New models and solvers: see the [Developer Guide](page/contributing.html).

# Installation and Building

The project uses **CMake** and **FYPP**. To build the library and the test suite:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target build_all_tests
ctest --test-dir build
```

`-DMUSCLE_REAL_KINDS=dp` generates only double precision (faster development builds); `-DMUSCLE_MATERIALS=<path>` replaces the test registry of materials with yours.

# Requirements

*   Fortran 2008 compliant compiler (GFortran 9+, Intel OneAPI).
*   Intel MKL or compatible BLAS/LAPACK library.
*   CMake 3.20+.
*   FYPP 3.1+ and Python 3 (preprocessor that generates the modules, material registry).

@Note
The tensor engine and the elasticity support 3D and 2D (plane strain/stress, axisymmetric) states with the same code; yield criteria, solvers and declared materials are 3D only.

@Bug
Numerical derivatives (`tests/common/muscle_fd.fypp`) require careful step size (`eps`) selection for highly non-linear yield surfaces like CPB06 to avoid round-off errors.
```