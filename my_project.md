---
project: MUSCLE
summary: Mechanics User-defined Solids Constitutive Library & Extensions. A high-performance, object-oriented Fortran library for constitutive material modeling in FEA (LS-DYNA, ANSYS APDL).
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
---

# Introduction

**MUSCLE** is a specialized Fortran library designed for the development and implementation of advanced constitutive laws in Computational Solid Mechanics. It provides a robust, modular, and high-performance framework to define material behaviors that can be integrated into Finite Element Analysis (FEA) codes such as **LS-DYNA (User Defined Materials)** and **ANSYS APDL (USERMAT)**.

The library is built on modern Fortran standards (2003/2008), utilizing Object-Oriented Programming (OOP) to allow easy extension of yield criteria, hardening laws, and viscoplastic models.

The tensor engine is generated with [FYPP](https://fypp.readthedocs.io) from a single source for 3D and 2D and for single (`sp`) and double (`dp`) precision, with static dispatch only (no polymorphism in the operations), so that it can also be compiled for GPUs. The constitutive models are being migrated to the same scheme; for now they are double precision only.

# Key Features

*   **Tensor Engine:** Optimized 2nd and 4th order tensor algebra using Voigt notation, in 3D and 2D, single and double precision, with `pure` operators usable in GPU kernels.
*   **Advanced Yield Criteria:** Support for Von Mises, Hill48, and the asymmetric CPB06 (Cazacu-Plunkett-Barlat) model.
*   **Plasticity Solvers:** Robust implementation of the **Closest Point Projection** (Implicit Newton-Raphson) algorithm with Consistent Tangent Operator calculation.
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

### 2. Elasticity (`elasticity`)
Defines the elastic predictor. The base class `Base_elasticity` allows for both linear isotropic and transient (viscoelastic) behaviors.

### 3. Yield Criteria (`yield_criteria`)
Implements the surface that delimits the elastic-plastic transition. 
- [[muscle_yield_vonmises]]: Classic J2 plasticity.
- [[muscle_yield_cpb06]]: Advanced criteria for materials with strength differential effects (asymmetry between tension and compression).

### 4. Hardening Laws (`hardening_laws`)
Defines the evolution of the flow stress:
- **Swift / Ludwik:** Power law hardening.
- **Voce:** Exponential saturation hardening.

### 5. Return Mapping Solvers (`return_mapping`)
Numerical algorithms to project the trial stress back to the yield surface.
- **Closest Point Projection:** The default robust implicit solver.
- **Consistent Tangent:** Essential for maintaining quadratic convergence in the global FEA Newton iterations.

# Installation and Building

The project uses **CMake** and **FYPP**. To build the library and the test suite:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target build_all_tests
ctest --test-dir build
```

`-DMUSCLE_REAL_KINDS=dp` generates only the double precision tensors (faster development builds).

# Requirements

*   Fortran 2008 compliant compiler (GFortran 9+, Intel OneAPI).
*   Intel MKL or compatible BLAS/LAPACK library.
*   CMake 3.20+.
*   FYPP 3.1+ (Python preprocessor that generates the tensor modules).

@Note
The tensor engine supports 3D and 2D (plane strain/stress, axisymmetric) states with the same code; the constitutive models are currently optimized for 3D stress states and support 2D formulations only in specific modules.

@Bug
Numerical derivatives in the `muscle_math_derivatives` module require careful step size (`eps`) selection for highly non-linear yield surfaces like CPB06 to avoid round-off errors.
```