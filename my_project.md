---
project: UMatLib
summary: A high-performance, object-oriented Fortran library for constitutive material modeling in FEA (LS-DYNA, ANSYS APDL).
src_dir: ./src
output_dir: ./doc
author: Matias Pacheco-Alarcon
email: matias.pacheco.a@usach.cl
fpp_extensions: fpp
preprocess: true
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

**UMatLib** is a specialized Fortran library designed for the development and implementation of advanced constitutive laws in Computational Solid Mechanics. It provides a robust, modular, and high-performance framework to define material behaviors that can be integrated into Finite Element Analysis (FEA) codes such as **LS-DYNA (User Defined Materials)** and **ANSYS APDL (USERMAT)**.

The library is built on modern Fortran standards (2003/2008), utilizing Object-Oriented Programming (OOP) to allow easy extension of yield criteria, hardening laws, and viscoplastic models.

# Key Features

*   **Tensor Engine:** Optimized 2nd and 4th order tensor algebra using Voigt notation.
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

### 1. Tensor Types (`tensors_types`)
The foundation of the library. It defines symmetric and general tensors.
- `ten_3D2Osym`: 2nd order symmetric tensor (6 components).
- `ten_3D4O3sym`: 4th order fully symmetric tensor (21 components).
- Custom operators like `.ddot.` (double contraction) and `.tdotsym.` (symmetrized dyadic product).

### 2. Elasticity (`elasticity`)
Defines the elastic predictor. The base class `Base_elasticity` allows for both linear isotropic and transient (viscoelastic) behaviors.

### 3. Yield Criteria (`yield_criteria`)
Implements the surface that delimits the elastic-plastic transition. 
- [[mod_vonMises]]: Classic J2 plasticity.
- [[mod_CPB06]]: Advanced criteria for materials with strength differential effects (asymmetry between tension and compression).

### 4. Hardening Laws (`hardening_laws`)
Defines the evolution of the flow stress:
- **Swift / Ludwik:** Power law hardening.
- **Voce:** Exponential saturation hardening.

### 5. Return Mapping Solvers (`return_mapping`)
Numerical algorithms to project the trial stress back to the yield surface.
- **Closest Point Projection:** The default robust implicit solver.
- **Consistent Tangent:** Essential for maintaining quadratic convergence in the global FEA Newton iterations.

# Installation and Building

The project uses **CMake**. To build the library and the test suite:

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
ctest
```

# Requirements

*   Fortran 2008 compliant compiler (GFortran 9+, Intel OneAPI).
*   Intel MKL or compatible BLAS/LAPACK library.
*   CMake 3.10+.

@Note
This library is currently optimized for 3D stress states but includes support for 2D plane strain/stress formulations in specific modules.

@Bug
Numerical derivatives in the `muscle_math_derivatives` module require careful step size (`eps`) selection for highly non-linear yield surfaces like CPB06 to avoid round-off errors.
```