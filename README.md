# MUSCLE: Mechanics User-defined Solids Constitutive Library & Extensions

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Fortran](https://img.shields.io/badge/Language-Fortran_2008-orange.svg)](https://gcc.gnu.org/fortran/)

**MUSCLE** is a high-performance, object-oriented Fortran library (compliant with Fortran 2003/2008 standards) designed for the development, testing, and implementation of advanced constitutive material models in Computational Solid Mechanics. 

It provides a modular, robust, and extensible framework easily wrapped by standard user subroutines in Finite Element Analysis (FEA) solvers, such as **LS-DYNA (UMAT)** and **ANSYS APDL (USERMAT)**.

---

## Key Features

### 1. Optimized Tensor Engine (`muscle_tensors`)
A specialized algebraic engine using Voigt notation with manual loop unrolling for maximum performance:
*   **Symmetric 2nd Order Tensors** (`ten_3D2Osym`, `ten_2D2Osym`).
*   **General 2nd Order Tensors** (`ten_3D2O`).
*   **Fully Symmetric 4th Order Tensors** (`ten_3D4O3sym`, `ten_2D4O3sym`).
*   **Custom Overloaded Operators:** Double contraction (`.ddot.`), symmetrized dyadic product (`.tdotsym.`), deviatoric projector (`.dev.`), and analytical fourth-order inversion (`.inv.`).
*   **Identity Tensors:** Support for symbolic standard and scaled identities (`iden_2O`, `iden_4O4T`, etc.).

### 2. Advanced Yield Criteria (`yield_criteria`)
*   **Von Mises:** Classic $J_2$ plasticity with analytical first and second derivatives.
*   **Hill48:** Orthotropic yield criterion.
*   **CPB06 (Cazacu-Plunkett-Barlat):** Advanced orthotropic criterion for materials exhibiting Strength Differential (SD) effects (tension-compression asymmetry), featuring exact, eigenvector-free analytical derivatives.

### 3. Isotropic Hardening Laws (`hardening_laws`)
*   **Bilinear:** Linear hardening.
*   **Swift & Ludwik:** Power-law hardening models.
*   **Voce Modified:** Combined linear and exponential saturation hardening.

### 4. High-Strain-Rate Viscoplasticity (`viscoplastic_laws`)
*   **Johnson-Cook (JC / JC-Full):** Rate-dependence with optional thermal softening.
*   **Rusinek-Klepaczko (RK / MRK):** Standard and modified rate-dependent formulations.
*   **Nemat-Nasser-Li (NNL):** Microstructure-based thermally activated model with numerical stabilization.
*   **Voyiadjis-Almasri (VA):** Energy-based viscoplastic model.

### 5. Return-Mapping Solvers (`solvers`)
*   **Closest Point Projection (CPP):** Robust implicit Newton-Raphson return-mapping solver.
*   **Consistent Tangent Operator:** Analytical computation of the algorithmic tangent tensor, essential to preserve quadratic convergence in global FEA implicit iterations.
*   **Cutting Plane:** Alternative explicit/implicit projection scheme.

---

## Requirements

*   **Compiler:** Fortran 2008 compliant compiler (GCC GFortran 9+ or Intel oneAPI `ifx`/`ifort`).
*   **Linear Algebra:** BLAS & LAPACK (or Intel MKL).
*   **Build System:** CMake 3.20+.

---

## Installation & Building

To build the static library and compile the complete test suite:

```bash
# Clone the repository
git clone https://github.com/mpacheco62/muscle.git
cd muscle

# Create a build directory
mkdir build && cd build

# Configure the project (Release mode)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Compile
make

# Run the test suite
ctest --output-on-failure
```

---

## Quick Usage Example

The following minimal example demonstrates how to initialize a tensor, define a linear elastic material, and evaluate stress:

```fortran
program main
    use muscle_tensors
    use muscle_elasticity_linear
    implicit none

    type(Elasticity_linear) :: material
    type(ten_3D2Osym) :: strain, stress
    real(8) :: E, nu

    E = 210000.0D0  ! Young's Modulus [MPa]
    nu = 0.3D0      ! Poisson's Ratio

    ! Initialize material properties
    call material%set_parameters(young=E, poisson=nu)

    ! Define a trial strain tensor (Voigt: xx, yy, zz, xy, yz, xz)
    call strain%init((/ 0.001D0, -0.0003D0, -0.0003D0, 0.0D0, 0.0D0, 0.0D0 /))

    ! Compute stress using Hooke's Law
    stress = material%stress(strain)

    print *, "Evaluated Stress XX [MPa]:", stress%xx()
end program main
```

---

## FEA Integration

MUSCLE is structurally optimized to be integrated into commercial Finite Element codes. A typical implementation inside a user subroutine (like LS-DYNA `umat` or ANSYS `usermat`) involves:

1.  Mapping the historical state variables and strain increments to `ten_3D2Osym` variables.
2.  Calling the return-mapping solver: `call closest_point2(strain, elasticity, hardening, yield, stress, strain_pf, strain_p, status)`.
3.  Retrieving the consistent tangent matrix via `call solver%tangent(...)` and exporting it to the solver's Jacobian array.

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