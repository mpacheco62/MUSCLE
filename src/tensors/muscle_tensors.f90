! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensors
    !! author: MPacheco
    !! version: 1.0 - Initial documentation
    !!
    !! Tensors Types Module
    !! ====================
    !!
    !! This module serves as the central access point for various tensor types.
    !! 
    !! This module serves as the central access point for various tensor types and their associated
    !! operations within the MUSCLE library. It aggregates tensor definitions and operator
    !! implementations from specialized submodules, providing a unified interface for users.
    !!
    !! Instead of defining the tensor types and operations directly, this module uses other
    !! modules (`muscle_tensor_*`) and makes specific types and
    !! operators publicly available. This promotes modularity and organization within the library.
    !!
    !! Nomenclature
    !! ------------
    !!
    !! The module exposes two main categories of tensor types: Regular Tensors and Identity Tensors.
    !!
    !! ### Regular Tensors:
    !!
    !! Nomenclature `ten_xyz`:
    !!
    !! - `x`: Dimensions of the tensor (e.g., 3D)
    !! - `y`: Order of the tensor (e.g., 2O, 4O)
    !! - `z`: Number of symmetries of the tensor (e.g., sym, 2sym, 3sym)
    !!
    !! Examples:
    !!
    !! - `ten_3D2O`: 3D tensor, second order, no specific symmetry assumed in the type name (likely stored as 3x3).
    !! - `ten_3D2Osym`: 3D tensor, second order, symmetric (likely stored efficiently, e.g., 6 components).
    !! - `ten_3D4O2sym`: 3D tensor, fourth order, with minor symmetries (\(C_{ijkl} = C_{jikl} = C_{ijlk}\)).
    !! - `ten_3D4O3sym`: 3D tensor, fourth order, with major and minor symmetries (\(C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij}\)).
    !!
    !! ### Identity Tensors:
    !!
    !! Nomenclature `iden_xyz` or `iden_xyzS`:
    !!
    !! - `x`: Dimensions of the tensor (e.g., 3D)
    !! - `y`: Order of the tensor (e.g., 2O, 4O)
    !! - `z`: Type of identity (e.g., None, 1T, 2T, 3T, 4T)
    !!     * `None`: Standard second-order identity tensor \(\delta_{ij}\). (`iden_2O`)
    !!     * `1T`: \(\delta_{ik}\delta_{jl}\), only valid for a fourth order tensor.
    !!     * `2T`: \(\delta_{il}\delta_{jk}\), only valid for a fourth order tensor.
    !!     * `3T`: \(\delta_{ij}\delta_{kl}\), only valid for a fourth order tensor. (`iden_4O3T`)
    !!     * `4T`: \(\frac{\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}}{2}\), symmetric fourth-order identity. (`iden_4O4T`)
    !! - `S`: Suffix indicating a *scaled* identity tensor (e.g., `val * identity`). (`iden_2OS`, `iden_4O3TS`, `iden_4O4TS`)
    !!
    !! Examples:
    !!
    !! - `iden_2O`: 3D identity tensor, second order (\(\delta_{ij}\)).
    !! - `iden_2OS`: Scaled 3D identity tensor, second order (\(c \delta_{ij}\)).
    !! - `iden_4O3T`: 3D identity tensor, fourth order, type 3 (\(\delta_{ij}\delta_{kl}\)).
    !! - `iden_4O4T`: 3D identity tensor, fourth order, type 4 (symmetric identity).
    !! - `iden_4O4TS`: Scaled 3D identity tensor, fourth order, type 4 (\(c \cdot (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})/2\)).
    !!
    !! Public Entities
    !! ---------------
    !!
    !! This module makes the following entities publicly available:
    !!
    !! ### Derived Types:
    !!
    !! - [[ten_3D2O]]:       General 3D second-order tensor (from [[muscle_tensor_3d2o]]).
    !! - [[ten_3D2Osym]]:    Symmetric 3D second-order tensor (from [[muscle_tensor_3d2osym]]).
    !! - [[ten_3D4O2sym]]:   3D fourth-order tensor with minor symmetries (from [[muscle_tensor_3d4o2sym]]).
    !! - [[ten_3D4O3sym]]:   3D fourth-order tensor with major and minor symmetries (from [[muscle_tensor_3d4o3sym]]).
    !! - [[iden_2O]]:      Standard 3D second-order identity tensor (from [[muscle_tensor_iden_2o]]).
    !! - [[iden_2OS]]:   Scaled 3D second-order identity tensor (from [[muscle_tensor_iden_2os]]).
    !! - [[iden_4O3T]]:    Type 3 3D fourth-order identity tensor (from [[muscle_tensor_iden_4o3t]]).
    !! - [[iden_4O3TS]]: Scaled Type 3 3D fourth-order identity tensor (from [[muscle_tensor_iden_4o3ts]]).
    !! - [[iden_4O4T]]:    Type 4 (symmetric) 3D fourth-order identity tensor (from [[muscle_tensor_iden_4o4t]]).
    !! - [[iden_4O4TS]]: Scaled Type 4 3D fourth-order identity tensor (from [[muscle_tensor_iden_4o4ts]]).
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded for various combinations of the public tensor types
    !! and intrinsic types (like `real(real64)`). The implementations are provided by the
    !! `mod_operator_*` modules.
    !!
    !! - `+`: Addition (e.g., `tensor + tensor`).
    !! - `-`: Subtraction (e.g., `tensor - tensor`).
    !! - `*`: Multiplication (scalar * tensor, tensor * scalar, potentially tensor * tensor depending on definitions).
    !! - `/`: Division (tensor / scalar).
    !! - `.approx.`: Custom equality comparison for tensors (checks for approximate equality of components).
    !! - `.dev.`: Deviatoric part of a second-order tensor.
    !! - `.ddot.`: Double dot product (e.g., `tensor4 : tensor2`, `tensor2 : tensor2`).
    !! - `.tdot.`: Tensor dot product (specific definition depends on implementation, often \( (A \otimes B)_{ijkl} = A_{ij} B_{kl} \) for second order).
    !! - `.inv.`: Inverse of a 3D fourth-order tensor, with major and minor symmetries.
    !!
    !! ### Assignment:
    !!
    !! - `=`: Overloaded assignment allows copying between compatible tensor types.
    !!
    !! Usage
    !! -----
    !!
    !! To use the tensor types and operations defined here, simply add `use muscle_tensors`
    !! to your Fortran code.
    !!
    !! ```fortran
    !! program example_usage
    !!   use muscle_tensors
    !!   implicit none
    !!
    !!   type(ten_3D2Osym) :: stress, strain, stress_dev
    !!   type(ten_3D4O3sym) :: C_elastic
    !!   real(real64) :: scalar_val
    !!
    !!   ! ... initialize tensors ...
    !!
    !!   stress = C_elastic .ddot. strain  ! Double dot product via overloaded operator
    !!   stress_dev = .dev. stress   ! Deviatoric part
    !!   stress = stress + stress_dev * scalar_val ! Addition, multiplication
    !!
    !!   ! ... etc ...
    !!
    !! end program example_usage
    !! ```

    use, intrinsic :: iso_fortran_env
    use muscle_tensor_3d2o
    use muscle_tensor_3d2osym
    use muscle_tensor_2d2osym
    use muscle_tensor_3d4o3sym
    use muscle_tensor_2d4o3sym
    use muscle_tensor_3d4o2sym
    use muscle_tensor_iden_2os
    use muscle_tensor_iden_2o
    use muscle_tensor_iden_4o3ts
    use muscle_tensor_iden_4o3t
    use muscle_tensor_iden_4o4ts
    use muscle_tensor_iden_4o4t
    use muscle_tensor_ops_addition_subtraction
    use muscle_tensor_ops_contraction_double
    use muscle_tensor_ops_dyadic
    use muscle_tensor_ops_contraction_single
    use muscle_tensor_ops_assignment
    use muscle_tensor_ops_transform
    implicit None

    public :: ten_3D2O
    public :: ten_3D2Osym
    public :: ten_2D2Osym
    public :: ten_3D4O3sym
    public :: ten_2D4O3sym
    public :: ten_3D4O2sym
    public :: iden_2OS
    public :: iden_2O
    public :: iden_4O3TS
    public :: iden_4O3T
    public :: iden_4O4TS
    public :: iden_4O4T

    private 

    public :: operator(.approx.)
    public :: operator(+)
    public :: operator(-)
    public :: operator(*)
    public :: operator( / )
    
    public :: operator(.dev.)

    public :: operator(.ddot.)

    public :: operator(.tdot.)
    public :: operator(.tdotsym.)
    public :: operator(.inv.)
    public :: operator(.transform.)

    public :: transpose
    

    public :: assignment (=)

#ifdef ENABLE_UDTIO
    public :: write(formatted)
#endif
end module muscle_tensors