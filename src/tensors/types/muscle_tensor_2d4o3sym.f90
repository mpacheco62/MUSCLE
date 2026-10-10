! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_2d4o3sym
    !! Module muscle_tensor_2d4o3sym
    !! =============================
    !!
    !! Defines the type for fully symmetric 2D fourth-order tensors and associated operations.
    !!
    !! This module provides the derived type `ten_2D4O3sym` to represent a fourth-order
    !! tensor in a 2D setting (plane or axisymmetric, with the out-of-plane normal
    !! direction kept) possessing both major and minor symmetries:
    !! \[ C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij} \]
    !! Such tensors are common in linear elasticity (the elasticity tensor).
    !!
    !! The tensor is stored internally using a compressed format with 10 components,
    !! corresponding to the independent elements of the symmetric 4x4 Voigt matrix representation.
    !! The storage order follows a specific convention (see type description).
    !!
    !! The module overloads standard arithmetic operators (+, -, *, /), a custom
    !! equality comparison operator (.approx.), an inverse operator (.inv.),
    !! and provides methods for initialization.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ten_2D4O3sym`: Represents a 2D fourth-order tensor with major and minor symmetries.
    !!     - Component: `vals(10) :: real(real64)` - Stores the 10 independent components.
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       10-element array or from 10 individual components.
    !!     - Procedures: `norm`, `set`, `is_approx`.
    !!
    !! ### Operators:
    !!
    !! - `.approx.`: Compares two `ten_2D4O3sym` tensors for approximate equality.
    !! - `.inv.`: Computes the inverse of the tensor based on its 4x4 symmetric Voigt matrix representation.
    !! - `+`: Adds two `ten_2D4O3sym` tensors.
    !! - `-`: Subtracts two `ten_2D4O3sym` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_2D4O3sym` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_2D4O3sym` tensor by a `real(real64)` scalar.
    !!
    !! ### Assignment:
    !!
    !! - `=`: Assigns a single `real(real64)` scalar to all 10 components.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_2D4O3sym_usage
    !!   use muscle_tensors ! Includes muscle_tensor_2d4o3sym and others
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_2D4O3sym) :: C_iso, C_inv
    !!   real(real64) :: lambda, mu
    !!   logical :: are_equal
    !!
    !!   ! Initialize isotropic elasticity tensor using individual components
    !!   lambda = 120.0D0 ! First Lame parameter
    !!   mu = 80.0D0      ! Shear modulus (second Lame parameter)
    !!
    !!   ! Components based on init2 order:
    !!   ! xxxx, yyyy, zzzz, xyxy, xxyy, yyzz, zzxy,
    !!   ! xxzz, yyxy, xxxy
    !!   call C_iso%init( &
    !!       lambda + 2*mu, lambda + 2*mu, lambda + 2*mu,  & ! xxxx, yyyy, zzzz
    !!       mu,                                           & ! xyxy
    !!       lambda, lambda,                               & ! xxyy, yyzz
    !!       0.0D0,                                        & ! zzxy
    !!       lambda,                                       & ! xxzz
    !!       0.0D0,                                        & ! yyxy
    !!       0.0D0                                        )  ! xxxy
    !!
    !!   ! Calculate the inverse (compliance tensor S = C^-1)
    !!   C_inv = .inv. C_iso
    !!
    !!   ! Check a component of the inverse
    !!   ! S_1111 = C_inv%vals(1) should be (lambda+mu)/(mu*(3*lambda+2*mu))
    !!   print *, "C_iso(1) (C_1111):", C_iso%vals(1)
    !!   print *, "C_inv(1) (S_1111):", C_inv%vals(1)
    !!   print *, "Expected S_1111:", (lambda+mu)/(mu*(3*lambda+2*mu))
    !!
    !!   ! Comparison
    !!   are_equal = (C_iso .approx. C_iso)
    !!   print *, "Is C_iso equal to itself?", are_equal
    !!
    !! end program example_ten_2D4O3sym_usage
    !! ```
    !!
    !! For more information see [[muscle_tensors]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: ten_2D4O3sym
        !! Fully Symmetric 2D Fourth-Order Tensor (10-Component Storage)
        !! ==============================================================
        !!
        !! Represents a fourth-order tensor \(C_{ijkl}\) in a 2D setting that possesses
        !! both minor and major symmetries:
        !! \[ C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij} \]
        !! This is the standard symmetry for linear elastic constitutive tensors.
        !!
        !! Storage:
        !! --------
        !! Due to the symmetries, only 10 independent components exist. These are stored
        !! internally in a 1D array `vals` of size 10. The storage order corresponds
        !! to the upper triangle of the 4x4 symmetric Voigt matrix representation,
        !! stored by diagonals:
        !!
        !! Voigt Matrix (Indices IJ):
        !! ```
        !! | 11 12 13 14 |
        !! |    22 23 24 |
        !! |       33 34 |
        !! |          44 |
        !! ```
        !! Storage order in `vals(1:10)`:
        !! (11, 22, 33, 44, 12, 23, 34, 13, 24, 14)
        !!
        !! where the Voigt index mapping is:
        !! - 1 <-> (1,1) or xx
        !! - 2 <-> (2,2) or yy
        !! - 3 <-> (3,3) or zz
        !! - 4 <-> (1,2) or (2,1) or xy
        !!
        !! Access:
        !! -------
        !! Components are typically accessed directly via the `vals` array using the
        !! appropriate index based on the storage order, or modified with `set`.
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 10-element `real(real64)`
        !! array (following the storage order above) or by providing the 10 components
        !! individually (see `init2_ten_2D4O3sym` for the required input order).
        !!
        !! For more information see [[muscle_tensors]]

        real(real64), dimension(10) :: vals
            !! Stores the 10 independent components following the compressed Voigt storage order.
        contains
            generic, public :: init => init_ten_2D4O3sym, init2_ten_2D4O3sym
                !! Generic interface for initialization.
            procedure, private :: init_ten_2D4O3sym, init2_ten_2D4O3sym
            procedure, public :: norm => norm_2D4O3sym
                !! Modified L1 norm of the tensor.
            procedure, public :: set => set_ten2D4O3sym
                !! Sets selected components by name.
            procedure, public :: is_approx => is_approx_2D4O3sym
                !! Compares two tensors for approximate equality.
    end type ten_2D4O3sym

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_2D4O3sym
    end interface

    public :: operator(.inv.)
    interface operator (.inv.)
        module procedure inv_2D4O3sym
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_2D4O3sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D4O3sym
        module procedure subU_2D4O3sym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_2D4O3sym
        module procedure mul_2D4O3sym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_2D4O3sym_real64
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_ten_2D4O3sym_real64
    end interface
contains

    pure subroutine init_ten_2D4O3sym(self, vals)
        !! Initializes a ten_2D4O3sym tensor from a 10-element array (compressed Voigt order).
        implicit none
        class(ten_2D4O3sym), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in) :: vals(10)
            !! Components in the storage order (11, 22, 33, 44, 12, 23, 34, 13, 24, 14).
        self%vals = vals
    end subroutine init_ten_2D4O3sym

    pure subroutine init2_ten_2D4O3sym(self,             &
                                              xxxx, yyyy, zzzz, &
                                              xyxy, xxyy, yyzz, &
                                              zzxy, xxzz, yyxy, &
                                              xxxy              &
                                              )
        !! Initializes a ten_2D4O3sym tensor from its 10 individual components.
        !! Input arguments correspond to the independent Voigt matrix components C(I,J):
        !!
        !! ```
        !!  | (1,1) (1,2) (1,3) (1,4) |   <- xxxx, xxyy, xxzz, xxxy
        !!  | (2,1) (2,2) (2,3) (2,4) |   <- yyxx, yyyy, yyzz, yyxy
        !!  | (3,1) (3,2) (3,3) (3,4) |   <- zzxx, zzyy, zzzz, zzxy
        !!  | (4,1) (4,2) (4,3) (4,4) |   <- xyxx, xyyy, xyzz, xyxy
        !! ```
        implicit none
        class(ten_2D4O3sym), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in) :: xxxx
            !! \(C_{1111}\)
        real(real64), intent(in) :: yyyy
            !! \(C_{2222}\)
        real(real64), intent(in) :: zzzz
            !! \(C_{3333}\)
        real(real64), intent(in) :: xyxy
            !! \(C_{1212}\)
        real(real64), intent(in) :: xxyy
            !! \(C_{1122}\)
        real(real64), intent(in) :: yyzz
            !! \(C_{2233}\)
        real(real64), intent(in) :: zzxy
            !! \(C_{3312}\)
        real(real64), intent(in) :: xxzz
            !! \(C_{1133}\)
        real(real64), intent(in) :: yyxy
            !! \(C_{2212}\)
        real(real64), intent(in) :: xxxy
            !! \(C_{1112}\)
        self%vals = (/xxxx, yyyy, zzzz, xyxy, &
                      xxyy, yyzz, zzxy, &
                      xxzz, yyxy, &
                      xxxy &
                      /)
    end subroutine

    pure function is_approx_2D4O3sym(a, b, tol) result(res)
        !! Component-wise approximate equality with a relative tolerance:
        !! \( |a_I - b_I| \le tol \cdot \max(\max|a|, \max|b|, 10^{-30}) \) for all I.
        implicit none
        class(ten_2D4O3sym), intent(in) :: a
            !! First tensor.
        class(ten_2D4O3sym), intent(in) :: b
            !! Second tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default `1E-12`).
        logical :: res
            !! `.TRUE.` if both tensors are approximately equal.
        real(real64), parameter :: EPS_ABS=1e-30
        real(real64) :: eps_check
        real(real64) :: tol2
        real(real64) :: max_val

        tol2 = 1E-12
        if (present(tol)) tol2 = tol
        max_val = max(maxval(abs(a%vals)), maxval(abs(b%vals)), EPS_ABS)
        eps_check = tol2 * max_val

        res = all(abs(a%vals - b%vals) .le. eps_check)
    end function is_approx_2D4O3sym

    pure function approx_2D4O3sym(a, b) result(res)
        !! `.approx.` Compares two ten_2D4O3sym tensors for approximate equality.
        !! Delegates to `is_approx` with its default tolerance.
        implicit none
        type(ten_2D4O3sym), intent(in) :: a
            !! First tensor.
        type(ten_2D4O3sym), intent(in) :: b
            !! Second tensor.
        logical :: res
            !! `.TRUE.` if both tensors are approximately equal.

        res = a%is_approx(b)
    end function approx_2D4O3sym

    pure function norm_2D4O3sym(a) result(norms)
        !! Modified L1 norm based on the 10 stored components, where components
        !! corresponding to off-diagonal Voigt matrix entries are weighted by 2:
        !! norm(A) = sum(|A_diag|) + 2*sum(|A_offdiag|) based on the 4x4 Voigt matrix.
        implicit none
        class(ten_2D4O3sym), intent(in) :: a
            !! Tensor.
        real(real64) :: norms
            !! Norm.
        norms = sum(abs(a%vals(1:4))) + 2.0D0 * sum(abs(a%vals(5:10)))
    end function norm_2D4O3sym

    pure function sum_2D4O3sym(a, b) result(res)
        !! Component-wise sum.
        implicit none
        type(ten_2D4O3sym), intent(in) :: a
            !! First operand.
        type(ten_2D4O3sym), intent(in) :: b
            !! Second operand.
        type(ten_2D4O3sym) :: res
            !! Sum.
        res%vals = a%vals + b%vals
    end function sum_2D4O3sym

    pure function sub_2D4O3sym(a, b) result(res)
        !! Component-wise difference.
        implicit none
        type(ten_2D4O3sym), intent(in) :: a
            !! Minuend.
        type(ten_2D4O3sym), intent(in) :: b
            !! Subtrahend.
        type(ten_2D4O3sym) :: res
            !! Difference.
        res%vals = a%vals - b%vals
    end function sub_2D4O3sym

    pure function subU_2D4O3sym(a) result(res)
        !! Unary negation.
        implicit none
        type(ten_2D4O3sym), intent(in) :: a
            !! Operand.
        type(ten_2D4O3sym) :: res
            !! Negated tensor.
        res%vals = -a%vals
    end function subU_2D4O3sym

    pure function mul_real64_2D4O3sym(a, b) result(res)
        !! Scalar product \(a\,\mathbb{B}\).
        implicit none
        real(real64), intent(in) :: a
            !! Scalar factor.
        type(ten_2D4O3sym), intent(in) :: b
            !! Tensor.
        type(ten_2D4O3sym) :: res
            !! Scaled tensor.
        res%vals = a * b%vals
    end function mul_real64_2D4O3sym

    pure function mul_2D4O3sym_real64(b, a) result(res)
        !! Scalar product \(\mathbb{B}\,a\).
        implicit none
        type(ten_2D4O3sym), intent(in) :: b
            !! Tensor.
        real(real64), intent(in) :: a
            !! Scalar factor.
        type(ten_2D4O3sym) :: res
            !! Scaled tensor.
        res%vals = a * b%vals
    end function mul_2D4O3sym_real64

    pure function div_2D4O3sym_real64(b, a) result(res)
        !! Scalar division \(\mathbb{B}/a\).
        implicit none
        type(ten_2D4O3sym), intent(in) :: b
            !! Tensor.
        real(real64), intent(in) :: a
            !! Scalar divisor (must be non-zero).
        type(ten_2D4O3sym) :: res
            !! Scaled tensor.
        res%vals = b%vals/a
    end function div_2D4O3sym_real64

    pure function inv_2D4O3sym(a) result(res)
        !! `.inv.` Computes the inverse of a ten_2D4O3sym tensor, i.e. the tensor
        !! \(\mathbb{A}^{-1}\) such that \(\mathbb{A}^{-1} : \mathbb{A} = \mathbb{I}^S\).
        !!
        !! Reconstructs the 4x4 symmetric Voigt matrix, inverts it with `M44INV`
        !! (LAPACK LU) and extracts the 10 independent components of the inverse,
        !! scaling the shear row/column by 1/2 and the shear-shear term by 1/4.
        !! This keeps the result consistent with the tensorial-shear Voigt convention
        !! used by `.ddot.` (same scheme as `inv_3D4O3sym`).
        !! A zero tensor is returned if the matrix is singular.
        use muscle_math_inverses, only : M44INV
        implicit none
        type(ten_2D4O3sym), intent(in) :: a
            !! Tensor to invert.
        type(ten_2D4O3sym) :: res
            !! Inverse tensor.
        real(real64) :: mat_a(4,4), mat_b(4,4), v(10)
        logical :: ok
        v = a%vals
        mat_a = reshape((/  v(1),  v(5),  v(8), v(10), &
                            v(5),  v(2),  v(6),  v(9), &
                            v(8),  v(6),  v(3),  v(7), &
                           v(10),  v(9),  v(7),  v(4)  &
                        /), (/4,4/))

        call M44INV(mat_a, mat_b, ok)

        v = (/ mat_b(1,1), mat_b(2,2), mat_b(3,3), mat_b(4,4)/4D0, &
               mat_b(1,2), mat_b(2,3), mat_b(3,4)/2D0,             &
               mat_b(1,3), mat_b(2,4)/2D0,                         &
               mat_b(1,4)/2D0                                      &
             /)

        call res%init(v)
    end function inv_2D4O3sym

    pure subroutine assign_ten_2D4O3sym_real64(a, b)
        !! Assigns the scalar `b` to every stored component of `a`.
        implicit none
        type(ten_2D4O3sym), intent(out) :: a
            !! Target tensor.
        real(real64), intent(in) :: b
            !! Scalar value.
        a%vals = b
    end subroutine assign_ten_2D4O3sym_real64

    pure subroutine set_ten2D4O3sym(a, &
                                    xxxx, yyyy, zzzz, xyxy, &
                                    xxyy, yyzz, zzxy, xxzz, &
                                    yyxy, xxxy              &
                                    )
        !! Sets only the components passed as (optional, named) arguments.
        implicit none
        class(ten_2D4O3sym), intent(inout) :: a
            !! Tensor to modify.
        real(real64), optional, intent(in) :: xxxx
            !! \(C_{1111}\)
        real(real64), optional, intent(in) :: yyyy
            !! \(C_{2222}\)
        real(real64), optional, intent(in) :: zzzz
            !! \(C_{3333}\)
        real(real64), optional, intent(in) :: xyxy
            !! \(C_{1212}\)
        real(real64), optional, intent(in) :: xxyy
            !! \(C_{1122}\)
        real(real64), optional, intent(in) :: yyzz
            !! \(C_{2233}\)
        real(real64), optional, intent(in) :: zzxy
            !! \(C_{3312}\)
        real(real64), optional, intent(in) :: xxzz
            !! \(C_{1133}\)
        real(real64), optional, intent(in) :: yyxy
            !! \(C_{2212}\)
        real(real64), optional, intent(in) :: xxxy
            !! \(C_{1112}\)

        if (present(xxxx)) a%vals(1) = xxxx
        if (present(yyyy)) a%vals(2) = yyyy
        if (present(zzzz)) a%vals(3) = zzzz
        if (present(xyxy)) a%vals(4) = xyxy
        if (present(xxyy)) a%vals(5) = xxyy
        if (present(yyzz)) a%vals(6) = yyzz
        if (present(zzxy)) a%vals(7) = zzxy
        if (present(xxzz)) a%vals(8) = xxzz
        if (present(yyxy)) a%vals(9) = yyxy
        if (present(xxxy)) a%vals(10) = xxxy
    end subroutine

end module muscle_tensor_2d4o3sym
