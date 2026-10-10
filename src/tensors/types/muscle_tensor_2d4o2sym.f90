! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_2d4o2sym
    !! Module muscle_tensor_2d4o2sym
    !! =============================
    !!
    !! Defines the type for 2D fourth-order tensors with minor symmetries and associated operations.
    !!
    !! This module provides the derived type `ten_2D4O2sym` to represent a fourth-order
    !! tensor in a 2D setting (plane or axisymmetric, out-of-plane normal direction kept)
    !! possessing minor symmetries, i.e., \(C_{ijkl} = C_{jikl} = C_{ijlk}\).
    !! It is the 2D counterpart of `ten_3D4O2sym` and appears, e.g., as the result of
    !! the dyadic product of two different symmetric tensors or of the double
    !! contraction of two fully symmetric fourth-order tensors.
    !!
    !! The tensor is stored internally using a 4x4 matrix based on the 2D Voigt index mapping
    !! (11->1, 22->2, 33->3, 12->4). Due to only having minor symmetries, this 4x4 matrix is
    !! generally **not** symmetric (\(C_{IJ} \neq C_{JI}\)).
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ten_2D4O2sym`: Represents a 2D fourth-order tensor with minor symmetries.
    !!     - Component: `vals(4,4) :: real(real64)` - Stores the 16 independent components
    !!       in a 4x4 matrix using Voigt index mapping.
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       4x4 array or from 16 individual components.
    !!     - Procedures: `is_approx`, `norm`, `convert_3sym`.
    !!
    !! ### Operators:
    !!
    !! - `.approx.`: Compares two `ten_2D4O2sym` tensors for approximate equality.
    !! - `+`: Adds two `ten_2D4O2sym` tensors.
    !! - `-`: Subtracts two `ten_2D4O2sym` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_2D4O2sym` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_2D4O2sym` tensor by a `real(real64)` scalar.
    !! - `.inv.`: Inverse with respect to \(\mathbb{I}^S\) (4x4 Voigt inversion with shear scaling).
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_2d4o2sym_usage
    !!   use muscle_tensors
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_2D2Osym)  :: a, b
    !!   type(ten_2D4O2sym) :: C
    !!
    !!   call a%init(1.0D0, 2.0D0, 3.0D0, 4.0D0)
    !!   call b%init(5.0D0, 6.0D0, 7.0D0, 8.0D0)
    !!   C = a .tdot. b          ! C(I,J) = a_I * b_J
    !!   print *, C%vals(1, 4)   ! C_1112 = a_11 * b_12 = 8.0
    !!
    !! end program example_ten_2d4o2sym_usage
    !! ```
    !!
    !! For more information see [[muscle_tensors]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: ten_2D4O2sym
        !! 2D Fourth-Order Tensor with Minor Symmetries (4x4 Voigt Storage)
        !! ================================================================
        !!
        !! Represents a fourth-order tensor \(C_{ijkl}\) in a 2D setting that possesses
        !! minor symmetries:
        !! \[ C_{ijkl} = C_{jikl} = C_{ijlk} \]
        !! It does **not** imply major symmetry (\(C_{ijkl} = C_{klij}\)).
        !!
        !! Storage:
        !! --------
        !! The 16 components are stored internally in a 4x4 `real(real64)` matrix `vals`
        !! using the 2D Voigt index mapping convention:
        !! - 1 <-> (1,1) or xx
        !! - 2 <-> (2,2) or yy
        !! - 3 <-> (3,3) or zz
        !! - 4 <-> (1,2) or (2,1) or xy
        !!
        !! The component \(C_{ijkl}\) is mapped to `vals(I, J)`, where \(I\) corresponds to the
        !! pair \((ij)\) and \(J\) corresponds to the pair \((kl)\).
        !! For example, \(C_{1122}\) is stored in `vals(1, 2)`, and \(C_{3312}\) in `vals(3, 4)`.
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 4x4 `real(real64)`
        !! array (corresponding directly to the `vals` matrix) or by providing the 16 components
        !! individually, row-wise (see `init2_ten_2D4O2sym`).
        !!
        !! For more information see [[muscle_tensors]]
        real(real64), dimension(4,4) :: vals
            !! Stores the 16 components as a 4x4 matrix using 2D Voigt index mapping.
        contains
            generic, public :: init => init_ten_2D4O2sym, init2_ten_2D4O2sym
                !! Generic interface for initialization.
            procedure, private :: init_ten_2D4O2sym, init2_ten_2D4O2sym
            procedure, public :: is_approx => is_approx_2D4O2sym
                !! Compares two tensors for approximate equality.
            procedure, public :: norm => norm_2D4O2sym
                !! L1 norm of the 4x4 Voigt matrix.
            procedure, public :: convert_3sym
                !! Projects onto a fully symmetric `ten_2D4O3sym` (upper triangle).
    end type ten_2D4O2sym

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_2D4O2sym
    end interface

    public :: operator(.inv.)
    interface operator (.inv.)
        module procedure inv_2D4O2sym
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_2D4O2sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D4O2sym
        module procedure subU_2D4O2sym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_2D4O2sym
        module procedure mul_2D4O2sym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_2D4O2sym_real64
    end interface
contains

    pure subroutine init_ten_2D4O2sym(self, vals)
        !! Initializes a ten_2D4O2sym tensor from a 4x4 array (Voigt matrix).
        implicit none
        class(ten_2D4O2sym), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in) :: vals(4,4)
            !! 4x4 Voigt matrix \(C_{IJ}\).
        self%vals = vals
    end subroutine init_ten_2D4O2sym

    pure subroutine init2_ten_2D4O2sym(self,                   &
                                       xxxx, xxyy, xxzz, xxxy, &
                                       yyxx, yyyy, yyzz, yyxy, &
                                       zzxx, zzyy, zzzz, zzxy, &
                                       xyxx, xyyy, xyzz, xyxy  &
                                       )
        !! Initializes a ten_2D4O2sym tensor from its 16 individual components,
        !! given row-wise:
        !!
        !! ```
        !!  | (1,1) (1,2) (1,3) (1,4) |   <- xxxx, xxyy, xxzz, xxxy
        !!  | (2,1) (2,2) (2,3) (2,4) |   <- yyxx, yyyy, yyzz, yyxy
        !!  | (3,1) (3,2) (3,3) (3,4) |   <- zzxx, zzyy, zzzz, zzxy
        !!  | (4,1) (4,2) (4,3) (4,4) |   <- xyxx, xyyy, xyzz, xyxy
        !! ```
        implicit none
        class(ten_2D4O2sym), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in) :: xxxx, xxyy, xxzz, xxxy
            !! Row 1: \(C_{11kl}\).
        real(real64), intent(in) :: yyxx, yyyy, yyzz, yyxy
            !! Row 2: \(C_{22kl}\).
        real(real64), intent(in) :: zzxx, zzyy, zzzz, zzxy
            !! Row 3: \(C_{33kl}\).
        real(real64), intent(in) :: xyxx, xyyy, xyzz, xyxy
            !! Row 4: \(C_{12kl}\).
        self%vals(:,1) = (/ xxxx, yyxx, zzxx, xyxx /)
        self%vals(:,2) = (/ xxyy, yyyy, zzyy, xyyy /)
        self%vals(:,3) = (/ xxzz, yyzz, zzzz, xyzz /)
        self%vals(:,4) = (/ xxxy, yyxy, zzxy, xyxy /)
    end subroutine

    pure function is_approx_2D4O2sym(a, b, tol) result(res)
        !! Component-wise approximate equality with a relative tolerance:
        !! \( |a_{IJ} - b_{IJ}| \le tol \cdot \max(\max|a|, \max|b|, 10^{-30}) \).
        implicit none
        class(ten_2D4O2sym), intent(in) :: a
            !! First tensor.
        class(ten_2D4O2sym), intent(in) :: b
            !! Second tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default `1E-8`, as in `ten_3D4O2sym`).
        logical :: res
            !! `.TRUE.` if both tensors are approximately equal.
        real(real64), parameter :: EPS_ABS=1e-30
        real(real64) :: max_val, eps_check, tol2

        tol2 = 1E-8
        if (present(tol)) tol2 = tol
        max_val = max(maxval(abs(a%vals)), maxval(abs(b%vals)), EPS_ABS)
        eps_check = tol2 * max_val

        res = all(abs(a%vals - b%vals) .le. eps_check)
    end function is_approx_2D4O2sym

    pure function approx_2D4O2sym(a, b) result(res)
        !! `.approx.` Compares two ten_2D4O2sym tensors for approximate equality.
        !! Delegates to `is_approx` with its default tolerance.
        implicit none
        type(ten_2D4O2sym), intent(in) :: a
            !! First tensor.
        type(ten_2D4O2sym), intent(in) :: b
            !! Second tensor.
        logical :: res
            !! `.TRUE.` if both tensors are approximately equal.

        res = a%is_approx(b)
    end function approx_2D4O2sym

    pure function norm_2D4O2sym(a) result(norm)
        !! L1 norm of the 4x4 Voigt matrix: \(\sum_{I,J} |A_{IJ}|\).
        implicit none
        class(ten_2D4O2sym), intent(in) :: a
            !! Tensor.
        real(real64) :: norm
            !! Norm.
        norm = sum(abs(a%vals))
    end function norm_2D4O2sym

    pure function sum_2D4O2sym(a, b) result(res)
        !! Component-wise sum.
        implicit none
        type(ten_2D4O2sym), intent(in) :: a
            !! First operand.
        type(ten_2D4O2sym), intent(in) :: b
            !! Second operand.
        type(ten_2D4O2sym) :: res
            !! Sum.
        res%vals = a%vals + b%vals
    end function sum_2D4O2sym

    pure function sub_2D4O2sym(a, b) result(res)
        !! Component-wise difference.
        implicit none
        type(ten_2D4O2sym), intent(in) :: a
            !! Minuend.
        type(ten_2D4O2sym), intent(in) :: b
            !! Subtrahend.
        type(ten_2D4O2sym) :: res
            !! Difference.
        res%vals = a%vals - b%vals
    end function sub_2D4O2sym

    pure function subU_2D4O2sym(a) result(res)
        !! Unary negation.
        implicit none
        type(ten_2D4O2sym), intent(in) :: a
            !! Operand.
        type(ten_2D4O2sym) :: res
            !! Negated tensor.
        res%vals = -a%vals
    end function subU_2D4O2sym

    pure function mul_real64_2D4O2sym(a, b) result(res)
        !! Scalar product \(a\,\mathbb{B}\).
        implicit none
        real(real64), intent(in) :: a
            !! Scalar factor.
        type(ten_2D4O2sym), intent(in) :: b
            !! Tensor.
        type(ten_2D4O2sym) :: res
            !! Scaled tensor.
        res%vals = a * b%vals
    end function mul_real64_2D4O2sym

    pure function mul_2D4O2sym_real64(b, a) result(res)
        !! Scalar product \(\mathbb{B}\,a\).
        implicit none
        type(ten_2D4O2sym), intent(in) :: b
            !! Tensor.
        real(real64), intent(in) :: a
            !! Scalar factor.
        type(ten_2D4O2sym) :: res
            !! Scaled tensor.
        res%vals = a * b%vals
    end function mul_2D4O2sym_real64

    pure function div_2D4O2sym_real64(b, a) result(res)
        !! Scalar division \(\mathbb{B}/a\).
        implicit none
        type(ten_2D4O2sym), intent(in) :: b
            !! Tensor.
        real(real64), intent(in) :: a
            !! Scalar divisor (must be non-zero).
        type(ten_2D4O2sym) :: res
            !! Scaled tensor.
        res%vals = b%vals/a
    end function div_2D4O2sym_real64

    pure function inv_2D4O2sym(a) result(res)
        !! `.inv.` Computes the inverse of a ten_2D4O2sym tensor, i.e. \(\mathbb{A}^{-1}\) such that
        !! \(\mathbb{A}^{-1} : \mathbb{A} = \mathbb{A} : \mathbb{A}^{-1} = \mathbb{I}^S\).
        !!
        !! The 4x4 Voigt matrix is inverted with `M44INV` (LAPACK LU) and the shear
        !! row/column are scaled by 1/2 (shear-shear term by 1/4) to stay consistent
        !! with the tensorial-shear convention of `.ddot.` (same scheme as `inv_3D4O2sym`).
        !! A zero tensor is returned if the matrix is singular.
        use muscle_math_inverses, only : M44INV
        implicit none
        type(ten_2D4O2sym), intent(in) :: a
            !! Tensor to invert.
        type(ten_2D4O2sym) :: res
            !! Inverse tensor.
        real(real64) :: mat_b(4,4)
        logical :: ok

        call M44INV(a%vals, mat_b, ok)

        mat_b(1:3, 4) = mat_b(1:3, 4)/2D0
        mat_b(4, 1:3) = mat_b(4, 1:3)/2D0
        mat_b(4, 4)   = mat_b(4, 4)/4D0

        call res%init(mat_b)
    end function inv_2D4O2sym

    pure function convert_3sym(self) result(res)
        !! Converts to a fully symmetric `ten_2D4O3sym` by keeping the upper triangle
        !! of the Voigt matrix (no symmetrization is performed; meaningful only when
        !! the tensor already has major symmetry).
        use muscle_tensor_2d4o3sym, only : ten_2D4O3sym
        implicit none
        class(ten_2D4O2sym), intent(in) :: self
            !! Tensor to convert.
        type(ten_2D4O3sym) :: res
            !! Fully symmetric tensor.

        call res%init((/ self%vals(1,1), self%vals(2,2), self%vals(3,3), self%vals(4,4), &
                         self%vals(1,2), self%vals(2,3), self%vals(3,4),                 &
                         self%vals(1,3), self%vals(2,4),                                 &
                         self%vals(1,4)                                                  &
                      /))
    end function convert_3sym

end module muscle_tensor_2d4o2sym
