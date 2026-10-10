! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_2d2o
    !! Module muscle_tensor_2d2o
    !! =========================
    !!
    !! Defines the type for general (non-symmetric) 2D second-order tensors and associated operations.
    !!
    !! This module provides the derived type `ten_2D2O` to represent a general
    !! second-order tensor in a 2D setting (plane or axisymmetric), for example the
    !! deformation gradient of a plane-strain or axisymmetric problem:
    !! \[ \mathbf{A} = \begin{bmatrix} A_{11} & A_{12} & 0 \\ A_{21} & A_{22} & 0 \\ 0 & 0 & A_{33} \end{bmatrix} \]
    !! Unlike `ten_2D2Osym`, no symmetry is assumed, so the 5 non-zero components are
    !! stored explicitly: the in-plane 2x2 block in column-major order followed by
    !! the out-of-plane component, (11, 21, 12, 22, 33).
    !!
    !! The module overloads standard arithmetic operators (+, -, *, /), a custom
    !! equality comparison operator (.approx.), the deviatoric operator (.dev.),
    !! the double dot product operator (.ddot.) and the intrinsic `transpose`
    !! for this tensor type. It also provides methods for initialization and
    !! assignment from a scalar.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ten_2D2O`: Represents a general 2D second-order tensor.
    !!     - Component: `vals(5) :: real(real64)` - Stores the 5 components in
    !!       the order (11, 21, 12, 22, 33).
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       5-element array or from individual xx, xy, yx, yy, zz components.
    !!     - Procedures: `norm`, `is_approx`, `det`, `transpose`.
    !!
    !! ### Operators:
    !!
    !! - `.approx.`: Compares two `ten_2D2O` tensors for approximate equality.
    !! - `+`: Adds two `ten_2D2O` tensors.
    !! - `-`: Subtracts two `ten_2D2O` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_2D2O` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_2D2O` tensor by a `real(real64)` scalar.
    !! - `.dev.`: Computes the deviatoric part of a `ten_2D2O` tensor.
    !! - `.ddot.`: Computes the double dot product (Frobenius inner product, scalar result) of two `ten_2D2O` tensors.
    !! - `transpose`: Returns the transposed tensor.
    !!
    !! ### Assignment:
    !!
    !! - `=`: Allows assigning a single `real(real64)` scalar value to all components
    !!        of a `ten_2D2O` tensor.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_2d2o_usage
    !!   use muscle_tensors
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_2D2O) :: F, Ft
    !!   real(real64) :: J
    !!
    !!   ! Simple shear plus out-of-plane stretch (xx, xy, yx, yy, zz)
    !!   call F%init(xx=1.0D0, xy=0.3D0, yx=0.0D0, yy=1.0D0, zz=1.1D0)
    !!
    !!   J  = F%det()          ! = 1.1
    !!   Ft = transpose(F)     ! Ft%vals = (1.0, 0.3, 0.0, 1.0, 1.1)
    !!
    !! end program example_ten_2d2o_usage
    !! ```
    !!
    !! For more information see [[muscle_tensors]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: ten_2D2O
        !! General 2D Second-Order Tensor (5-Component Storage)
        !! =====================================================
        !!
        !! Represents a general (potentially non-symmetric) second-order tensor in a
        !! 2D setting, \(A_{ij}\), whose out-of-plane coupling terms vanish
        !! (\(A_{13} = A_{23} = A_{31} = A_{32} = 0\)).
        !!
        !! Storage:
        !! --------
        !! The tensor components are stored internally in a 1D array `vals` of size 5:
        !! the in-plane 2x2 block in **column-major** ordering (consistent with
        !! `ten_3D2O`), followed by the out-of-plane component:
        !!
        !! - `vals(1)`: Component (1,1) or xx
        !! - `vals(2)`: Component (2,1) or yx
        !! - `vals(3)`: Component (1,2) or xy
        !! - `vals(4)`: Component (2,2) or yy
        !! - `vals(5)`: Component (3,3) or zz
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 5-element `real(real64)`
        !! array (storage order above) or by providing the 5 components individually
        !! in the order: (xx, xy, yx, yy, zz). Note that the `init2` procedure
        !! internally rearranges these into the storage order.
        !!
        !! For more information see [[muscle_tensors]]

        real(real64), dimension(5) :: vals
            !! Stores the 5 components in the order: (11, 21, 12, 22, 33).
    contains
        generic, public :: init => init_ten_2D2O, init2_ten_2D2O
            !! Generic interface for initialization.
        procedure, private :: init_ten_2D2O, init2_ten_2D2O
        procedure, public :: norm => norm_2D2O
            !! Computes the L1 norm of the tensor.
        procedure, public :: is_approx => is_approx_2D2O
            !! Compares two tensors for approximate equality.
        procedure, public :: det => det_2D2O
            !! Computes the determinant of the tensor.
        procedure, public :: transpose => transpose_2D2O
            !! Returns the transposed tensor.
    end type ten_2D2O

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_2D2O
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_2D2O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D2O
        module procedure subU_2D2O
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_2D2O
        module procedure mul_2D2O_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_2D2O_real64
    end interface

    public :: operator(.dev.)
    interface operator (.dev.)
        module procedure dev_2D2O
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_2D2O_2D2O
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure ten_2D2O_real64_assign
    end interface

#ifdef ENABLE_UDTIO
    public :: write(formatted)
    interface write(formatted)
        module procedure print_ten_2D2O
    end interface
#endif

    public :: transpose
    interface transpose
        module procedure transpose_2D2O
    end interface
contains

    pure subroutine ten_2D2O_real64_assign(a, b)
        !! Assigns the scalar `b` to every stored component of `a`.
        implicit none
        type(ten_2D2O), intent(out) :: a
            !! Target tensor.
        real(real64), intent(in) :: b
            !! Scalar value.
        a%vals = b
    end subroutine ten_2D2O_real64_assign

    pure subroutine init_ten_2D2O(self, vals)
        !! Initializes a ten_2D2O tensor from a 5-element array (storage order).
        implicit none
        class(ten_2D2O), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in) :: vals(5)
            !! Components in the order (11, 21, 12, 22, 33).
        self%vals = vals
    end subroutine init_ten_2D2O

    pure subroutine init2_ten_2D2O(self, xx, xy, yx, yy, zz)
        !! Initializes a ten_2D2O tensor from its 5 individual components.
        !! Input order is (xx, xy, yx, yy, zz) (row-wise, as in `ten_3D2O`).
        !! Internal storage is (xx, yx, xy, yy, zz).
        implicit none
        class(ten_2D2O), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in) :: xx
            !! Component (1,1).
        real(real64), intent(in) :: xy
            !! Component (1,2).
        real(real64), intent(in) :: yx
            !! Component (2,1).
        real(real64), intent(in) :: yy
            !! Component (2,2).
        real(real64), intent(in) :: zz
            !! Component (3,3).
        self%vals = (/xx, yx, xy, yy, zz/)
    end subroutine init2_ten_2D2O

    pure function norm_2D2O(a) result(res)
        !! L1 norm of the 5 stored components: \(\sum |A_{ij}|\).
        implicit none
        class(ten_2D2O), intent(in) :: a
            !! Tensor.
        real(real64) :: res
            !! Norm.
        res = sum(abs(a%vals))
    end function norm_2D2O

    pure function is_approx_2D2O(a, b, tol) result(res)
        !! Component-wise approximate equality with a relative tolerance:
        !! \( |a_I - b_I| \le tol \cdot \max(\max|a|, \max|b|, 10^{-30}) \) for all I.
        implicit none
        class(ten_2D2O), intent(in) :: a
            !! First tensor.
        class(ten_2D2O), intent(in) :: b
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
    end function is_approx_2D2O

    pure function approx_2D2O(a, b) result(res)
        !! `.approx.` Compares two ten_2D2O tensors for approximate equality.
        !! Delegates to `is_approx` with its default tolerance.
        implicit none
        type(ten_2D2O), intent(in) :: a
            !! First tensor.
        type(ten_2D2O), intent(in) :: b
            !! Second tensor.
        logical :: res
            !! `.TRUE.` if both tensors are approximately equal.

        res = a%is_approx(b)
    end function approx_2D2O

    pure function sum_2D2O(a, b) result(res)
        !! Component-wise sum.
        implicit none
        type(ten_2D2O), intent(in) :: a
            !! First operand.
        type(ten_2D2O), intent(in) :: b
            !! Second operand.
        type(ten_2D2O) :: res
            !! Sum.
        res%vals = a%vals + b%vals
    end function sum_2D2O

    pure function sub_2D2O(a, b) result(res)
        !! Component-wise difference.
        implicit none
        type(ten_2D2O), intent(in) :: a
            !! Minuend.
        type(ten_2D2O), intent(in) :: b
            !! Subtrahend.
        type(ten_2D2O) :: res
            !! Difference.
        res%vals = a%vals - b%vals
    end function sub_2D2O

    pure function subU_2D2O(a) result(res)
        !! Unary negation.
        implicit none
        type(ten_2D2O), intent(in) :: a
            !! Operand.
        type(ten_2D2O) :: res
            !! Negated tensor.
        res%vals = -a%vals
    end function subU_2D2O

    pure function mul_real64_2D2O(a, b) result(res)
        !! Scalar product \(a\,\mathbf{B}\).
        implicit none
        real(real64), intent(in) :: a
            !! Scalar factor.
        type(ten_2D2O), intent(in) :: b
            !! Tensor.
        type(ten_2D2O) :: res
            !! Scaled tensor.
        res%vals = a * b%vals
    end function mul_real64_2D2O

    pure function mul_2D2O_real64(a, b) result(res)
        !! Scalar product \(\mathbf{A}\,b\).
        implicit none
        type(ten_2D2O), intent(in) :: a
            !! Tensor.
        real(real64), intent(in) :: b
            !! Scalar factor.
        type(ten_2D2O) :: res
            !! Scaled tensor.
        res%vals = a%vals * b
    end function mul_2D2O_real64

    pure function div_2D2O_real64(a, b) result(res)
        !! Scalar division \(\mathbf{A}/b\).
        implicit none
        type(ten_2D2O), intent(in) :: a
            !! Tensor.
        real(real64), intent(in) :: b
            !! Scalar divisor (must be non-zero).
        type(ten_2D2O) :: res
            !! Scaled tensor.
        res%vals = a%vals/b
    end function div_2D2O_real64

    pure function ddot_2D2O_2D2O(a, b) result(res)
        !! Double contraction (Frobenius inner product) \(\mathbf{A}:\mathbf{B} = A_{ij}B_{ij}\).
        implicit none
        type(ten_2D2O), intent(in) :: a
            !! First operand.
        type(ten_2D2O), intent(in) :: b
            !! Second operand.
        real(real64) :: res
            !! Scalar result.
        res = sum(a%vals * b%vals)
    end function ddot_2D2O_2D2O

    pure function dev_2D2O(a) result(res)
        !! Deviatoric part \(\mathbf{A} - \frac{1}{3}\mathrm{tr}(\mathbf{A})\mathbf{I}\).
        implicit none
        type(ten_2D2O), intent(in) :: a
            !! Operand.
        type(ten_2D2O) :: res
            !! Deviatoric tensor.
        real(real64) :: hydro
        hydro = (a%vals(1) + a%vals(4) + a%vals(5))/3D0
        res%vals(1) = a%vals(1) - hydro
        res%vals(2) = a%vals(2)
        res%vals(3) = a%vals(3)
        res%vals(4) = a%vals(4) - hydro
        res%vals(5) = a%vals(5) - hydro
    end function dev_2D2O

    pure function det_2D2O(a) result(res)
        !! Determinant \(\det\mathbf{A} = (A_{11}A_{22} - A_{12}A_{21})\,A_{33}\).
        implicit none
        class(ten_2D2O), intent(in) :: a
            !! Tensor.
        real(real64) :: res
            !! Determinant.
        res = (a%vals(1)*a%vals(4) - a%vals(2)*a%vals(3)) * a%vals(5)
    end function det_2D2O

#ifdef ENABLE_UDTIO
    subroutine print_ten_2D2O(dtv, unit, iotype, v_list, iostat, iomsg)
        !! Custom I/O formatting for ten_2D2O.
        !! List mode prints the 5 stored components; matrix mode prints the full 3x3 matrix.
        class(ten_2D2O), intent(in)     :: dtv
            !! Tensor to print.
        integer, intent(in)             :: unit
            !! Output unit.
        character(len=*), intent(in)    :: iotype
            !! Edit descriptor type.
        integer, intent(in)             :: v_list(:)
            !! Optional (width, decimals).
        integer, intent(out)            :: iostat
            !! I/O status.
        character(len=*), intent(inout) :: iomsg
            !! I/O error message.

        character(len=100) :: fmt_string
        integer :: w, d

        w = 11
        d = 4
        iostat = 0

        if (size(v_list) >= 1) w = v_list(1)
        if (size(v_list) >= 2) d = v_list(2)

        if (iotype == "DTLIST" .or. iotype == "LIST") then
            write(fmt_string, "('(A, 5(ES', I0, '.', I0, ', 1X), A)')") w, d
            write(unit, fmt_string, iostat=iostat) &
                "[", dtv%vals(1), dtv%vals(2), dtv%vals(3), dtv%vals(4), dtv%vals(5), "]"
        else
            write(fmt_string, "('(/, 3(ES', I0, '.', I0, ', 2X), /, 3(ES', I0, '.', I0, ', 2X), /, 3(ES', I0, '.', I0, ', 2X))')") &
                w, d, w, d, w, d

            write(unit, fmt_string, iostat=iostat) &
                dtv%vals(1), dtv%vals(3), 0D0,         & ! Row 1
                dtv%vals(2), dtv%vals(4), 0D0,         & ! Row 2
                0D0,         0D0,         dtv%vals(5)    ! Row 3
        end if

        if (iostat /= 0) then
            iomsg = "Error in print_ten_2D2O: Failed to write to the specified unit."
        end if
    end subroutine print_ten_2D2O
#endif

    pure function transpose_2D2O(self) result(res)
        !! Returns the transpose of a general 2D second-order tensor (\(A^T_{ij} = A_{ji}\)).
        implicit none
        class(ten_2D2O), intent(in) :: self
            !! Tensor.
        type(ten_2D2O) :: res
            !! Transposed tensor.

        res%vals(1) = self%vals(1) ! xx -> xx
        res%vals(2) = self%vals(3) ! yx <- xy
        res%vals(3) = self%vals(2) ! xy <- yx
        res%vals(4) = self%vals(4) ! yy -> yy
        res%vals(5) = self%vals(5) ! zz -> zz
    end function transpose_2D2O

end module muscle_tensor_2d2o
