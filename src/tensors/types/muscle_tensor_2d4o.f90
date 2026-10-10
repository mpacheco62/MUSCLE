! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_2d4o
    !! Module muscle_tensor_2d4o
    !! =========================
    !! Defines a general (unsymmetric) 2D fourth-order tensor, the 2D counterpart of
    !! `ten_3D4O`. It acts on 2D general second-order tensors (`ten_2D2O`, 5 components:
    !! 11, 21, 12, 22, 33), so it is stored as a 5x5 matrix (25 components).
    !!
    !! Useful for tangents lacking minor/major symmetries in plane or axisymmetric
    !! problems, such as the 1st Piola-Kirchhoff tangent stiffness
    !! \( \mathbb{A}_{iAjB} = \partial P_{iA} / \partial F_{jB} \).
    !!
    !! Index mapping
    !! -------------
    !! The 5D index \(I(i,j)\) of a pair \((i,j)\) follows the `ten_2D2O` storage:
    !! (1,1)->1, (2,1)->2, (1,2)->3, (2,2)->4, (3,3)->5. Pairs that couple the plane with
    !! the out-of-plane direction ((1,3), (3,1), (2,3), (3,2)) are identically zero.
    !! \(C_{ijkl}\) is stored in `vals(I(i,j), I(k,l))`.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! - `ten_2D4O`: type with `vals(5,5)`, `init` (5x5 matrix or `C4(3,3,3,3)`),
    !!   `get`, `set`, `is_approx`, `norm`, `convert_2sym`.
    !! - Operators: `.approx.`, `+`, `-`, `*`, `/`, `.inv.`, assignment from `real(real64)`.
    !!
    !! For more information see [[muscle_tensors]]

    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensor_2d4o2sym, only : ten_2D4O2sym
    use muscle_math_inverses, only : M55INV
    implicit none
    private

    type, public :: ten_2D4O
        !! Unsymmetric 2D Fourth-Order Tensor (25 components)
        !! ==================================================
        !! Stored as a 5x5 matrix whose row/column indices follow the `ten_2D2O`
        !! storage: \(C_{ijkl}\) lives in `vals(I(i,j), I(k,l))`, with
        !! I(1,1)=1, I(2,1)=2, I(1,2)=3, I(2,2)=4, I(3,3)=5.
        real(real64), dimension(5,5) :: vals = 0.0D0
            !! 5x5 matrix representation.
    contains
        generic, public :: init => init_matrix_5x5, init_array_3333
            !! Generic interface for initialization.
        procedure, private :: init_matrix_5x5
        procedure, private :: init_array_3333

        procedure, public :: get => get_component
            !! Returns \(C_{ijkl}\) (zero for out-of-plane coupling indices).
        procedure, public :: set => set_component
            !! Sets \(C_{ijkl}\) (ignored for out-of-plane coupling indices).
        procedure, public :: is_approx => is_approx_2D4O
            !! Compares two tensors for approximate equality.
        procedure, public :: norm => norm_2D4O
            !! L1 norm of the 5x5 matrix.
        procedure, public :: convert_2sym => convert_to_2D4O2sym
            !! Projects onto a minor-symmetric 4x4 Voigt `ten_2D4O2sym`.
    end type ten_2D4O

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_2D4O
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_2D4O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D4O
        module procedure subU_2D4O
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_2D4O
        module procedure mul_2D4O_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_2D4O_real64
    end interface

    public :: operator(.inv.)
    interface operator (.inv.)
        module procedure inv_2D4O
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_2D4O_real64
    end interface

#ifdef ENABLE_UDTIO
    public :: write(formatted)
    interface write(formatted)
        module procedure print_ten_2D4O
    end interface
#endif

contains

    pure function idx5(i, j) result(res)
        !! 5D index of the pair (i,j) in the `ten_2D2O` storage, or 0 for the
        !! (identically zero) out-of-plane coupling pairs.
        implicit none
        integer, intent(in) :: i
            !! First index (1..3).
        integer, intent(in) :: j
            !! Second index (1..3).
        integer :: res
            !! Storage index (1..5) or 0.
        integer, parameter :: map(3,3) = reshape((/ 1, 2, 0, &
                                                    3, 4, 0, &
                                                    0, 0, 5 /), (/3,3/))
        res = map(i, j)
    end function idx5

    pure subroutine init_matrix_5x5(self, vals)
        !! Initializes from a 5x5 matrix.
        implicit none
        class(ten_2D4O), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in)       :: vals(5,5)
            !! 5x5 matrix representation.
        self%vals = vals
    end subroutine init_matrix_5x5

    pure subroutine init_array_3333(self, C4)
        !! Initializes from a full 3x3x3x3 array, keeping only the 2D components
        !! (out-of-plane coupling components of `C4` are discarded).
        !!
        !! The 25 components are copied explicitly, column by column, using the
        !! mapping I(1,1)=1, I(2,1)=2, I(1,2)=3, I(2,2)=4, I(3,3)=5. This avoids the
        !! 81-iteration loop with index look-ups and branches (about 20x faster at -O0,
        !! same speed at -O3).
        implicit none
        class(ten_2D4O), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in)       :: C4(3,3,3,3)
            !! Full fourth-order array \(C_{ijkl}\).

        ! Column 1: (k,l) = (1,1)
        self%vals(1,1) = C4(1,1,1,1); self%vals(2,1) = C4(2,1,1,1); self%vals(3,1) = C4(1,2,1,1)
        self%vals(4,1) = C4(2,2,1,1); self%vals(5,1) = C4(3,3,1,1)
        ! Column 2: (k,l) = (2,1)
        self%vals(1,2) = C4(1,1,2,1); self%vals(2,2) = C4(2,1,2,1); self%vals(3,2) = C4(1,2,2,1)
        self%vals(4,2) = C4(2,2,2,1); self%vals(5,2) = C4(3,3,2,1)
        ! Column 3: (k,l) = (1,2)
        self%vals(1,3) = C4(1,1,1,2); self%vals(2,3) = C4(2,1,1,2); self%vals(3,3) = C4(1,2,1,2)
        self%vals(4,3) = C4(2,2,1,2); self%vals(5,3) = C4(3,3,1,2)
        ! Column 4: (k,l) = (2,2)
        self%vals(1,4) = C4(1,1,2,2); self%vals(2,4) = C4(2,1,2,2); self%vals(3,4) = C4(1,2,2,2)
        self%vals(4,4) = C4(2,2,2,2); self%vals(5,4) = C4(3,3,2,2)
        ! Column 5: (k,l) = (3,3)
        self%vals(1,5) = C4(1,1,3,3); self%vals(2,5) = C4(2,1,3,3); self%vals(3,5) = C4(1,2,3,3)
        self%vals(4,5) = C4(2,2,3,3); self%vals(5,5) = C4(3,3,3,3)
    end subroutine init_array_3333

    pure function get_component(self, i, j, k, l) result(res)
        !! Returns \(C_{ijkl}\); zero for out-of-plane coupling index pairs.
        implicit none
        class(ten_2D4O), intent(in) :: self
            !! Tensor.
        integer, intent(in)         :: i, j, k, l
            !! Tensor indices (1..3).
        real(real64)                :: res
            !! Component value.
        integer :: row, col
        row = idx5(i, j)
        col = idx5(k, l)
        res = 0.0D0
        if (row > 0 .and. col > 0) res = self%vals(row, col)
    end function get_component

    pure subroutine set_component(self, i, j, k, l, val)
        !! Sets \(C_{ijkl}\); silently ignored for out-of-plane coupling index pairs.
        implicit none
        class(ten_2D4O), intent(inout) :: self
            !! Tensor.
        integer, intent(in)            :: i, j, k, l
            !! Tensor indices (1..3).
        real(real64), intent(in)       :: val
            !! Value to set.
        integer :: row, col
        row = idx5(i, j)
        col = idx5(k, l)
        if (row > 0 .and. col > 0) self%vals(row, col) = val
    end subroutine set_component

    pure function norm_2D4O(self) result(res)
        !! L1 norm \(\sum |A_{IJ}|\) of the 5x5 matrix.
        implicit none
        class(ten_2D4O), intent(in) :: self
            !! Tensor.
        real(real64)                :: res
            !! Norm.
        res = sum(abs(self%vals))
    end function norm_2D4O

    pure function is_approx_2D4O(self, other, tol) result(res)
        !! Component-wise approximate equality with a relative tolerance:
        !! \( |a_{IJ} - b_{IJ}| \le tol \cdot \max(\max|a|, \max|b|, 10^{-30}) \).
        implicit none
        class(ten_2D4O), intent(in)         :: self
            !! First tensor.
        class(ten_2D4O), intent(in)         :: other
            !! Second tensor.
        real(real64), optional, intent(in)  :: tol
            !! Relative tolerance (default `1E-12`).
        logical                             :: res
            !! `.TRUE.` if both tensors are approximately equal.
        real(real64), parameter             :: EPS_ABS = 1.0D-30
        real(real64)                        :: tol2, max_val, eps_check

        tol2 = 1.0D-12
        if (present(tol)) tol2 = tol
        max_val = max(maxval(abs(self%vals)), maxval(abs(other%vals)), EPS_ABS)
        eps_check = tol2 * max_val

        res = all(abs(self%vals - other%vals) <= eps_check)
    end function is_approx_2D4O

    pure function approx_2D4O(a, b) result(res)
        !! `.approx.` Delegates to `is_approx` with its default tolerance.
        implicit none
        type(ten_2D4O), intent(in) :: a
            !! First tensor.
        type(ten_2D4O), intent(in) :: b
            !! Second tensor.
        logical                    :: res
            !! `.TRUE.` if both tensors are approximately equal.
        res = a%is_approx(b)
    end function approx_2D4O

    pure function sum_2D4O(a, b) result(res)
        !! Component-wise sum.
        implicit none
        type(ten_2D4O), intent(in) :: a
            !! First operand.
        type(ten_2D4O), intent(in) :: b
            !! Second operand.
        type(ten_2D4O)             :: res
            !! Sum.
        res%vals = a%vals + b%vals
    end function sum_2D4O

    pure function sub_2D4O(a, b) result(res)
        !! Component-wise difference.
        implicit none
        type(ten_2D4O), intent(in) :: a
            !! Minuend.
        type(ten_2D4O), intent(in) :: b
            !! Subtrahend.
        type(ten_2D4O)             :: res
            !! Difference.
        res%vals = a%vals - b%vals
    end function sub_2D4O

    pure function subU_2D4O(a) result(res)
        !! Unary negation.
        implicit none
        type(ten_2D4O), intent(in) :: a
            !! Operand.
        type(ten_2D4O)             :: res
            !! Negated tensor.
        res%vals = -a%vals
    end function subU_2D4O

    pure function mul_real64_2D4O(a, b) result(res)
        !! Scalar product \(a\,\mathbb{B}\).
        implicit none
        real(real64), intent(in)   :: a
            !! Scalar factor.
        type(ten_2D4O), intent(in) :: b
            !! Tensor.
        type(ten_2D4O)             :: res
            !! Scaled tensor.
        res%vals = a * b%vals
    end function mul_real64_2D4O

    pure function mul_2D4O_real64(b, a) result(res)
        !! Scalar product \(\mathbb{B}\,a\).
        implicit none
        type(ten_2D4O), intent(in) :: b
            !! Tensor.
        real(real64), intent(in)   :: a
            !! Scalar factor.
        type(ten_2D4O)             :: res
            !! Scaled tensor.
        res%vals = b%vals * a
    end function mul_2D4O_real64

    pure function div_2D4O_real64(b, a) result(res)
        !! Scalar division \(\mathbb{B}/a\).
        implicit none
        type(ten_2D4O), intent(in) :: b
            !! Tensor.
        real(real64), intent(in)   :: a
            !! Scalar divisor (must be non-zero).
        type(ten_2D4O)             :: res
            !! Scaled tensor.
        res%vals = b%vals / a
    end function div_2D4O_real64

    pure subroutine assign_2D4O_real64(a, b)
        !! Assigns the scalar `b` to every component of `a`.
        implicit none
        type(ten_2D4O), intent(out) :: a
            !! Target tensor.
        real(real64), intent(in)    :: b
            !! Scalar value.
        a%vals = b
    end subroutine assign_2D4O_real64

    pure function convert_to_2D4O2sym(self) result(res)
        !! Projects onto the 4x4 Voigt `ten_2D4O2sym` by sampling the (xx, yy, zz, xy)
        !! rows/columns (no symmetrization; meaningful when minor symmetries hold).
        implicit none
        class(ten_2D4O), intent(in) :: self
            !! Tensor to convert.
        type(ten_2D4O2sym)          :: res
            !! Minor-symmetric tensor.
        integer, parameter          :: map4(4) = (/1, 4, 5, 3/) ! xx,yy,zz,xy in 5D

        integer :: I, J
        do J = 1, 4
            do I = 1, 4
                res%vals(I, J) = self%vals(map4(I), map4(J))
            end do
        end do
    end function convert_to_2D4O2sym

    pure function inv_2D4O(a) result(res)
        !! `.inv.` Inverse of the 5x5 matrix representation (delegates to `M55INV`).
        !! A zero tensor is returned if the matrix is singular.
        implicit none
        type(ten_2D4O), intent(in) :: a
            !! Tensor to invert.
        type(ten_2D4O)             :: res
            !! Inverse tensor.
        logical                    :: ok

        call M55INV(a%vals, res%vals, ok)
        if (.not. ok) res%vals = 0.0D0
    end function inv_2D4O

#ifdef ENABLE_UDTIO
    subroutine print_ten_2D4O(dtv, unit, iotype, v_list, iostat, iomsg)
        !! Custom I/O formatting for ten_2D4O (prints the 5x5 matrix row by row).
        class(ten_2D4O), intent(in)     :: dtv
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
        integer :: w, d, r

        w = 11; d = 4; iostat = 0
        if (size(v_list) >= 1) w = v_list(1)
        if (size(v_list) >= 2) d = v_list(2)

        write(fmt_string, "('(/, 5(ES', I0, '.', I0, ', 1X))')") w, d
        do r = 1, 5
            write(unit, fmt_string, iostat=iostat) dtv%vals(r, :)
        end do
        if (iostat /= 0) then
            iomsg = "Error in print_ten_2D4O: Failed to write to the specified unit."
        end if
    end subroutine print_ten_2D4O
#endif

end module muscle_tensor_2d4o
