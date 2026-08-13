! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_3d4o
    !! Module muscle_tensor_3d4o
    !! ==========================
    !! Defines a general (unsymmetric) 3D fourth-order tensor (81 components)
    !! stored as a 9x9 matrix in column-major 9D vector space.
    !!
    !! Useful for tangents lacking minor/major symmetries, such as the 1st Piola-Kirchhoff 
    !! tangent stiffness A_iAjB = dP_iA / dF_jB.

    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensor_3d2o, only : ten_3D2O
    use muscle_tensor_3d4o2sym, only : ten_3D4O2sym
    use muscle_math_inverses, only : M99INV
    implicit none
    private

    type, public :: ten_3D4O
        !! Unsymmetric 3D Fourth-Order Tensor (81 components)
        !! Stored as a 9x9 matrix corresponding to column-major 3x3 tensor indices:
        !! Row I = 3*(j-1) + i,   Col J = 3*(l-1) + k
        real(real64), dimension(9,9) :: vals = 0.0D0
    contains
        generic, public :: init => init_matrix_9x9, init_array_3333
        procedure, private :: init_matrix_9x9
        procedure, private :: init_array_3333

        procedure, public :: get => get_component
        procedure, public :: set => set_component
        procedure, public :: is_approx => is_approx_3D4O
        procedure, public :: norm => norm_3D4O
        procedure, public :: convert_2sym => convert_to_3D4O2sym
    end type ten_3D4O

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_3D4O
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O
        module procedure subU_3D4O
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_3D4O
        module procedure mul_3D4O_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_3D4O_real64
    end interface

    public :: operator(.inv.)
    interface operator (.inv.)
        module procedure inv_3D4O
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_3D4O_real64
    end interface

#ifdef ENABLE_UDTIO
    public :: write(formatted)
    interface write(formatted)
        module procedure print_ten_3D4O
    end interface
#endif

contains

    pure subroutine init_matrix_9x9(self, vals)
        implicit none
        class(ten_3D4O), intent(inout) :: self
        real(real64), intent(in)       :: vals(9,9)
        self%vals = vals
    end subroutine init_matrix_9x9

    pure subroutine init_array_3333(self, C4)
        implicit none
        class(ten_3D4O), intent(inout) :: self
        real(real64), intent(in)       :: C4(3,3,3,3)
        integer :: i, j, k, l, row, col

        do l = 1, 3
            do k = 1, 3
                col = 3*(l-1) + k
                do j = 1, 3
                    do i = 1, 3
                        row = 3*(j-1) + i
                        self%vals(row, col) = C4(i,j,k,l)
                    end do
                end do
            end do
        end do
    end subroutine init_array_3333

    pure function get_component(self, i, j, k, l) result(res)
        implicit none
        class(ten_3D4O), intent(in) :: self
        integer, intent(in)         :: i, j, k, l
        real(real64)                :: res
        res = self%vals(3*(j-1)+i, 3*(l-1)+k)
    end function get_component

    pure subroutine set_component(self, i, j, k, l, val)
        implicit none
        class(ten_3D4O), intent(inout) :: self
        integer, intent(in)            :: i, j, k, l
        real(real64), intent(in)       :: val
        self%vals(3*(j-1)+i, 3*(l-1)+k) = val
    end subroutine set_component

    pure function norm_3D4O(self) result(res)
        implicit none
        class(ten_3D4O), intent(in) :: self
        real(real64)                :: res
        res = sum(abs(self%vals))
    end function norm_3D4O

    pure function is_approx_3D4O(self, other, tol) result(res)
        implicit none
        class(ten_3D4O), intent(in)         :: self, other
        real(real64), optional, intent(in)  :: tol
        logical                             :: res
        real(real64), parameter             :: EPS_ABS = 1.0D-30
        real(real64)                        :: tol2, max_val, eps_check

        tol2 = 1.0D-12
        if (present(tol)) tol2 = tol
        max_val = max(maxval(abs(self%vals)), maxval(abs(other%vals)), EPS_ABS)
        eps_check = tol2 * max_val

        res = all(abs(self%vals - other%vals) <= eps_check)
    end function is_approx_3D4O

    pure function approx_3D4O(a, b) result(res)
        implicit none
        type(ten_3D4O), intent(in) :: a, b
        logical                    :: res
        res = a%is_approx(b)
    end function approx_3D4O

    pure function sum_3D4O(a, b) result(res)
        implicit none
        type(ten_3D4O), intent(in) :: a, b
        type(ten_3D4O)             :: res
        res%vals = a%vals + b%vals
    end function sum_3D4O

    pure function sub_3D4O(a, b) result(res)
        implicit none
        type(ten_3D4O), intent(in) :: a, b
        type(ten_3D4O)             :: res
        res%vals = a%vals - b%vals
    end function sub_3D4O

    pure function subU_3D4O(a) result(res)
        implicit none
        type(ten_3D4O), intent(in) :: a
        type(ten_3D4O)             :: res
        res%vals = -a%vals
    end function subU_3D4O

    pure function mul_real64_3D4O(a, b) result(res)
        implicit none
        real(real64), intent(in)   :: a
        type(ten_3D4O), intent(in) :: b
        type(ten_3D4O)             :: res
        res%vals = a * b%vals
    end function mul_real64_3D4O

    pure function mul_3D4O_real64(b, a) result(res)
        implicit none
        type(ten_3D4O), intent(in) :: b
        real(real64), intent(in)   :: a
        type(ten_3D4O)             :: res
        res%vals = b%vals * a
    end function mul_3D4O_real64

    pure function div_3D4O_real64(b, a) result(res)
        implicit none
        type(ten_3D4O), intent(in) :: b
        real(real64), intent(in)   :: a
        type(ten_3D4O)             :: res
        res%vals = b%vals / a
    end function div_3D4O_real64

    pure subroutine assign_3D4O_real64(a, b)
        implicit none
        type(ten_3D4O), intent(out) :: a
        real(real64), intent(in)    :: b
        a%vals = b
    end subroutine assign_3D4O_real64

    pure function convert_to_3D4O2sym(self) result(res)
        !! Symmetrizes minor symmetries and projects to 6x6 Voigt ten_3D4O2sym
        implicit none
        class(ten_3D4O), intent(in) :: self
        type(ten_3D4O2sym)          :: res
        integer, parameter          :: map6(6) = (/1, 5, 9, 4, 8, 7/) ! xx,yy,zz,xy,yz,xz in 9D

        integer :: I, J
        do I = 1, 6
            do J = 1, 6
                res%vals(I, J) = self%vals(map6(I), map6(J))
            end do
        end do
    end function convert_to_3D4O2sym

    pure function inv_3D4O(a) result(res)
        !! Computes the inverse of ten_3D4O by delegating to M99INV in muscle_math_inverses.
        implicit none
        type(ten_3D4O), intent(in) :: a
        type(ten_3D4O)             :: res
        logical                    :: ok

        call M99INV(a%vals, res%vals, ok)
        if (.not. ok) res%vals = 0.0D0
    end function inv_3D4O

#ifdef ENABLE_UDTIO
    subroutine print_ten_3D4O(dtv, unit, iotype, v_list, iostat, iomsg)
        class(ten_3D4O), intent(in)     :: dtv
        integer, intent(in)             :: unit
        character(len=*), intent(in)    :: iotype
        integer, intent(in)             :: v_list(:)
        integer, intent(out)            :: iostat
        character(len=*), intent(inout) :: iomsg
        character(len=100) :: fmt_string
        integer :: w, d, r

        w = 11; d = 4; iostat = 0
        if (size(v_list) >= 1) w = v_list(1)
        if (size(v_list) >= 2) d = v_list(2)

        write(fmt_string, "('(/, 9(ES', I0, '.', I0, ', 1X))')") w, d
        do r = 1, 9
            write(unit, fmt_string, iostat=iostat) dtv%vals(r, :)
        end do
    end subroutine print_ten_3D4O
#endif

end module muscle_tensor_3d4o