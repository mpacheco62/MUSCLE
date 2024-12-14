module mod_ten_3D2Osym
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: ten_3D2Osym
        !! Type of a second order (O2), symmetric (sym) whith three dimensions (3D)
        real(real64), dimension(6) :: vals
        contains
            generic, public :: init => init_ten_3D2Osym, init2_ten_3D2Osym 
            procedure, private :: init_ten_3D2Osym, init2_ten_3D2Osym
    end type ten_3D2Osym

    public :: operator(.isequal.)
    interface operator (.isequal.)
        module procedure isequal_3D2Osym
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D2Osym
        module procedure subU_3D2Osym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_3D2Osym
        module procedure mul_3D2Osym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_3D2Osym_real64
    end interface

    public :: operator(.dev.)
    interface operator (.dev.)
        module procedure dev_3D2Osym
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D2Osym_3D2Osym
    end interface

    ! interface operator (.tdot.)
    !     module procedure tdot_3D2Osym_3D2Osym
    ! end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure ten_3D2Osym_real64_assign
    end interface

contains

    pure module subroutine ten_3D2Osym_real64_assign(a, b)
        implicit none
        type(ten_3D2Osym), intent(out) :: a
        real(real64), intent(in) :: b
        a%vals = b
    end subroutine

    module subroutine init_ten_3D2Osym(self, vals)
        implicit none
        class(ten_3D2Osym), intent(inout) :: self
        real(real64), intent(in) :: vals(6)
        self%vals = vals
    end subroutine

    module subroutine init2_ten_3D2Osym(self, xx, yy, zz, xy, yz, xz)
        ! voigt notation used: 11, 22, 33, 12, 23, 13
        implicit none
        class(ten_3D2Osym), intent(inout) :: self
        real(real64), intent(in) :: xx, yy, zz, xy, yz, xz
        self%vals = (/xx, yy, zz, xy, yz, xz/)
    end subroutine

    pure module function isequal_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        logical :: res

        real(real64), parameter :: EPS=1e-7, EPS_ABS=1e-30
        real(real64) :: norm_a, norm_b, norm_max, norm

        norm_a =   abs(a%vals(1)) + abs(a%vals(2)) + abs(a%vals(3)) &
                 + 2*abs(a%vals(4)) + 2*abs(a%vals(5)) + 2*abs(a%vals(6))
        norm_b =   abs(b%vals(1)) + abs(b%vals(2)) + abs(b%vals(3)) &
                 + 2*abs(b%vals(4)) + 2*abs(b%vals(5)) + 2*abs(b%vals(6))

        norm_max = max(max(norm_a, norm_b), EPS_ABS)

        norm =     abs(a%vals(1)-b%vals(1)) +   abs(a%vals(2)-b%vals(2)) +   abs(a%vals(3)-b%vals(3)) &
               + 2*abs(a%vals(4)-b%vals(4)) + 2*abs(a%vals(5)-b%vals(5)) + 2*abs(a%vals(6)-b%vals(6))

        
        if (norm/norm_max .gt. EPS) res=.false.
        if (norm/norm_max .le. EPS) res=.true.
        
    end function isequal_3D2Osym

    pure module function sum_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        type(ten_3D2Osym) :: res
        res%vals = a%vals + b%vals
    end function sum_3D2Osym

    pure module function sub_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        type(ten_3D2Osym) :: res
        res%vals = a%vals - b%vals
    end function sub_3D2Osym

    pure module function subU_3D2Osym(a) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = -a%vals
    end function subU_3D2Osym

    pure module function mul_real64_3D2Osym(a, b) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D2Osym), intent(in) :: b
        type(ten_3D2Osym) :: res
        res%vals = a * b%vals
    end function mul_real64_3D2Osym

    pure module function mul_3D2Osym_real64(a, b) result(res)
        implicit none
        real(real64), intent(in) :: b
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals =  a%vals * b
    end function mul_3D2Osym_real64

    pure module function div_3D2Osym_real64(a, b) result(res)
        implicit none
        real(real64), intent(in) :: b
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals/b
    end function div_3D2Osym_real64

    pure module function ddot_3D2Osym_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        real(real64) :: res
        res =   a%vals(1)*b%vals(1) &
                             + a%vals(2)*b%vals(2) &
                             + a%vals(3)*b%vals(3) &
                             + 2*a%vals(4)*b%vals(4) &
                             + 2*a%vals(5)*b%vals(5) &
                             + 2*a%vals(6)*b%vals(6)
    end function ddot_3D2Osym_3D2Osym

    pure module function dev_3D2Osym(a) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        real(real64) :: hydro
        hydro = (a%vals(1) + a%vals(2) + a%vals(3))/3D0
        res%vals(1:3) = a%vals(1:3) - hydro
        res%vals(4:6) = a%vals(4:6)
    end function dev_3D2Osym

end module mod_ten_3D2Osym