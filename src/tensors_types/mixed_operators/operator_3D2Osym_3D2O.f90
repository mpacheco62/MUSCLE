module mod_operator_3D2Osym_3D2O
    use, intrinsic :: iso_fortran_env
    use mod_ten_3D2Osym
    use mod_ten_3D2O
    implicit none
    private
    
    public :: operator(*)
    interface operator (*)
        module procedure dot_3D2Osym_3D2Osym
    end interface

!     public :: operator(-)
!     interface operator (-)
!         module procedure sub_3D2O_3D2Osym
!         module procedure sub_3D2Osym_3D2O
!     end interface

!     public :: operator(.ddot.)
!     interface operator (.ddot.)
!         module procedure ddot_3D2O_3D2Osym
!         module procedure ddot_3D2Osym_3D2O
!     end interface

!     public :: assignment (=)
!     interface assignment (=)
!         module procedure assign_3D2O_3D2Osym
!     end interface

contains
    pure function dot_3D2Osym_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        class(ten_3D2Osym), intent(in) :: b
        type(ten_3D2O) :: res
        ! xx
        res%vals(1) = a%vals(1)*b%vals(1) + a%vals(4)*b%vals(4) + a%vals(6)*b%vals(6)
        ! yx
        res%vals(2) = a%vals(4)*b%vals(1) + a%vals(2)*b%vals(4) + a%vals(5)*b%vals(6)
        ! zx
        res%vals(3) = a%vals(6)*b%vals(1) + a%vals(5)*b%vals(4) + a%vals(3)*b%vals(6)
        ! xy
        res%vals(4) = a%vals(1)*b%vals(4) + a%vals(4)*b%vals(2) + a%vals(6)*b%vals(5)
        ! yy
        res%vals(5) = a%vals(4)*b%vals(4) + a%vals(2)*b%vals(2) + a%vals(5)*b%vals(5)
        ! zy
        res%vals(6) = a%vals(6)*b%vals(4) + a%vals(5)*b%vals(2) + a%vals(3)*b%vals(5)
        ! xz
        res%vals(7) = a%vals(1)*b%vals(6) + a%vals(4)*b%vals(5) + a%vals(6)*b%vals(3)
        ! yz
        res%vals(8) = a%vals(4)*b%vals(6) + a%vals(2)*b%vals(5) + a%vals(5)*b%vals(3)
        ! zz
        res%vals(9) = a%vals(6)*b%vals(6) + a%vals(5)*b%vals(5) + a%vals(3)*b%vals(3)
    end function dot_3D2Osym_3D2Osym

end module mod_operator_3D2Osym_3D2O