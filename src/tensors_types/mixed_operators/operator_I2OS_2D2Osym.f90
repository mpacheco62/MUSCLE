module mod_operator_I2OS_2D2Osym
    use, intrinsic :: iso_fortran_env
    use mod_iden_2OS
    use mod_ten_2D2Osym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_2D2Osym_I2OS
        module procedure sum_I2OS_2D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D2Osym_I2OS
        module procedure sub_I2OS_2D2Osym
    end interface
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2OS_2D2Osym
        module procedure ddot_2D2Osym_I2OS
    end interface

    contains

    pure function ddot_I2OS_2D2Osym(I2, b) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2
        class(ten_2D2Osym), intent(in) :: b
        real(real64) :: res
        res = I2%val*(b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_I2OS_2D2Osym

    pure function ddot_2D2Osym_I2OS(b, I2) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2
        class(ten_2D2Osym), intent(in) :: b
        real(real64) :: res
        res = I2%val*(b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_2D2Osym_I2OS

    pure function sum_I2OS_2D2Osym(I2, a) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_I2OS_2D2Osym

    pure function sum_2D2Osym_I2OS(a, I2) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_2D2Osym_I2OS

    pure function sub_I2OS_2D2Osym(I2, a) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = -a%vals
        res%vals(1:3) = I2%val + res%vals(1:3)
    end function sub_I2OS_2D2Osym

    pure function sub_2D2Osym_I2OS(a, I2) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - I2%val
    end function sub_2D2Osym_I2OS
end module mod_operator_I2OS_2D2Osym