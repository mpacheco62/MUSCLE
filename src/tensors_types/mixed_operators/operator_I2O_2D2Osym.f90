module mod_operator_I2O_2D2Osym
    use, intrinsic :: iso_fortran_env
    use mod_iden_2O
    use mod_ten_2D2Osym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_2D2Osym_I2O
        module procedure sum_I2O_2D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D2Osym_I2O
        module procedure sub_I2O_2D2Osym
    end interface
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_2D2Osym_I2O
        module procedure mul_I2O_2D2Osym
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2O_2D2Osym
        module procedure ddot_2D2Osym_I2O
    end interface

    contains

    pure module function ddot_I2O_2D2Osym(I2, b) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        class(ten_2D2Osym), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_I2O_2D2Osym

    pure module function ddot_2D2Osym_I2O(b, I2) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        class(ten_2D2Osym), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_2D2Osym_I2O

    pure module function sum_I2O_2D2Osym(I2, a) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1D0
    end function sum_I2O_2D2Osym

    pure module function sum_2D2Osym_I2O(a, I2) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1D0
    end function sum_2D2Osym_I2O

    pure module function sub_I2O_2D2Osym(I2, a) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = -a%vals
        res%vals(1:3) = 1d0 + res%vals(1:3)
    end function sub_I2O_2D2Osym

    pure module function sub_2D2Osym_I2O(a, I2) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - 1D0
    end function sub_2D2Osym_I2O

    pure module function mul_I2O_2D2Osym(I2, a) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
    end function mul_I2O_2D2Osym

    pure module function mul_2D2Osym_I2O(a, I2) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
    end function mul_2D2Osym_I2O
end module mod_operator_I2O_2D2Osym