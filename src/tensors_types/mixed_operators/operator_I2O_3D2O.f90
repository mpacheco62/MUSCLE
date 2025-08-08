module mod_operator_I2O_3D2O
    use, intrinsic :: iso_fortran_env
    use mod_iden_2O
    use mod_ten_3D2O
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_3D2O_I2O
        module procedure sum_I2O_3D2O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D2O_I2O
        module procedure sub_I2O_3D2O
    end interface
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_3D2O_I2O
        module procedure mul_I2O_3D2O
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2O_3D2O
        module procedure ddot_3D2O_I2O
    end interface

    contains

    pure module function ddot_I2O_3D2O(I2, b) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        class(ten_3D2O), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(5) + b%vals(9)
    end function ddot_I2O_3D2O

    pure module function ddot_3D2O_I2O(b, I2) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        class(ten_3D2O), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(5) + b%vals(9)
    end function ddot_3D2O_I2O

    pure module function sum_I2O_3D2O(I2, a) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1D0
    end function sum_I2O_3D2O

    pure module function sum_3D2O_I2O(a, I2) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) + 1D0
        res%vals(5) = res%vals(5) + 1D0
        res%vals(9) = res%vals(9) + 1D0
    end function sum_3D2O_I2O

    pure module function sub_I2O_3D2O(I2, a) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = -a%vals
        res%vals(1) = 1d0 + res%vals(1)
        res%vals(5) = 1d0 + res%vals(5)
        res%vals(9) = 1d0 + res%vals(9)
    end function sub_I2O_3D2O

    pure module function sub_3D2O_I2O(a, I2) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) - 1D0
        res%vals(5) = res%vals(5) - 1D0
        res%vals(9) = res%vals(9) - 1D0
    end function sub_3D2O_I2O

    pure module function mul_I2O_3D2O(I2, a) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
    end function mul_I2O_3D2O

    pure module function mul_3D2O_I2O(a, I2) result(res)
        implicit none
        class(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
    end function mul_3D2O_I2O

end module mod_operator_I2O_3D2O