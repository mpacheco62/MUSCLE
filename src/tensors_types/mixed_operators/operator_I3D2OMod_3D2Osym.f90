module mod_operator_I3D2OMod_3D2Osym
    use, intrinsic :: iso_fortran_env
    use mod_iden_3D2OMod
    use mod_ten_3D2Osym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_3D2Osym_I3D2OMod
        module procedure sum_I3D2OMod_3D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D2Osym_I3D2OMod
        module procedure sub_I3D2OMod_3D2Osym
    end interface
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I3D2OMod_3D2Osym
        module procedure ddot_3D2Osym_I3D2OMod
    end interface

    contains

    pure module function ddot_I3D2OMod_3D2Osym(I2, b) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2
        class(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = I2%val*(b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_I3D2OMod_3D2Osym

    pure module function ddot_3D2Osym_I3D2OMod(b, I2) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2
        class(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = I2%val*(b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_3D2Osym_I3D2OMod

    pure module function sum_I3D2OMod_3D2Osym(I2, a) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_I3D2OMod_3D2Osym

    pure module function sum_3D2Osym_I3D2OMod(a, I2) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_3D2Osym_I3D2OMod

    pure module function sub_I3D2OMod_3D2Osym(I2, a) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = -a%vals
        res%vals(1:3) = I2%val + res%vals(1:3)
    end function sub_I3D2OMod_3D2Osym

    pure module function sub_3D2Osym_I3D2OMod(a, I2) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - I2%val
    end function sub_3D2Osym_I3D2OMod
end module mod_operator_I3D2OMod_3D2Osym