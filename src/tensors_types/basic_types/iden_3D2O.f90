module mod_iden_3D2O
    use, intrinsic :: iso_fortran_env
    use mod_iden_3D2OS
    implicit none
    private

    type, public :: iden_3D2O  ! val*\delta_ij
    end type iden_3D2O


    public :: operator(+)
    interface operator (+)
        module procedure sum_I3D2O_I3D2O
        module procedure sum_I3D2O_I3D2OS
        module procedure sum_I3D2OS_I3D2O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure subU_I3D2O
        module procedure sub_I3D2O_I3D2O
        module procedure sub_I3D2O_I3D2OS
        module procedure sub_I3D2OS_I3D2O
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D2O_real64
        module procedure mul_real64_I3D2O
        module procedure mul_I3D2O_I3D2O
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D2O_real64
    end interface

contains

    pure module function sum_I3D2O_I3D2O(I2a, I2b) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2a, I2b
        type(iden_3D2OS) :: res
        res%val = 2D0
    end function sum_I3D2O_I3D2O

    pure module function sum_I3D2O_I3D2OS(I2a, I2b) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2a
        class(iden_3D2OS), intent(in) :: I2b
        type(iden_3D2OS) :: res
        res%val = 1D0 + I2b%val
    end function sum_I3D2O_I3D2OS

    pure module function sum_I3D2OS_I3D2O(I2a, I2b) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: I2a
        class(iden_3D2O), intent(in) :: I2b
        type(iden_3D2OS) :: res
        res%val = I2a%val + 1D0
    end function sum_I3D2OS_I3D2O

    pure module function subU_I3D2O(I2a) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2a
        type(iden_3D2OS) :: res
        res%val = -1D0
    end function subU_I3D2O

    pure module function sub_I3D2O_I3D2O(I2a, I2b) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2a, I2b
        type(iden_3D2OS) :: res
        res%val = 0D0
    end function sub_I3D2O_I3D2O

    pure module function sub_I3D2O_I3D2OS(I2a, I2b) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2a
        class(iden_3D2OS), intent(in) :: I2b
        type(iden_3D2OS) :: res
        res%val = 1D0 - I2b%val
    end function sub_I3D2O_I3D2OS

    pure module function sub_I3D2OS_I3D2O(I2a, I2b) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: I2a
        class(iden_3D2O), intent(in) :: I2b
        type(iden_3D2OS) :: res
        res%val = I2a%val - 1D0
    end function sub_I3D2OS_I3D2O

    pure module function mul_I3D2O_real64(I2, a) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = a 
    end function mul_I3D2O_real64

    pure module function mul_real64_I3D2O(a, I2) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = a 
    end function mul_real64_I3D2O

    pure module function mul_I3D2O_I3D2O(I2a, I2b) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2a, I2b
        type(iden_3D2O) :: res
    end function mul_I3D2O_I3D2O

    pure module function div_I3D2O_real64(I2, a) result(res)
        implicit none
        class(iden_3D2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = 1D0/a
    end function div_I3D2O_real64

end module mod_iden_3D2O