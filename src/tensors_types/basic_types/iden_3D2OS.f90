module mod_iden_3D2OS
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: iden_3D2OS  ! val*\delta_ij
        real(real64) :: val
        contains
            generic, public :: init => init_iden_3D2OS 
            procedure, private :: init_iden_3D2OS
    end type iden_3D2OS


    public :: operator(+)
    interface operator (+)
        module procedure sum_I3D2OS_I3D2OS
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure subU_I3D2OS
        module procedure sub_I3D2OS_I3D2OS
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D2OS_real64
        module procedure mul_real64_I3D2OS
        module procedure mul_I3D2OS_I3D2OS
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D2OS_real64
    end interface

    contains

    module subroutine init_iden_3D2OS(self, val)
        implicit none
        class(iden_3D2OS), intent(inout) :: self
        real(real64), intent(in) :: val
        self%val = val
    end subroutine init_iden_3D2OS

    pure module function sum_I3D2OS_I3D2OS(I2a, I2b) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: I2a, I2b
        type(iden_3D2OS) :: res
        res%val = I2a%val + I2b%val
    end function sum_I3D2OS_I3D2OS

    pure module function subU_I3D2OS(I2a) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: I2a
        type(iden_3D2OS) :: res
        res%val = -I2a%val
    end function subU_I3D2OS

    pure module function sub_I3D2OS_I3D2OS(I2a, I2b) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: I2a, I2b
        type(iden_3D2OS) :: res
        res%val = I2a%val - I2b%val
    end function sub_I3D2OS_I3D2OS

    pure module function mul_I3D2OS_I3D2OS(I2a, I2b) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: I2a, I2b
        type(iden_3D2OS) :: res
        res%val = I2a%val*I2b%val
    end function mul_I3D2OS_I3D2OS

    pure module function mul_I3D2OS_real64(IMod, a) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = IMod%val * a 
    end function mul_I3D2OS_real64

    pure module function mul_real64_I3D2OS(a, IMod) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = IMod%val * a 
    end function mul_real64_I3D2OS

    pure module function div_I3D2OS_real64(IMod, a) result(res)
        implicit none
        class(iden_3D2OS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = IMod%val/a 
    end function div_I3D2OS_real64

end module mod_iden_3D2OS