module mod_iden_3D2OMod
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: iden_3D2Omod  ! val*\delta_ij
        real(real64) :: val
        contains
            generic, public :: init => init_iden_3D2OMod 
            procedure, private :: init_iden_3D2OMod
    end type iden_3D2Omod


    public :: operator(+)
    interface operator (+)
        module procedure sum_I3D2OMod_I3D2OMod
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure subU_I3D2OMod
        module procedure sub_I3D2OMod_I3D2OMod
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D2OMod_real64
        module procedure mul_real64_I3D2OMod
        module procedure mul_I3D2OMod_I3D2OMod
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D2OMod_real64
    end interface

    contains

    module subroutine init_iden_3D2OMod(self, val)
        implicit none
        class(iden_3D2Omod), intent(inout) :: self
        real(real64), intent(in) :: val
        self%val = val
    end subroutine init_iden_3D2OMod

    pure module function sum_I3D2OMod_I3D2OMod(I2a, I2b) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2a, I2b
        type(iden_3D2Omod) :: res
        res%val = I2a%val + I2b%val
    end function sum_I3D2OMod_I3D2OMod

    pure module function subU_I3D2OMod(I2a) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2a
        type(iden_3D2Omod) :: res
        res%val = -I2a%val
    end function subU_I3D2OMod

    pure module function sub_I3D2OMod_I3D2OMod(I2a, I2b) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2a, I2b
        type(iden_3D2Omod) :: res
        res%val = I2a%val - I2b%val
    end function sub_I3D2OMod_I3D2OMod

    pure module function mul_I3D2OMod_I3D2OMod(I2a, I2b) result(res)
        implicit none
        class(iden_3D2Omod), intent(in) :: I2a, I2b
        type(iden_3D2Omod) :: res
        res%val = I2a%val*I2b%val
    end function mul_I3D2OMod_I3D2OMod

    pure module function mul_I3D2OMod_real64(IMod, a) result(res)
        implicit none
        class(iden_3D2OMod), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D2OMod) :: res
        res%val = IMod%val * a 
    end function mul_I3D2OMod_real64

    pure module function mul_real64_I3D2OMod(a, IMod) result(res)
        implicit none
        class(iden_3D2OMod), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D2OMod) :: res
        res%val = IMod%val * a 
    end function mul_real64_I3D2OMod

    pure module function div_I3D2OMod_real64(IMod, a) result(res)
        implicit none
        class(iden_3D2OMod), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D2OMod) :: res
        res%val = IMod%val/a 
    end function div_I3D2OMod_real64

end module mod_iden_3D2OMod