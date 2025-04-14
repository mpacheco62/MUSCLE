! (\delta_ik\delta_jl + \delta_il\delta_jk)/2
module mod_iden_3D4O4T
    use, intrinsic :: iso_fortran_env
    use mod_iden_3D4O4TS
    implicit none
    private

    type, public :: iden_3D4O4T  ! (\delta_ik\delta_jl + \delta_il\delta_jk)/2
    end type iden_3D4O4T


    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D4O4T_real64
        module procedure mul_real64_I3D4O4T
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D4O4T_real64
    end interface

contains

    pure module function mul_I3D4O4T_real64(I2, a) result(res)
        implicit none
        class(iden_3D4O4T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D4O4TS) :: res
        res%val = a 
    end function mul_I3D4O4T_real64

    pure module function mul_real64_I3D4O4T(a, I2) result(res)
        implicit none
        class(iden_3D4O4T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D4O4TS) :: res
        res%val = a 
    end function mul_real64_I3D4O4T

    pure module function div_I3D4O4T_real64(I2, a) result(res)
        implicit none
        class(iden_3D4O4T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D4O4TS) :: res
        res%val = 1D0/a
    end function div_I3D4O4T_real64

end module mod_iden_3D4O4T