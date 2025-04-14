! \delta_ij\delta_kl
module mod_iden_3D4O3TS
    ! \delta_ij\delta_kl
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: iden_3D4O3TS  ! val*\delta_ij\delta_kl
        real(real64) :: val
    end type iden_3D4O3TS

    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D4O3TS_real64
        module procedure mul_real64_I3D4O3TS
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D4O3TS_real64
    end interface

    contains

    pure module function mul_I3D4O3TS_real64(IMod, a) result(res)
        implicit none
        class(iden_3D4O3TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TS) :: res
        res%val = IMod%val * a 
    end function mul_I3D4O3TS_real64

    pure module function mul_real64_I3D4O3TS(a, IMod) result(res)
        implicit none
        class(iden_3D4O3TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TS) :: res
        res%val = IMod%val * a 
    end function mul_real64_I3D4O3TS

    pure module function div_I3D4O3TS_real64(IMod, a) result(res)
        implicit none
        class(iden_3D4O3TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TS) :: res
        res%val = IMod%val/a 
    end function div_I3D4O3TS_real64

end module mod_iden_3D4O3TS