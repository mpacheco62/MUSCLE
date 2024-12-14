! \delta_ij\delta_kl
module mod_iden_3D4O3TMod
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: iden_3D4O3TMod  ! val*\delta_ij
        real(real64) :: val
    end type iden_3D4O3TMod

    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D4O3TMod_real64
        module procedure mul_real64_I3D4O3TMod
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D4O3TMod_real64
    end interface

    contains

    pure module function mul_I3D4O3TMod_real64(IMod, a) result(res)
        implicit none
        class(iden_3D4O3TMod), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TMod) :: res
        res%val = IMod%val * a 
    end function mul_I3D4O3TMod_real64

    pure module function mul_real64_I3D4O3TMod(a, IMod) result(res)
        implicit none
        class(iden_3D4O3TMod), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TMod) :: res
        res%val = IMod%val * a 
    end function mul_real64_I3D4O3TMod

    pure module function div_I3D4O3TMod_real64(IMod, a) result(res)
        implicit none
        class(iden_3D4O3TMod), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TMod) :: res
        res%val = IMod%val/a 
    end function div_I3D4O3TMod_real64

end module mod_iden_3D4O3TMod