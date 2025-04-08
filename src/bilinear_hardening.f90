module mod_bilinear_hardening
    use, intrinsic :: iso_fortran_env
    use mod_hardening_law
    implicit none

    type, extends(Base_hardening_law) :: Bilinear_hardening
        real(real64) :: y0, K
    contains
        procedure :: stress => stress_bilinear
        procedure :: dstress_dep => dstress_dep_bilinear
        procedure :: ddstress_ddep => ddstress_ddep_bilinear
    end type Bilinear_hardening

    contains
    pure function stress_bilinear(self, ep) result(res)
        implicit none
        Class(Bilinear_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = self%y0 + self%K*ep
    end function stress_bilinear

    pure function dstress_dep_bilinear(self, ep) result(res)
        implicit none
        Class(Bilinear_hardening), intent(in) :: self  ! parameters
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = self%K
    end function dstress_dep_bilinear

    pure function ddstress_ddep_bilinear(self, ep) result(res)
        implicit none
        Class(Bilinear_hardening), intent(in) :: self  ! parameters
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = 0
    end function ddstress_ddep_bilinear
end module