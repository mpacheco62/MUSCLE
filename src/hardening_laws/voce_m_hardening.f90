module mod_voce_m_hardening
    use, intrinsic :: iso_fortran_env
    use mod_hardening_law
    implicit none
    private
    public :: Voce_modified_hardening

    type, extends(Base_hardening_law) :: Voce_modified_hardening
        real(real64) :: k 
        real(real64) :: q 
        real(real64) :: n 
    contains
        procedure :: stress        => stress_voce_mod
        procedure :: dstress_dep   => dstress_dep_voce_mod
        procedure :: ddstress_ddep => ddstress_ddep_voce_mod
    end type Voce_modified_hardening

contains

    pure function stress_voce_mod(self, ep) result(res)
        class(Voce_modified_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = self%k*ep + self%q * (1.0 - exp(-self%n*ep))
    end function stress_voce_mod

    pure function dstress_dep_voce_mod(self, ep) result(res)
        class(Voce_modified_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = self%k + self%q*self%n*exp(-self%n*ep)
    end function dstress_dep_voce_mod

    pure function ddstress_ddep_voce_mod(self, ep) result(res)
        class(Voce_modified_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = -self%q*self%n*self%n*exp(-self%n*ep)
    end function ddstress_ddep_voce_mod

end module mod_voce_m_hardening
