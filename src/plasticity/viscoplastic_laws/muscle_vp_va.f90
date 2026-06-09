! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_vp_va
    use, intrinsic :: iso_fortran_env
    use muscle_vp_base
    implicit none
    private
    public :: VA_viscoplastic

    type, extends(Base_viscoplastic_law) :: VA_viscoplastic
        real(real64) :: B    
        real(real64) :: n    
        real(real64) :: B1   
        real(real64) :: B2   
        real(real64) :: m    
        real(real64) :: sig_u
    contains
        procedure :: flow_stress => flow_VA
        procedure :: dstress_dep => dstress_dep_VA
    end type VA_viscoplastic

contains

pure function flow_VA(self, ep, epd) result(res)
    class(VA_viscoplastic), intent(in) :: self
    real(real64), intent(in) :: ep, epd
    real(real64) :: res
    real(real64) :: factor_rate, ep_eff, epd_eff

    ep_eff = max(ep, 1.0d-10) 
    epd_eff = max(epd, 1.0d-12)

    if (self%m /= 0.0d0) then
        factor_rate = 1.0d0 + self%B1 * (epd_eff**(1.0d0/self%m)) - self%B2
    else
        factor_rate = 1.0d0 - self%B2
    end if

    res = self%B * (ep_eff**self%n) * factor_rate + self%sig_u
end function flow_VA

pure function dstress_dep_VA(self, ep, epd, dt) result(res)
    class(VA_viscoplastic), intent(in) :: self
    real(real64), intent(in) :: ep, epd, dt
    real(real64) :: res
    real(real64) :: ep_eff, epd_eff, factor_rate
    real(real64) :: dsig_dep, dsig_depd, d_factor_depd, hardening_part

    ep_eff = max(ep, 1.0d-10)
    epd_eff = max(epd, 1.0d-12)

    if (self%m /= 0.0d0) then
        factor_rate = 1.0d0 + self%B1 * (epd_eff**(1.0d0/self%m)) - self%B2
        d_factor_depd = (self%B1 / self%m) * (epd_eff**(1.0d0/self%m - 1.0d0))
    else
        factor_rate = 1.0d0 - self%B2
        d_factor_depd = 0.0d0
    end if

    hardening_part = self%B * (ep_eff**self%n)
    
    dsig_dep = self%B * self%n * (ep_eff**(self%n - 1.0d0)) * factor_rate
    
    dsig_depd = hardening_part * d_factor_depd
    
    res = dsig_dep + (dsig_depd / dt)
end function dstress_dep_VA

end module muscle_vp_va
