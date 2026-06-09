! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_vp_mrk
    !! Module muscle_vp_mrk
    !! ===========================
    !!
    !! Implements the **Modified Rusinek-Klepaczko (MRK)** Viscoplastic Flow Stress Law.
    !!
    !! **CONTEXT NOTE:** The MRK model is a complex constitutive formulation, often referred to as a 
    !! **Modified Rusinek-Klepaczko Model**, developed as an extension of the original RK model 
    !! to better capture the nonlinear rate and strain dependence in materials like polymers and 
    !! glassy systems under high strain rates.
    !!
    !! This module defines the flow stress ($\sigma_{flow}$) based on the accumulated plastic
    !! strain ($\epsilon_p$) and the equivalent plastic strain rate ($\dot{\epsilon}_p$).
    !!
    !! The final flow stress is calculated as: $\sigma_{flow} = \sigma_{a} + \sigma_u$, where
    !! $\sigma_u$ is the ultimate stress (or back stress) and $\sigma_{a}$ is the rate-modified
    !! hardening component.
    use, intrinsic :: iso_fortran_env
    use muscle_vp_base
    implicit none
    private
    public :: MRK_viscoplastic

    type, extends(Base_viscoplastic_law) :: MRK_viscoplastic
    !! MRK Viscoplastic Law Implementation (Modified Rusinek-Klepaczko)
    !! ===============================================================
    !!
    !! Implements the complex viscoplastic model MRK. This model calculates its own
    !! hardening internally, hence it does not rely on the inherited `hard_law` pointer
    !! for the main hardening component.
        real(real64) :: B01    
        real(real64) :: B02    
        real(real64) :: epdmax 
        real(real64) :: nu1    
        real(real64) :: nu2    
        real(real64) :: n0     
        real(real64) :: D2     
        real(real64) :: epdmin 
        real(real64) :: chi1  
        real(real64) :: chi2  
        real(real64) :: sig_u 
    contains
        procedure :: flow_stress => flow_MRK
        procedure :: dstress_dep => dstress_dep_MRK
    end type MRK_viscoplastic

contains

    pure function flow_MRK(self, ep, epd) result(res)
        !! flow_MRK - Calculates the viscoplastic flow stress according to the MRK model.
        !!
        !! Computes $\sigma_{flow} = (\text{Rate\_Factor})^{\frac{1}{\chi_2}} \cdot \sigma_{0,a}(\epsilon_p, \dot{\epsilon}_p) + \sigma_u$.
        class(MRK_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        real(real64) :: n_val, B1_val, B2_val
        real(real64) :: sig0_a, sig_a
        real(real64) :: factor_n, factor_rate, log_arg

        if (epd <= 0.0) then
            res = self%sig_u
            return
        end if

        factor_n = 1.0 - self%D2 * log10(epd / self%epdmin)
        if (factor_n > 0.0) then
            n_val = self%n0 * factor_n
        else
            n_val = 0.0
        end if

        log_arg = log10(self%epdmax / epd)
        B1_val  = self%B01 * log_arg**(-self%nu1)
        B2_val  = self%B02 * log_arg**(-self%nu2)

        sig0_a = B1_val*ep + B2_val*(1.0 - exp(-ep*n_val))

        factor_rate = 1.0 - self%chi1 * log10(self%epdmax / epd)
        if (factor_rate > 0.0) then
            sig_a = factor_rate**(1.0/self%chi2) * sig0_a
        else
            sig_a = 0.0
        end if

        res = sig_a + self%sig_u
    end function flow_MRK

    pure function dstress_dep_MRK(self, ep, epd, dt) result(res)
        class(MRK_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd, dt
        real(real64) :: res
        real(real64) :: n_val, B1_val, B2_val, sig0_a, factor_rate
        real(real64) :: dsig0_dep, dsig0_depd, dfactor_depd, dn_depd, dB1_depd, dB2_depd
        real(real64) :: log_arg, inv_ln10, epd_eff

        if (epd <= 1.0d-12) then
            res = 0.0d0 
            return
        end if

        epd_eff = epd
        inv_ln10 = 1.0d0 / log(10.0d0)
        log_arg = max(1.0d-10, log10(self%epdmax / epd_eff))

       
        n_val  = self%n0 * max(0.0d0, 1.0d0 - self%D2 * log10(epd_eff / self%epdmin))
        B1_val = self%B01 * log_arg**(-self%nu1)
        B2_val = self%B02 * log_arg**(-self%nu2)
        sig0_a = B1_val * ep + B2_val * (1.0d0 - exp(-ep * n_val))
        factor_rate = max(1.0d-10, 1.0d0 - self%chi1 * log_arg) 
       
        dsig0_dep = B1_val + B2_val * n_val * exp(-ep * n_val)
        
        dn_depd  = -self%n0 * self%D2 * inv_ln10 / epd_eff
        dB1_depd = self%nu1 * B1_val / (log_arg * epd_eff * log(10.0d0))
        dB2_depd = self%nu2 * B2_val / (log_arg * epd_eff * log(10.0d0))
        
        dsig0_depd = dB1_depd * ep + dB2_depd * (1.0d0 - exp(-ep * n_val)) &
                     + B2_val * (ep * exp(-ep * n_val) * dn_depd)

        dfactor_depd = (1.0d0 / self%chi2) * (factor_rate**(1.0d0 / self%chi2 - 1.0d0)) &
                       * (self%chi1 * inv_ln10 / epd_eff)

        res = (factor_rate**(1.0d0 / self%chi2)) * dsig0_dep + &
              ((dfactor_depd * sig0_a) + (factor_rate**(1.0d0 / self%chi2)) * dsig0_depd) / dt
    end function dstress_dep_MRK
end module muscle_vp_mrk
