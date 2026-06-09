! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_vp_rk
    !! Module muscle_vp_rk
    !! ==========================
    !!
    !! Implements the **Rusinek-Klepaczko (RK)** Viscoplastic Flow Stress Law.
    !!
    !! The RK model is a widely used phenomenological thermo-viscoplastic constitutive model 
    !! designed to describe the mechanical response of materials (typically metals) over a wide 
    !! range of strain rates ($\dot{\epsilon}$) and temperatures. It effectively separates the 
    !! stress response into components dependent on strain, strain rate, and temperature.
    !!
    !! This implementation calculates the flow stress ($\sigma_{flow}$) based on the accumulated 
    !! plastic strain ($\epsilon_p$) and the equivalent plastic strain rate ($\dot{\epsilon}_p$).
    !!
    !! Key features of this implementation:
    !! 1. **Internal Hardening:** The hardening function is calculated internally based on 
    !!     a power-law type saturation, $B(\dot{\epsilon}_p) \cdot (\epsilon_0 + \epsilon_p)^{n(\dot{\epsilon}_p)}$.
    !!     Therefore, this model does **not** utilize the inherited `hard_law` pointer 
    !!     from `Base_viscoplastic_law`.
    !! 2. **Viscoplasticity:** Both the power-law exponent ($n_{val}$) and the scaling 
    !!     coefficient ($B_{val}$) are explicitly modified by the strain rate ($\dot{\epsilon}_p$).
    !!
    !! The final flow stress structure is:
    !!
    !! $$
    !! \sigma_{flow} = B(\dot{\epsilon}_p) \cdot (\epsilon_0 + \epsilon_p)^{n(\dot{\epsilon}_p)} + \sigma_{extra}(\dot{\epsilon}_p)
    !! $$
    use, intrinsic :: iso_fortran_env
    use muscle_vp_base
    implicit none
    private
    public :: RK_viscoplastic

    type, extends(Base_viscoplastic_law) :: RK_viscoplastic
        !! Rusinek-Klepaczko (RK) Viscoplastic Law Implementation
        !! ========================================================
        !!
        !! Extends the abstract base type to implement the full RK constitutive law.
        !!
        real(real64) :: B0     
        real(real64) :: ep0    
        real(real64) :: n0     
        real(real64) :: D2     
        real(real64) :: D1     
        real(real64) :: nu     
        real(real64) :: epdmin 
        real(real64) :: epdmax 
        real(real64) :: sig0   
        real(real64) :: m      
        real(real64) :: epmin = 1D-6 
        real(real64) :: epmax = 10000 

    contains
        procedure :: flow_stress => flow_RK
        procedure :: dstress_dep => dstress_dep_RK
    end type RK_viscoplastic

contains

    pure function flow_RK(self, ep, epd) result(res)
        !! flow_RK - Calculates the viscoplastic flow stress according to the Rusinek-Klepaczko (RK) model.
        !!
        !! Computes the flow stress by combining the rate-dependent power-law hardening term 
        !! with the extra rate-dependent stress component ($\sigma_{extra}$).
        !!
        class(RK_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        real(real64) :: n_val, B_val, sig_extra
        real(real64) :: factor_n, factor_sig, log_arg

        if (epd <= 0.0) then
            res = 0.0
            return
        end if

        factor_n = 1.0 - self%D2 * log10(epd / self%epdmin)
        if (factor_n > 0.0) then
            n_val = self%n0 * factor_n
        else
            n_val = 0.0
        end if

        log_arg = log10(self%epdmax / epd)
        B_val   = self%B0 * log_arg**(-self%nu)

        factor_sig = 1.0 - self%D1 * log10(self%epdmax / epd)
        if (factor_sig > 0.0) then
            sig_extra = self%sig0 * factor_sig**(self%m)
        else
            sig_extra = 0.0
        end if

        res = B_val * (self%ep0 + ep)**(n_val) + sig_extra
    end function flow_RK
    pure function dstress_dep_RK(self, ep, epd, dt) result(res)
        class(RK_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd, dt  ! ep = ep_old + dep
        real(real64) :: res
        
        real(real64) :: n_val, B_val, factor_sig, log_arg_B, log_arg_n, epd_tmp
        real(real64) :: dn_depd, dB_depd, dsig_extra_depd, inv_ln10
        real(real64) :: K_base, term1_deriv, term2_deriv

        epd_tmp = max(self%epmin, min(epd, self%epdmax * 0.99d0))
        inv_ln10 = 1.0d0 / log(10.0d0)
        
        log_arg_n = log10(epd_tmp / self%epdmin)
        n_val = self%n0 * max(0.0d0, (1.0d0 - self%D2 * log_arg_n))
        
        log_arg_B = max(log10(self%epdmax / epd_tmp), 0.01d0)
        B_val = self%B0 * (log_arg_B**(-self%nu))
        
        factor_sig = 1.0d0 - self%D1 * log10(self%epdmax / epd_tmp)

        K_base = self%ep0 + ep

        dn_depd = self%n0 * (-self%D2) * (inv_ln10 / (epd_tmp * dt))

        dB_depd = (self%B0 * self%nu * (log_arg_B**(-self%nu - 1.0d0))) * &
                  (inv_ln10 / (epd_tmp * dt))

        term1_deriv = (dB_depd * (K_base**n_val)) + &
                      (B_val * (K_base**n_val) * log(K_base) * dn_depd)

        if (factor_sig > 1.0d-10) then
            term2_deriv = self%sig0 * self%m * (factor_sig**(self%m - 1.0d0)) * &
                          (self%D1 * inv_ln10 / (epd_tmp * dt))
        else
            term2_deriv = 0.0d0
        end if

        res = term1_deriv + term2_deriv

    end function dstress_dep_RK
    

end module muscle_vp_rk
