module mod_RK_viscoplastic
    !! Module mod_RK_viscoplastic
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
    !!     from `basis_viscoplastic_law`.
    !! 2. **Viscoplasticity:** Both the power-law exponent ($n_{val}$) and the scaling 
    !!     coefficient ($B_{val}$) are explicitly modified by the strain rate ($\dot{\epsilon}_p$).
    !!
    !! The final flow stress structure is:
    !!
    !! $$
    !! \sigma_{flow} = B(\dot{\epsilon}_p) \cdot (\epsilon_0 + \epsilon_p)^{n(\dot{\epsilon}_p)} + \sigma_{extra}(\dot{\epsilon}_p)
    !! $$
    use, intrinsic :: iso_fortran_env
    use mod_basis_viscoplastic_law
    implicit none
    private
    public :: RK_viscoplastic

    type, extends(basis_viscoplastic_law) :: RK_viscoplastic
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
    contains
        procedure :: flow_stress => flow_RK
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

end module mod_RK_viscoplastic
