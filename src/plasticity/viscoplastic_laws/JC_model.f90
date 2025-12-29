module mod_JC_viscoplastic
    !! Module mod_JC_viscoplastic
    !! ===========================
    !!
    !! Implements the strain-rate dependent component of the **Johnson-Cook (JC) Viscoplastic Flow Stress Model**.
    !!
    !! This module defines the concrete derived type `JC_viscoplastic`, which extends the
    !! abstract base type `basis_viscoplastic_law`. The Johnson-Cook model defines the flow stress
    !! ($\sigma_{flow}$) as a product of three independent factors:
    !!
    !! 1. Hardening (strain dependency)
    !! 2. Strain Rate dependency
    !! 3. Temperature dependency (omitted in this basic implementation)
    !!
    !! The overall flow stress is computed using the associated hardening law ($\sigma_{hard}$)
    !! and the internal strain rate parameters ($C$ and $\dot{\epsilon}_{p0}$):
    !!
    !! $$
    !! \sigma_{flow} = \sigma_{hard}(\epsilon_p) \cdot \left[ 1 + C \ln \left( \frac{\dot{\epsilon}_p}{\dot{\epsilon}_{p0}} \right) \right]
    !! $$
    !!
    !! Where $\sigma_{hard}(\epsilon_p) = \text{self}\% \text{hard\_law}\% \text{stress}(ep)$.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `JC_viscoplastic`: Concrete type implementing the Johnson-Cook strain-rate dependency.
    !!     - Inherits Component: `hard_law` (from `basis_viscoplastic_law`) - The delegated hardening law ($\sigma_{hard}$).
    !!     - Component: `C :: real(real64)` - The Johnson-Cook strain rate sensitivity parameter.
    !!     - Component: `epdmax :: real(real64)` - The reference equivalent plastic strain rate ($\dot{\epsilon}_{p0}$).
    !!     - Implements Procedure: `flow_stress => flow_JC` - Computes the total flow stress including hardening and rate effects.
    !!
    !! Procedures
    !! ----------
    !!
    !! * `flow_JC(self, ep, epd)`:
    !!     - **Purpose:** Calculates the viscoplastic flow stress using the JC rate formulation.
    !!     - **Input Arguments:**
    !!       - `ep` (`real(real64)`): Equivalent plastic strain ($\epsilon_p$).
    !!       - `epd` (`real(real64)`): Equivalent plastic strain rate ($\dot{\epsilon}_p$).
    !!     - **Output:** `res` (`real(real64)`): The flow stress $\sigma_{flow}$.
    use, intrinsic :: iso_fortran_env
    use mod_basis_viscoplastic_law
    implicit none
    private
    public :: JC_viscoplastic

    type, extends(basis_viscoplastic_law) :: JC_viscoplastic
    !! Johnson-Cook Viscoplastic Law Implementation
    !! ===========================================
    !!
    !! Extends the `basis_viscoplastic_law` to implement the specific rate-dependent
    !! component of the Johnson-Cook model.
        real(real64) :: C     
        real(real64) :: epdmax
    contains
        procedure :: flow_stress => flow_JC
    end type JC_viscoplastic

contains

    pure function flow_JC(self, ep, epd) result(res)
    !! flow_JC - Calculates the flow stress using the Johnson-Cook rate term.
    !!
    !! Calculates the total viscoplastic flow stress:
    !! $\sigma_{flow} = \sigma_{hard} \cdot (1 + C \ln(\dot{\epsilon}_p / \dot{\epsilon}_{p0}))$.
        class(JC_viscoplastic), intent(in) :: self
        !! self The JC_viscoplastic object (containing C and epdmax)
        real(real64), intent(in) :: ep, epd
        !! ep: Equivalent plastic strain.
        !! epd: Equivalent plastic strain rate.
        real(real64) :: res
        !! Output flow stress $\sigma_{flow}$.
        real(real64) :: sigma0, factor_rate

        if (epd <= 0.0) then
            !! Handles zero or negative strain rate to avoid issues with log(0).
            res = 0.0
            return
        end if

        if (.not. associated(self%hard_law)) then
            !! Hardening law must be defined (associated).
            res = 0.0
            return
        end if

        ! 1. Get Hardening Stress ($\sigma_{hard}$)
        sigma0 = self%hard_law%stress(ep)
        ! 2. Calculate Strain Rate Factor
        factor_rate = 1.0 + self%C * log(epd / self%epdmax)
        ! 3. Combine to get Flow Stress
        res = sigma0 * factor_rate
    end function flow_JC

end module mod_JC_viscoplastic
