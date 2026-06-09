module mod_voce_m_hardening
    !! author: [Your Name/Initials]
    !! version: 1.0 - Initial documentation (Stress Function Only)
    !!
    !! Module mod_voce_m_hardening
    !! ==========================
    !!
    !! Implements the **Modified Voce Hardening Law** (also known as Voce-Linear or Voce-Saturated hardening).
    !!
    !! This law defines the hardening stress ($\sigma$) as a function of the equivalent plastic strain 
    !! ($\epsilon_p$) by combining a linear term and an exponential saturation term (Voce term). 
    !! This structure is commonly used to model materials that exhibit an initial rapid decrease in 
    !! the hardening rate, followed by a constant, non-zero hardening rate at large strains.
    !!
    !! La función de esfuerzo es:
    !!
    !! $$
    !! \sigma(\epsilon_p) = k \epsilon_p + q \left[ 1 - e^{-n \epsilon_p} \right]
    !! $$
    !! 
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `Voce_modified_hardening`: Concrete type extending `Base_hardening_laws`.
    !!     - Component: `k :: real(real64)` - The **linear hardening modulus** (asymptotic slope).
    !!     - Component: `q :: real(real64)` - The **saturation stress range** of the non-linear component.
    !!     - Component: `n :: real(real64)` - The **saturation rate exponent** (controls the speed of saturation).
    !!
    !! Procedures
    !! ----------
    !!
    !! * `stress => stress_voce_mod`: Calculates the hardening stress $\sigma(\epsilon_p)$.
    !!
    use, intrinsic :: iso_fortran_env
    use mod_hardening_laws
    implicit none
    private
    public :: Voce_modified_hardening

    type, extends(Base_hardening_laws) :: Voce_modified_hardening
        !! Voce Modified Hardening Law
        !! ===========================
        !! Implements the Voce hardening law with an additional linear term.
        real(real64) :: Sy
            !! Linear hardening modulus ($Sy$).
        real(real64) :: k 
            !! Linear hardening modulus ($k$).
        real(real64) :: q 
            !! Saturation stress range ($q$).
        real(real64) :: n 
            !! Saturation rate exponent ($n$).
        contains
            procedure :: stress => stress_voce_mod
            procedure :: dstress_dep   => dstress_dep_voce_mod
            procedure :: ddstress_ddep => ddstress_ddep_voce_mod
            !! Implements the stress calculation $\sigma(\epsilon_p)$.
    end type Voce_modified_hardening

contains

    pure function stress_voce_mod(self, ep) result(res)
       !! stress_voce_mod - Calculates the hardening stress $\sigma(\epsilon_p)$.
       implicit none
       class(Voce_modified_hardening), intent(in) :: self
           !! self The hardening object.
       real(real64), intent(in) :: ep
           !! ep Equivalent plastic strain ($\epsilon_p$).
       real(real64) :: res
           !! Output hardening stress $\sigma$.
       res = self%Sy + self%k*ep + self%q * (1.0 - exp(-self%n*ep))
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