module mod_JohnsonCook_visco_hardening
    !! Module mod_JohnsonCook_visco_hardening
    !! ==========================
    !!
    !! Defines a concrete implementation for the Johnson-Cook model 
    !!
    !! This module provides the `JohnsonCook_visco_hardening` derived type, which represents the
    !! Johnson Cook model for hardening in visco-plasticity. It extends the abstract
    !! `Base_visco_hardening_law` type defined in `mod_JohnsonCook_visco_hardening_law`.
    !!
    !! The JohnsonCook law describes the evolution of the flow stress (yield stress) as a
    !! function of the equivalent plastic strain.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `JohnsonCook_visco_hardening`: Concrete type representing the JohnsonCook isotropic visco hardening law.
    !!     - Extends: `Base_visco_hardening_law`.
    !!     - Components:
    !!         - `k :: real(real64)`: Strength coefficient.
    !!         - `e0 :: real(real64)`: Initial strain offset.
    !!         - `n :: real(real64)`: Hardening exponent.
    !!     - Procedure: `stress(ep)` - Calculates the flow stress using the JohnsonCook equation.
    !!     - Procedure: `dstress_dep(ep)` - Calculates the first derivative of the flow stress
    !!       with respect to equivalent plastic strain.
    !!     - Procedure: `ddstress_ddep(ep)` - Calculates the second derivative of the flow stress
    !!       with respect to equivalent plastic strain.
    !!
    !! Mathematical Background
    !! -----------------------
    !!
    !! The flow stress (\(\sigma_{flow}\)) is calculated as a function of the equivalent plastic
    !! strain (\(\epsilon_p\)) using the JohnsonCook equation:
    !! \[ \sigma_{flow}(\epsilon_p) = K (\epsilon_0 + \epsilon_p)^n \]
    !! where \(K\) is the strength coefficient, \(\epsilon_0\) is the initial strain offset,
    !! and \(n\) is the hardening exponent.
    !!
    !! The first derivative (hardening modulus \(H\)) is:
    !! \[ H = \frac{d\sigma_{flow}}{d\epsilon_p} = K n (\epsilon_0 + \epsilon_p)^{n-1} \]
    !!
    !! The second derivative is:
    !! \[ \frac{d^2\sigma_{flow}}{d\epsilon_p^2} = K n (n-1) (\epsilon_0 + \epsilon_p)^{n-2} \]
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_JohnsonCook_visco_hardening_usage
    !!   use mod_JohnsonCook_visco_hardening
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(JohnsonCook_visco_hardening) :: material_JohnsonCook
    !!   real(real64) :: eq_plastic_strain, yield_stress, hardening_H, H_prime
    !!
    !!   ! Define JohnsonCook hardening parameters
    !!   material_JohnsonCook%k = 500.0D0  ! Strength coefficient (e.g., MPa)
    !!   material_JohnsonCook%e0 = 0.001D0 ! Initial strain offset
    !!   material_JohnsonCook%n = 0.2D0    ! Hardening exponent
    !!
    !!   ! Define an equivalent plastic strain value
    !!   eq_plastic_strain = 0.05D0
    !!
    !!   ! Calculate hardening properties at the given strain
    !!   yield_stress = material_JohnsonCook%stress(eq_plastic_strain)
    !!   hardening_H = material_JohnsonCook%dstress_dep(eq_plastic_strain)
    !!   H_prime = material_JohnsonCook%ddstress_ddep(eq_plastic_strain)
    !!
    !!   print *, "JohnsonCook Parameters:"
    !!   print *, "  K =", material_JohnsonCook%k
    !!   print *, "  e0 =", material_JohnsonCook%e0
    !!   print *, "  n =", material_JohnsonCook%n
    !!   print *, "At eq. plastic strain =", eq_plastic_strain
    !!   print *, "  Yield Stress =", yield_stress
    !!   print *, "  Hardening Modulus (H) =", hardening_H
    !!   print *, "  H' =", H_prime
    !!
    !! end program example_JohnsonCook_visco_hardening_usage
    !! ```
    !!
    !! For the base class definition see [[mod_visco_hardening_law]].

    use, intrinsic :: iso_fortran_env
    use mod_visco_hardening_law
    implicit none
    PRIVATE

    PUBLIC :: JohnsonCook_visco_hardening
    type, extends(Base_visco_hardening_law) :: JohnsonCook_visco_hardening
        !! Concrete type for JohnsonCook Isotropic Hardening Law.
        !! Implements the flow stress and its derivatives according to the Johnson-Cook equation:
        !! sigma_flow = K * (e0 + ep)^n
        real(real64) :: k  !! Strength coefficient (K)
        real(real64) :: e0 !! Initial strain offset (epsilon_0)
        real(real64) :: n  !! Hardening exponent (n)
    contains
        procedure :: stress => stress_JohnsonCook
            !! Calculates the flow stress sigma_flow = K*(e0 + ep)**n.
        procedure :: dstress_dep => dstress_dep_JohnsonCook
            !! Calculates the derivative d(sigma_flow)/d(ep).
        procedure :: ddstress_ddep => ddstress_ddep_JohnsonCook
            !! Calculates the second derivative d^2(sigma_flow)/d(ep)^2.
    end type JohnsonCook_visco_hardening

    contains
    pure function stress_JohnsonCook(self, ep) result(res)
        !! Calculates the flow stress using the JohnsonCook equation.
        implicit none
        class(JohnsonCook_visco_hardening), intent(in) :: self
            !! The JohnsonCook hardening law object containing parameters K, e0, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain.
        real(real64) :: res
            !! Output flow stress (yield stress) at the given `ep`.
        res = self%K*(self%e0 + ep)**self%n
    end function stress_JohnsonCook

    pure function dstress_dep_JohnsonCook(self, ep) result(res)
        !! Calculates the first derivative of the flow stress w.r.t. equivalent plastic strain (Hardening Modulus H).
        implicit none
        class(JohnsonCook_visco_hardening), intent(in) :: self
            !! The JohnsonCook hardening law object containing parameters K, e0, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain at which the derivative is evaluated.
        real(real64) :: res
            !! Output hardening modulus (d_stress / d_ep) at the given `ep`.
        res = self%K*self%n*(self%e0 + ep)**(self%n-1D0)
    end function dstress_dep_JohnsonCook

    pure function ddstress_ddep_JohnsonCook(self, ep) result(res)
        !! Calculates the second derivative of the flow stress w.r.t. equivalent visco-plastic strain.
        implicit none
        class(JohnsonCook_visco_hardening), intent(in) :: self
            !! The JohnsonCook hardening law object containing parameters K, e0, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain at which the second derivative is evaluated.
        real(real64) :: res
            !! Output second derivative (d^2_stress / d_ep^2) at the given `ep`.
        res = self%K*self%n*(self%n-1D0)*(self%e0 + ep)**(self%n-2D0)
    end function ddstress_ddep_JohnsonCook
end module