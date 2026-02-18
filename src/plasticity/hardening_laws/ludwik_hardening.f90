module mod_ludwik_hardening
    !! Module mod_ludwik_hardening
    !! ==========================
    !!
    !! Defines a concrete implementation for the Ludwik isotropic hardening law.
    !!
    !! This module provides the `ludwik_hardening` derived type, which represents the 
    !! Ludwik model for isotropic hardening in plasticity. It extends the abstract    
    !! `Base_hardening_law` type defined in `mod_hardening_law`.                     
    !!
    !! The Ludwik law describes the evolution of the flow stress (yield stress) as a
    !! function of the equivalent plastic strain.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ludwik_hardening`: Concrete type representing the Ludwik isotropic hardening law.
    !!     - Extends: `Base_hardening_law`.
    !!     - Components:
    !!         - `sigma0 :: real(real64)`: Initial yield stress.
    !!         - `k :: real(real64)`: Strength coefficient.
    !!         - `n :: real(real64)`: Hardening exponent.
    !!     - Procedure: `stress(ep)` - Calculates the flow stress using the Ludwik equation.
    !!     - Procedure: `dstress_dep(ep)` - Calculates the first derivative of the flow stress
    !!       with respect to equivalent plastic strain.
    !!     - Procedure: `ddstress_ddep(ep)` - Calculates the second derivative of the flow stress
    !!       with respect to equivalent plastic strain.
    !!
    !! Mathematical Background
    !! -----------------------
    !!
    !! The flow stress ($\sigma_{flow}$) is calculated as a function of the equivalent plastic
    !! strain ($\epsilon_p$) using the Ludwik equation:
    !! $$ \sigma_{flow}(\epsilon_p) = \sigma_0 + K ( \epsilon_p)^n $$
    !! where $\sigma_0$ is the initial yield stress, $K$ is the strength coefficient,
    !!  and $n$ is the hardening exponent.
    !!
    !! The first derivative (hardening modulus $H$) is:
    !! $$ H = \frac{d\sigma_{flow}}{d\epsilon_p} = K n (\epsilon_p)^{n-1} $$
    !!
    !! The second derivative is:
    !! $$ \frac{d^2\sigma_{flow}}{d\epsilon_p^2} = K n (n-1) (\epsilon_p)^{n-2} $$
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ludwik_hardening_usage
    !!   use mod_ludwik_hardening
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ludwik_hardening) :: material_ludwik
    !!   real(real64) :: eq_plastic_strain, yield_stress, hardening_H, H_prime
    !!
    !!   ! Define Ludwik hardening parameters
    !!   material_ludwik%sigma0 = 250.0D0 ! Initial yield stress
    !!   material_ludwik%k = 500.0D0      ! Strength coefficient (K)
    !!   material_ludwik%n = 0.2D0        ! Hardening exponent
    !!
    !!   ! Define an equivalent plastic strain value
    !!   eq_plastic_strain = 0.05D0
    !!
    !!   ! Calculate hardening properties
    !!   yield_stress = material_ludwik%stress(eq_plastic_strain)
    !!   hardening_H = material_ludwik%dstress_dep(eq_plastic_strain)
    !!   H_prime = material_ludwik%ddstress_ddep(eq_plastic_strain)
    !!
    !!   print *, "Ludwik Parameters:"
    !!   print *, "  sigma0 =", material_ludwik%sigma0
    !!   print *, "  K =", material_ludwik%k
    !!   print *, "  n =", material_ludwik%n
    !! end program example_ludwik_hardening_usage
    !! ```
    !! For the base class definition see [[mod_hardening_law]].

    use, intrinsic :: iso_fortran_env
    use mod_hardening_law
    implicit none
    PRIVATE

    PUBLIC :: ludwik_hardening
    type, extends(Base_hardening_law) :: ludwik_hardening
        !! Concrete type for Ludwik Isotropic Hardening Law.
        !! Implements the flow stress and its derivatives according to the Ludwik equation:
        !! sigma_flow = sigma0 + K * (ep)^n
        real(real64) :: sigma0 !! Initial yield stress
        real(real64) :: k      !! Strength coefficient (K)
        real(real64) :: n      !! Hardening exponent (n)
    contains
        procedure :: stress => stress_ludwik
            !! Calculates the flow stress sigma_flow = sigma0 + K*(ep)**n.
        procedure :: dstress_dep => dstress_dep_ludwik
            !! Calculates the derivative d(sigma_flow)/d(ep).
        procedure :: ddstress_ddep => ddstress_ddep_ludwik
            !! Calculates the second derivative d^2(sigma_flow)/d(ep)^2.
    end type ludwik_hardening

    contains
    pure function stress_ludwik(self, ep) result(res)
            !! Calculates the flow stress using the Ludwik equation.
        implicit none
            class(ludwik_hardening), intent(in) :: self
                !! The Ludwik hardening law object containing parameters sigma0, K, e0, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain.
        real(real64) :: res
            !! Output flow stress (yield stress) at the given `ep`.
            res = self%sigma0 + self%k*(ep)**self%n
        end function stress_ludwik

    pure function dstress_dep_ludwik(self, ep) result(res)
            !! Calculates the first derivative of the flow stress w.r.t. equivalent plastic strain (Hardening Modulus H).
        implicit none
            class(ludwik_hardening), intent(in) :: self
                !! The Ludwik hardening law object containing parameters sigma0, K, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain at which the derivative is evaluated.
        real(real64) :: res
            !! Output hardening modulus (d_stress / d_ep) at the given `ep`.
            res = self%k*self%n*(ep+1.0D-10)**(self%n-1D0)
        end function dstress_dep_ludwik

    pure function ddstress_ddep_ludwik(self, ep) result(res)
            !! Calculates the second derivative of the flow stress w.r.t. equivalent plastic strain.
        implicit none
            class(ludwik_hardening), intent(in) :: self
                !! The Ludwik hardening law object containing parameters sigma0, K, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain at which the second derivative is evaluated.
        real(real64) :: res
            !! Output second derivative (d^2_stress / d_ep^2) at the given `ep`.
            res = self%k*self%n*(self%n-1D0)*(ep+1.0D-10)**(self%n-2D0)
        end function ddstress_ddep_ludwik
    end module