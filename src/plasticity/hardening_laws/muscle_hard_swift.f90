! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco@usach.cl>

module muscle_hard_swift
    !! Module muscle_hard_swift
    !! ==========================
    !!
    !! Defines a concrete implementation for the Swift isotropic hardening law.
    !!
    !! This module provides the `Swift_hardening` derived type, which represents the
    !! Swift model for isotropic hardening in plasticity. It extends the abstract
    !! `Base_hardening_laws` type defined in `muscle_hard_base`.
    !!
    !! The Swift law describes the evolution of the flow stress (yield stress) as a
    !! function of the equivalent plastic strain.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `Swift_hardening`: Concrete type representing the Swift isotropic hardening law.
    !!     - Extends: `Base_hardening_laws`.
    !!     - Components:
    !!         - `k :: real(real64)`: Strength coefficient.
    !!         - `e0 :: real(real64)`: Initial strain offset.
    !!         - `n :: real(real64)`: Hardening exponent.
    !!     - Procedure: `stress(ep)` - Calculates the flow stress using the Swift equation.
    !!     - Procedure: `dstress_dep(ep)` - Calculates the first derivative of the flow stress
    !!       with respect to equivalent plastic strain.
    !!     - Procedure: `ddstress_ddep(ep)` - Calculates the second derivative of the flow stress
    !!       with respect to equivalent plastic strain.
    !!
    !! Mathematical Background
    !! -----------------------
    !!
    !! The flow stress (\(\sigma_{flow}\)) is calculated as a function of the equivalent plastic
    !! strain (\(\epsilon_p\)) using the Swift equation:
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
    !! program example_swift_hardening_usage
    !!   use muscle_hard_swift
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(Swift_hardening) :: material_swift
    !!   real(real64) :: eq_plastic_strain, yield_stress, hardening_H, H_prime
    !!
    !!   ! Define Swift hardening parameters
    !!   material_swift%k = 500.0D0  ! Strength coefficient (e.g., MPa)
    !!   material_swift%e0 = 0.001D0 ! Initial strain offset
    !!   material_swift%n = 0.2D0    ! Hardening exponent
    !!
    !!   ! Define an equivalent plastic strain value
    !!   eq_plastic_strain = 0.05D0
    !!
    !!   ! Calculate hardening properties at the given strain
    !!   yield_stress = material_swift%stress(eq_plastic_strain)
    !!   hardening_H = material_swift%dstress_dep(eq_plastic_strain)
    !!   H_prime = material_swift%ddstress_ddep(eq_plastic_strain)
    !!
    !!   print *, "Swift Parameters:"
    !!   print *, "  K =", material_swift%k
    !!   print *, "  e0 =", material_swift%e0
    !!   print *, "  n =", material_swift%n
    !!   print *, "At eq. plastic strain =", eq_plastic_strain
    !!   print *, "  Yield Stress =", yield_stress
    !!   print *, "  Hardening Modulus (H) =", hardening_H
    !!   print *, "  H' =", H_prime
    !!
    !! end program example_swift_hardening_usage
    !! ```
    !!
    !! For the base class definition see [[muscle_hard_base]].

    use, intrinsic :: iso_fortran_env
    use muscle_hard_base
    implicit none
    PRIVATE

    PUBLIC :: Swift_hardening
    type, extends(Base_hardening_laws) :: Swift_hardening
        !! Concrete type for Swift Isotropic Hardening Law.
        !! Implements the flow stress and its derivatives according to the Swift equation:
        !! sigma_flow = K * (e0 + ep)^n
        real(real64) :: k  !! Strength coefficient (K)
        real(real64) :: e0 !! Initial strain offset (epsilon_0)
        real(real64) :: n  !! Hardening exponent (n)
    contains
        procedure :: stress => stress_swift
            !! Calculates the flow stress sigma_flow = K*(e0 + ep)**n.
        procedure :: dstress_dep => dstress_dep_swift
            !! Calculates the derivative d(sigma_flow)/d(ep).
        procedure :: ddstress_ddep => ddstress_ddep_swift
            !! Calculates the second derivative d^2(sigma_flow)/d(ep)^2.
    end type Swift_hardening

    contains
    pure function stress_swift(self, ep) result(res)
        !! Calculates the flow stress using the Swift equation.
        implicit none
        class(Swift_hardening), intent(in) :: self
            !! The Swift hardening law object containing parameters K, e0, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain.
        real(real64) :: res
            !! Output flow stress (yield stress) at the given `ep`.
        res = self%K*(self%e0 + ep)**self%n
    end function stress_swift

    pure function dstress_dep_swift(self, ep) result(res)
        !! Calculates the first derivative of the flow stress w.r.t. equivalent plastic strain (Hardening Modulus H).
        implicit none
        class(Swift_hardening), intent(in) :: self
            !! The Swift hardening law object containing parameters K, e0, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain at which the derivative is evaluated.
        real(real64) :: res
            !! Output hardening modulus (d_stress / d_ep) at the given `ep`.
        res = self%K*self%n*(self%e0 + ep)**(self%n-1D0)
    end function dstress_dep_swift

    pure function ddstress_ddep_swift(self, ep) result(res)
        !! Calculates the second derivative of the flow stress w.r.t. equivalent plastic strain.
        implicit none
        class(Swift_hardening), intent(in) :: self
            !! The Swift hardening law object containing parameters K, e0, n.
        real(real64), intent(in) :: ep
            !! Input equivalent plastic strain at which the second derivative is evaluated.
        real(real64) :: res
            !! Output second derivative (d^2_stress / d_ep^2) at the given `ep`.
        res = self%K*self%n*(self%n-1D0)*(self%e0 + ep)**(self%n-2D0)
    end function ddstress_ddep_swift
end module