! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco@usach.cl>

module muscle_hard_base
    !! Module muscle_hard_base
    !! ========================
    !!
    !! Defines the abstract base type for hardening laws used in plasticity models.
    !!
    !! This module provides an abstract derived type, `Base_hardening_laws`, which serves as
    !! a blueprint for implementing various isotropic hardening models (e.g., linear, power law,
    !! saturation hardening). It defines the essential interface that any concrete hardening
    !! law must provide: procedures to calculate the flow stress (yield stress), the hardening
    !! modulus (first derivative of flow stress w.r.t. equivalent plastic strain), and the
    !! second derivative of the flow stress.
    !!
    !! Concrete implementations of hardening laws should extend this abstract type and
    !! provide specific implementations for the deferred procedures `stress`, `dstress_dep`,
    !! and `ddstress_ddep`.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Abstract Derived Type:
    !!
    !! - `Base_hardening_laws`: Abstract base type for isotropic hardening laws.
    !!     - Deferred Procedure: `stress(ep)` - Interface for a function that computes
    !!       the current flow stress (yield stress) `res` (`real(real64)`) as a function of the
    !!       accumulated equivalent plastic strain `ep` (`real(real64)`).
    !!     - Deferred Procedure: `dstress_dep(ep)` - Interface for a function that computes
    !!       the hardening modulus `res` (`real(real64)`), which is the derivative of the flow
    !!       stress with respect to the equivalent plastic strain (\( \frac{d\sigma_{flow}}{d\epsilon_p} \)),
    !!       evaluated at `ep`.
    !!     - Deferred Procedure: `ddstress_ddep(ep)` - Interface for a function that computes
    !!       the second derivative `res` (`real(real64)`) of the flow stress with respect to the
    !!       equivalent plastic strain (\( \frac{d^2\sigma_{flow}}{d\epsilon_p^2} \)), evaluated at `ep`.
    !!
    !! Usage (Conceptual)
    !! ------------------
    !!
    !! Concrete hardening laws will extend `Base_hardening_laws`:
    !!
    !! ```fortran
    !! module mod_linear_hardening
    !!   use muscle_hard_base
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type, extends(Base_hardening_laws) :: LinearHardening
    !!     real(real64) :: initial_yield_stress = 0.0D0 ! sigma_y0
    !!     real(real64) :: hardening_modulus    = 0.0D0 ! H
    !!   contains
    !!     procedure :: stress => calculate_linear_stress
    !!     procedure :: dstress_dep => calculate_linear_dstress_dep
    !!     procedure :: ddstress_ddep => calculate_linear_ddstress_ddep
    !!   end type LinearHardening
    !!
    !! contains
    !!
    !!   pure function calculate_linear_stress(self, ep) result(res)
    !!     class(LinearHardening), intent(in) :: self
    !!     real(real64), intent(in) :: ep
    !!     real(real64) :: res
    !!     res = self%initial_yield_stress + self%hardening_modulus * ep
    !!   end function calculate_linear_stress
    !!
    !!   pure function calculate_linear_dstress_dep(self, ep) result(res)
    !!     class(LinearHardening), intent(in) :: self
    !!     real(real64), intent(in) :: ep ! ep is unused for linear hardening derivative
    !!     real(real64) :: res
    !!     res = self%hardening_modulus
    !!   end function calculate_linear_dstress_dep
    !!
    !!   pure function calculate_linear_ddstress_ddep(self, ep) result(res)
    !!     class(LinearHardening), intent(in) :: self
    !!     real(real64), intent(in) :: ep ! ep is unused
    !!     real(real64) :: res
    !!     res = 0.0D0 ! Second derivative is zero for linear hardening
    !!   end function calculate_linear_ddstress_ddep
    !!
    !! end module mod_linear_hardening
    !!
    !! program use_hardening_law
    !!   use mod_linear_hardening
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(LinearHardening) :: material_hardening
    !!   class(Base_hardening_laws), allocatable :: law_ptr
    !!   real(real64) :: eq_plastic_strain, current_yield, current_H, current_H_prime
    !!
    !!   ! Initialize material hardening parameters
    !!   material_hardening%initial_yield_stress = 250.0D0 ! MPa
    !!   material_hardening%hardening_modulus = 1000.0D0   ! MPa
    !!
    !!   ! Point the abstract pointer to the concrete implementation
    !!   allocate(law_ptr, source=material_hardening)
    !!
    !!   ! Define an equivalent plastic strain value
    !!   eq_plastic_strain = 0.01D0
    !!
    !!   ! Use the abstract interface to call the concrete methods
    !!   current_yield = law_ptr%stress(eq_plastic_strain)
    !!   current_H = law_ptr%dstress_dep(eq_plastic_strain)
    !!   current_H_prime = law_ptr%ddstress_ddep(eq_plastic_strain)
    !!
    !!   print *, "At eq. plastic strain =", eq_plastic_strain
    !!   print *, "  Current Yield Stress =", current_yield
    !!   print *, "  Hardening Modulus (H) =", current_H
    !!   print *, "  H' =", current_H_prime
    !!
    !!   deallocate(law_ptr)
    !!
    !! end program use_hardening_law
    !! ```
    use, intrinsic :: iso_fortran_env
    implicit None
    PRIVATE
    
    PUBLIC :: Base_hardening_laws
    type, abstract :: Base_hardening_laws
        !! Abstract Base Type for Isotropic Hardening Laws
        !! ===============================================
        !!
        !! Serves as the fundamental interface for all isotropic hardening laws within the library.
        !! Any concrete hardening model (e.g., linear, power law) must extend this type
        !! and provide implementations for the deferred procedures `stress`, `dstress_dep`,
        !! and `ddstress_ddep`. This allows for polymorphic handling of different hardening behaviors
        !! in plasticity models. The independent variable is typically the equivalent plastic strain.
        contains
            procedure(stress_interface), deferred :: stress
                !! Computes the current flow stress (yield stress) for a given equivalent plastic strain.
            procedure(dstress_dep_interface), deferred :: dstress_dep
                !! Computes the hardening modulus (d_stress / d_ep) for a given equivalent plastic strain.
            procedure(ddstress_ddep_interface), deferred :: ddstress_ddep
                !! Computes the second derivative (d^2_stress / d_ep^2) for a given equivalent plastic strain.
    end type Base_hardening_laws

    interface
        pure function stress_interface(self, ep) result(res)
            !! Interface for the `stress` procedure.
            !! Must be implemented by concrete subtypes of `Base_hardening_laws`.
            use, intrinsic :: iso_fortran_env
            import Base_hardening_laws
            class(Base_hardening_laws), intent(in) :: self
                !! The hardening law object.
            real(real64), intent(in) :: ep
                !! Input equivalent plastic strain.
            real(real64) :: res
                !! Output flow stress (yield stress) at the given `ep`.
        end function stress_interface

        pure function dstress_dep_interface(self, ep) result(res)
            !! Interface for the `dstress_dep` procedure.
            !! Must be implemented by concrete subtypes of `Base_hardening_laws`.
            use, intrinsic :: iso_fortran_env
            import Base_hardening_laws
            class(Base_hardening_laws), intent(in) :: self
                !! The hardening law object.
            real(real64), intent(in) :: ep
                !! Input equivalent plastic strain at which the derivative is evaluated.
            real(real64) :: res
                !! Output hardening modulus (d_stress / d_ep) at the given `ep`.
        end function dstress_dep_interface

        pure function ddstress_ddep_interface(self, ep) result(res)
            !! Interface for the `ddstress_ddep` procedure.
            !! Must be implemented by concrete subtypes of `Base_hardening_laws`.
            use, intrinsic :: iso_fortran_env
            import Base_hardening_laws
            class(Base_hardening_laws), intent(in) :: self
                !! The hardening law object.
            real(real64), intent(in) :: ep
                !! Input equivalent plastic strain at which the second derivative is evaluated.
            real(real64) :: res
                !! Output second derivative (d^2_stress / d_ep^2) at the given `ep`.
        end function ddstress_ddep_interface

    end interface
end module