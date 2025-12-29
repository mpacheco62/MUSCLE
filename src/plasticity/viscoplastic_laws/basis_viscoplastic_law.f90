module mod_viscoplastic_law
    !! Module mod_basis_viscoplastic_law
    !! =================================
    !!
    !! Defines the abstract base type for viscoplastic laws.
    !!
    !! This module provides an abstract derived type, `basis_viscoplastic_law`, which serves
    !! as a blueprint for implementing various viscoplastic models (e.g., Perzyna,
    !! Norton-Hoff, or Power Law creep models). It requires the **composition** of a
    !! hardening law object (`Base_hardening_law` from `mod_hardening_law`) to define
    !! the strain-dependent component of the flow stress.
    !!
    !! It defines the essential interface that any concrete viscoplastic law must provide:
    !! a procedure to calculate the viscoplastic flow stress.
    !!
    !! Concrete implementations of viscoplastic laws must extend this abstract type and
    !! provide a specific implementation for the **deferred** procedure `flow_stress`.
    !!
    !! Public Entities
    !! ------------------
    !!
    !! ### Abstract Derived Type:
    !!
    !! - `basis_viscoplastic_law`: Abstract base type for viscoplastic laws.
    !!     - Component: `hard_law` (`class(Base_hardening_law), pointer`) - Pointer to the
    !!       hardening law object that defines the static (strain-dependent) yield stress component.
    !!     - Deferred Procedure: `flow_stress(ep, epd)` - Interface for a function that computes
    !!       the **viscoplastic** flow stress `res` (`real(real64)`) as a function of the
    !!       accumulated equivalent plastic strain `ep` (`real(real64)`) and the
    !!       equivalent plastic strain rate `epd` (`real(real64)`).
    !!
    !! Usage (Conceptual)
    !! ------------------
    !!
    !! A concrete viscoplastic law (`MyViscoLaw`) would extend `basis_viscoplastic_law`
    !! and then use its `hard_law` component to retrieve the hardened yield stress
    !! (the static component) inside the `flow_stress` implementation.
    !!
    !! ```fortran
    !! module mod_my_visco_law
    !!   use mod_basis_viscoplastic_law
    !!   use mod_linear_hardening ! Example: depends on a concrete hardening law
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type, extends(basis_viscoplastic_law) :: MyViscoLaw
    !!     real(real64) :: viscosity_parameter = 1.0D0 ! eta
    !!   end type MyViscoLaw
    !!
    !! contains
    !!   ...
    !!   ! Implementation of flow_stress would call self%hard_law%stress(ep)
    !!   ...
    !! end module mod_my_visco_law
    !! ```  
    use, intrinsic :: iso_fortran_env
    use mod_hardening_law
    implicit none
    private
    public :: Base_viscoplastic_law

    type, abstract :: Base_viscoplastic_law
    !! Abstract Base Type for Viscoplastic Laws
    !! ========================================
    !!
    !! Defines the base type that combines the **hardening law** (static component) with
    !! **strain-rate dependent** behavior (viscoplastic component). This enables
    !! polymorphic handling of different viscoplastic models.
        class(Base_hardening_law), pointer :: hard_law => null()
        !! Pointer to the hardening law object (from mod_hardening_law).
        !! Must be explicitly allocated and associated before use.
    contains
        procedure(flow_stress_interface), deferred :: flow_stress
    end type Base_viscoplastic_law

    abstract interface
        pure function flow_stress_interface(self, ep, epd) result(res)
        !! Interface for the `flow_stress` procedure.
        !! Must be implemented by concrete subtypes of `basis_viscoplastic_law`.
            use, intrinsic :: iso_fortran_env
            import Base_viscoplastic_law
            class(Base_viscoplastic_law), intent(in) :: self
            real(real64), intent(in) :: ep    
            real(real64), intent(in) :: epd   
            real(real64) :: res               
        end function flow_stress_interface
    end interface

end module mod_viscoplastic_law
