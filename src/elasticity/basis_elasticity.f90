module mod_base_elasticity
    !! Module mod_base_elasticity
    !! ==========================
    !!
    !! Defines the abstract base type for elasticity models.
    !!
    !! This module provides an abstract derived type, `Base_elasticity`, which serves as
    !! a blueprint for implementing various elasticity models (e.g., isotropic, anisotropic).
    !! It defines the essential interface that any concrete elasticity model must provide:
    !! procedures to calculate the stress tensor and the tangent modulus (material stiffness tensor)
    !! given a strain tensor.
    !!
    !! Concrete implementations of elasticity models should extend this abstract type and
    !! provide specific implementations for the deferred procedures `stress` and `dstress_dstrain`.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Abstract Derived Type:
    !!
    !! - `Base_elasticity`: Abstract base type for elasticity models.
    !!     - Deferred Procedure: `stress(strain)` - Interface for a function that computes
    !!       the second-order stress tensor (`ten_3D2Osym`) from a given second-order
    !!       strain tensor (`ten_3D2Osym`).
    !!     - Deferred Procedure: `dstress_dstrain(strain)` - Interface for a function that
    !!       computes the fourth-order tangent modulus tensor (`ten_3D4O3sym`) (material stiffness \( \frac{\partial \sigma}{\partial \epsilon} \))
    !!       at a given second-order strain tensor (`ten_3D2Osym`).
    !!
    !! Usage (Conceptual)
    !! ------------------
    !!
    !! Concrete elasticity models will extend `Base_elasticity`:
    !!
    !! ```fortran
    !! module mod_isotropic_elasticity
    !!   use mod_base_elasticity
    !!   use tensors_types
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type, extends(Base_elasticity) :: IsotropicElasticity
    !!     real(real64) :: E_mod = 0.0D0  ! Young's Modulus
    !!     real(real64) :: nu    = 0.0D0  ! Poisson's Ratio
    !!   contains
    !!     procedure :: stress => calculate_iso_stress
    !!     procedure :: dstress_dstrain => calculate_iso_stiffness
    !!   end type IsotropicElasticity
    !!
    !! contains
    !!
    !!   pure function calculate_iso_stress(self, strain) result(res)
    !!     class(IsotropicElasticity), intent(in) :: self
    !!     class(ten_3D2Osym), intent(in) :: strain
    !!     type(ten_3D2Osym) :: res
    !!     ! ... implementation using self%E_mod, self%nu ...
    !!   end function calculate_iso_stress
    !!
    !!   pure function calculate_iso_stiffness(self, strain) result(res)
    !!     class(IsotropicElasticity), intent(in) :: self
    !!     class(ten_3D2Osym), intent(in) :: strain
    !!     type(ten_3D4O3sym) :: res
    !!     ! ... implementation using self%E_mod, self%nu ...
    !!   end function calculate_iso_stiffness
    !!
    !! end module mod_isotropic_elasticity
    !!
    !! program use_elasticity_model
    !!   use mod_isotropic_elasticity
    !!   use tensors_types
    !!   implicit none
    !!
    !!   type(IsotropicElasticity) :: material
    !!   class(Base_elasticity), allocatable :: model_ptr
    !!   type(ten_3D2Osym) :: current_strain, resulting_stress
    !!   type(ten_3D4O3sym) :: tangent_modulus
    !!
    !!   ! Initialize material properties
    !!   material%E_mod = 210000.0D0
    !!   material%nu = 0.3D0
    !!
    !!   ! Point the abstract pointer to the concrete implementation
    !!   allocate(model_ptr, source=material)
    !!
    !!   ! Initialize strain
    !!   call current_strain%init(0.001, 0.0, 0.0, 0.0005, 0.0, 0.0)
    !!
    !!   ! Use the abstract interface to call the concrete methods
    !!   resulting_stress = model_ptr%stress(current_strain)
    !!   tangent_modulus = model_ptr%dstress_dstrain(current_strain)
    !!
    !!   print *, "Resulting Stress (xx):", resulting_stress%xx()
    !!   print *, "Tangent Modulus (1111):", tangent_modulus%vals(1)
    !!
    !!   deallocate(model_ptr)
    !!
    !! end program use_elasticity_model
    !! ```
    !!
    !! For tensor type definitions see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    ! use tensors_types
    implicit None
    PRIVATE

    PUBLIC :: Base_elasticity
    type, abstract :: Base_elasticity
        !! Abstract Base Type for Elasticity Models
        !! ========================================
        !!
        !! Serves as the fundamental interface for all elasticity models within the library.
        !! Any concrete elasticity model (e.g., isotropic, anisotropic) must extend this type
        !! and provide implementations for the deferred procedures `stress` and `dstress_dstrain`.
        !! This allows for polymorphic handling of different material models.
        contains
            procedure(stress_interface_3D), deferred :: stress_3D
                !! Computes the stress tensor for a given strain tensor.
            procedure(stress_interface_2D), deferred :: stress_2D
                !! Computes the stress tensor for a given strain tensor.
            procedure(dstress_dstrain_interface_3D), deferred :: dstress_dstrain_3D
                !! Computes the tangent modulus (stiffness tensor) for a given strain tensor.
            procedure(dstress_dstrain_interface_2D), deferred :: dstress_dstrain_2D
                !! Computes the tangent modulus (stiffness tensor) for a given strain tensor.
            generic, public :: stress => stress_3D, stress_2D
            generic, public :: dstress_dstrain => dstress_dstrain_3D, dstress_dstrain_2D
    end type Base_elasticity

    interface
        pure function stress_interface_3D(self, strain) result(res)
            !! Interface for the `stress` procedure.
            !! Must be implemented by concrete subtypes of `Base_elasticity`.
            use, intrinsic :: iso_fortran_env
            use tensors_types, only : ten_3D2Osym
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self 
                !! The elasticity model object.
            class(ten_3D2Osym), intent(in) :: strain
                !! Input strain tensor (`ten_3D2Osym`).
            type(ten_3D2Osym) :: res
                !! Output stress tensor (`ten_3D2Osym`).
        end function stress_interface_3D

        pure function stress_interface_2D(self, strain) result(res)
            !! Interface for the `stress` procedure.
            !! Must be implemented by concrete subtypes of `Base_elasticity`.
            use, intrinsic :: iso_fortran_env
            use tensors_types, only : ten_2D2Osym
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self 
                !! The elasticity model object.
            class(ten_2D2Osym), intent(in) :: strain
                !! Input strain tensor (`ten_2D2Osym`).
            type(ten_2D2Osym) :: res
                !! Output stress tensor (`ten_2D2Osym`).
        end function stress_interface_2D

        pure function dstress_dstrain_interface_3D(self, strain) result(res)
            !! Interface for the `dstress_dstrain` procedure.
            !! Must be implemented by concrete subtypes of `Base_elasticity`.
            use, intrinsic :: iso_fortran_env
            use tensors_types, only : ten_3D2Osym, ten_3D4O3sym
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object.
            class(ten_3D2Osym), intent(in) :: strain
                !! Input strain tensor (`ten_3D2Osym`) at which the tangent is evaluated.
            type(ten_3D4O3sym) :: res
                !! Output tangent modulus tensor (`ten_3D4O3sym`).
        end function dstress_dstrain_interface_3D

        pure function dstress_dstrain_interface_2D(self, strain) result(res)
            !! Interface for the `dstress_dstrain` procedure.
            !! Must be implemented by concrete subtypes of `Base_elasticity`.
            use, intrinsic :: iso_fortran_env
            use tensors_types, only : ten_2D2Osym, ten_2D4O3sym
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object.
            class(ten_2D2Osym), intent(in) :: strain
                !! Input strain tensor (`ten_3D2Osym`) at which the tangent is evaluated.
            type(ten_2D4O3sym) :: res
                !! Output tangent modulus tensor (`ten_3D4O3sym`).
        end function dstress_dstrain_interface_2D

    end interface

end module