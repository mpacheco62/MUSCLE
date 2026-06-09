module mod_base_elasticity
    !! Module mod_base_elasticity
    !! ==========================
    !!
    !! Defines the abstract base types for elasticity models, including support for
    !! transient (time-dependent) material behaviour such as viscoelasticity, and for
    !! deviatoric stress and stiffness projections.
    !!
    !! This module provides two abstract derived types:
    !!
    !! - `Base_internal_elasticity_variables`: A minimal abstract container for internal
    !!   state variables (e.g., viscous stress tensors in each Maxwell branch). Concrete
    !!   implementations carry the actual data needed by their respective material model.
    !!
    !! - `Base_elasticity`: The fundamental interface for all elasticity models. It exposes
    !!   four families of generic procedures:
    !!
    !!     - `stress(strain)` / `stress(strain, dstrain, dtime, internal)` — compute the
    !!       full Cauchy stress tensor without modifying any state (both are `pure` functions).
    !!     - `stress_dev(strain)` / `stress_dev(strain, dstrain, dtime, internal)` — compute
    !!       the deviatoric part of the Cauchy stress tensor (both are `pure` functions).
    !!     - `dstress_dstrain(strain)` / `dstress_dstrain(strain, dstrain, dtime, internal)` —
    !!       compute the algorithmic tangent modulus without modifying any state (`pure`).
    !!     - `dstress_dstrain_dev(strain)` / `dstress_dstrain_dev(strain, dstrain, dtime, internal)` —
    !!       compute the deviatoric projection of the algorithmic tangent modulus (`pure`).
    !!     - `stress_and_update(stress, strain, dstrain, dtime, internal)` — compute the
    !!       stress tensor AND advance the internal state variables to the end of the time
    !!       step. This is a `pure subroutine` and must only be called once per converged
    !!       increment.
    !!
    !! Design Rationale
    !! ----------------
    !!
    !! ### Pure vs. updating calls
    !!
    !! In an implicit finite-element solver the Newton–Raphson loop evaluates stress and
    !! tangent repeatedly for trial strains that may not converge. Modifying internal state
    !! variables during those iterations would corrupt the solution history. The split between
    !! `stress` / `dstress_dstrain` (pure, read-only) and `stress_and_update` (updates
    !! `internal`) enforces the correct call pattern at the type level:
    !!
    !! ```
    !! do iter = 1, max_iter                              ! Newton–Raphson
    !!     sigma = material%stress(strain, dstrain, dt, iv)          ! pure — iv unchanged
    !!     Calg  = material%dstress_dstrain(strain, dstrain, dt, iv) ! pure — iv unchanged
    !!     ! ... assemble residual, check convergence ...
    !! end do
    !! ! Step has converged — now advance internal variables
    !! call material%stress_and_update(sigma, strain, dstrain, dt, iv)
    !! ```
    !!
    !! ### Deviatoric variants
    !!
    !! `stress_dev` and `dstress_dstrain_dev` return the deviatoric part of the
    !! corresponding full quantities:
    !!
    !! \[
    !!   \boldsymbol{s} = \mathrm{dev}(\boldsymbol{\sigma})
    !!   \qquad
    !!   \mathbf{C}^{\mathrm{dev}} = \mathbb{P} : \mathbf{C}
    !! \]
    !!
    !! where \( \mathbb{P} \) is the fourth-order deviatoric projector. The default
    !! implementations delegate to `stress` / `dstress_dstrain` and apply the `.dev.`
    !! operator. Subtypes may override them if a more efficient direct evaluation
    !! is available.
    !!
    !! ### Default (elastic) behaviour
    !!
    !! The transient overloads are provided as non-deferred (concrete) procedures on
    !! `Base_elasticity`. Their default implementation simply delegates to the
    !! time-independent `stress` / `dstress_dstrain` procedures, which is correct for
    !! purely elastic models. Viscoelastic (or other history-dependent) subtypes override
    !! these procedures to supply time- and history-aware implementations.
    !!
    !! ### Internal variables
    !!
    !! `Base_internal_elasticity_variables` is intentionally empty at the base level so that
    !! the `stress_and_update` signature is uniform across all material models. A concrete
    !! viscoelastic type would extend it to carry, for example, the per-branch viscous stress
    !! tensors \( \mathbf{q}_k \).
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Abstract Derived Types
    !!
    !! - `Base_internal_elasticity_variables` — Abstract container for internal state variables.
    !!   Extend this type to store model-specific history data (e.g., viscous stresses).
    !!
    !! - `Base_elasticity` — Abstract base type for elasticity models.
    !!
    !!     **Deferred procedures** (must be implemented by every concrete subtype):
    !!
    !!     - `stress_3D(strain)` → `ten_3D2Osym` — full stress from 3-D strain (pure).
    !!     - `stress_2D(strain)` → `ten_2D2Osym` — full stress from 2-D strain (pure).
    !!     - `dstress_dstrain_3D(strain)` → `ten_3D4O3sym` — tangent modulus from 3-D strain (pure).
    !!     - `dstress_dstrain_2D(strain)` → `ten_2D4O3sym` — tangent modulus from 2-D strain (pure).
    !!
    !!     **Non-deferred procedures** (override in history-dependent or specialised subtypes):
    !!
    !!     - `stress_dev_3D(strain)` → `ten_3D2Osym` — deviatoric stress, 3-D (pure).
    !!     - `stress_dev_2D(strain)` → `ten_2D2Osym` — deviatoric stress, 2-D (pure).
    !!     - `stress_3D_transient(strain, dstrain, dtime, internal)` → `ten_3D2Osym` (pure).
    !!     - `stress_2D_transient(strain, dstrain, dtime, internal)` → `ten_2D2Osym` (pure).
    !!     - `stress_dev_3D_transient(strain, dstrain, dtime, internal)` → `ten_3D2Osym` (pure).
    !!     - `stress_dev_2D_transient(strain, dstrain, dtime, internal)` → `ten_2D2Osym` (pure).
    !!     - `dstress_dstrain_3D_transient(strain, dstrain, dtime, internal)` → `ten_3D4O3sym` (pure).
    !!     - `dstress_dstrain_2D_transient(strain, dstrain, dtime, internal)` → `ten_2D4O3sym` (pure).
    !!     - `dstress_dstrain_dev_3D(strain)` → `ten_3D4O3sym` — deviatoric tangent, 3-D (pure).
    !!     - `dstress_dstrain_dev_3D_transient(strain, dstrain, dtime, internal)` → `ten_3D4O3sym` (pure).
    !!     - `stress_and_update_3D(stress, strain, dstrain, dtime, internal)` — pure subroutine;
    !!       updates `internal` in-place after a converged increment.
    !!
    !!     **Generic interfaces** (resolve to the above based on argument signature):
    !!
    !!     - `stress` — full stress; dispatches to `stress_3D`, `stress_2D`,
    !!       `stress_3D_transient`, or `stress_2D_transient`.
    !!     - `stress_dev` — deviatoric stress; dispatches to `stress_dev_3D`,
    !!       `stress_dev_2D`, `stress_dev_3D_transient`, or `stress_dev_2D_transient`.
    !!     - `dstress_dstrain` — full tangent modulus; dispatches to `dstress_dstrain_3D`,
    !!       `dstress_dstrain_2D`, `dstress_dstrain_3D_transient`, or
    !!       `dstress_dstrain_2D_transient`.
    !!     - `dstress_dstrain_dev` — deviatoric tangent; dispatches to
    !!       `dstress_dstrain_dev_3D` or `dstress_dstrain_dev_3D_transient`.
    !!     - `stress_and_update` — dispatches to `stress_and_update_3D`.
    !!
    !! Usage Example
    !! -------------
    !!
    !! ### Purely elastic model (no history)
    !!
    !! ```fortran
    !! module mod_isotropic_elasticity
    !!   use mod_base_elasticity
    !!   use muscle_tensors
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type, extends(Base_elasticity) :: IsotropicElasticity
    !!     real(real64) :: E_mod = 0.0D0   ! Young's modulus
    !!     real(real64) :: nu    = 0.0D0   ! Poisson's ratio
    !!   contains
    !!     procedure :: stress_3D          => iso_stress_3D
    !!     procedure :: stress_2D          => iso_stress_2D
    !!     procedure :: dstress_dstrain_3D => iso_stiffness_3D
    !!     procedure :: dstress_dstrain_2D => iso_stiffness_2D
    !!     ! All remaining procedures (transient, deviatoric, stress_and_update)
    !!     ! are inherited from Base_elasticity and delegate to the above.
    !!   end type IsotropicElasticity
    !! end module mod_isotropic_elasticity
    !! ```
    !!
    !! ### Viscoelastic model (generalised Maxwell)
    !!
    !! ```fortran
    !! module mod_generalized_maxwell
    !!   use mod_base_elasticity
    !!   use muscle_tensors
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type, extends(Base_internal_elasticity_variables) :: MaxwellInternalVars
    !!     type(ten_3D2Osym), allocatable :: q(:)   ! per-branch viscous stress tensors
    !!   end type MaxwellInternalVars
    !!
    !!   type, extends(Base_elasticity) :: GeneralizedMaxwell3D
    !!     real(real64) :: E_inf = 0.0D0
    !!     real(real64) :: nu    = 0.0D0
    !!     real(real64), allocatable :: alpha(:)
    !!     real(real64), allocatable :: tau(:)
    !!   contains
    !!     procedure :: stress_3D                    => gm_stress_3D
    !!     procedure :: stress_2D                    => gm_stress_2D
    !!     procedure :: dstress_dstrain_3D            => gm_dstress_dstrain_3D
    !!     procedure :: dstress_dstrain_2D            => gm_dstress_dstrain_2D
    !!     procedure :: stress_3D_transient           => gm_stress_3D_transient
    !!     procedure :: dstress_dstrain_3D_transient   => gm_tangent_3D_transient
    !!     procedure :: stress_and_update_3D          => gm_stress_and_update_3D
    !!   end type GeneralizedMaxwell3D
    !! end module mod_generalized_maxwell
    !! ```
    !!
    !! For tensor type definitions see [[muscle_tensors]].

    use, intrinsic :: iso_fortran_env
    implicit none
    PRIVATE

    PUBLIC :: Base_internal_elasticity_variables
    type, abstract :: Base_internal_elasticity_variables
        !! Abstract Base Type for Internal Elasticity State Variables
        !! ===========================================================
        !!
        !! Serves as a uniform container for the internal (history) variables of any
        !! elasticity model. For purely elastic models no concrete extension is needed.
        !! History-dependent models (e.g., viscoelastic, viscoplastic) must extend this
        !! type to store their specific state data, such as per-branch viscous stress
        !! tensors \( \mathbf{q}_k \) in a generalised Maxwell model.
        !!
        !! An instance of a concrete subtype is passed to `stress_and_update` and to the
        !! transient overloads of `stress` and `dstress_dstrain`, ensuring a single,
        !! polymorphic interface across all material models.
    end type Base_internal_elasticity_variables


    PUBLIC :: Base_elasticity
    type, abstract :: Base_elasticity
        !! Abstract Base Type for Elasticity Models
        !! ==========================================
        !!
        !! Defines the fundamental interface that every concrete elasticity model must satisfy.
        !! The type exposes five families of generic procedures:
        !!
        !! - `stress` — full Cauchy stress tensor (pure, no side effects).
        !! - `stress_dev` — deviatoric part of the Cauchy stress tensor (pure).
        !! - `dstress_dstrain` — algorithmic tangent modulus (pure, no side effects).
        !! - `dstress_dstrain_dev` — deviatoric projection of the tangent modulus (pure).
        !! - `stress_and_update` — stress tensor plus in-place advance of internal state
        !!   variables at the end of a converged increment (pure subroutine).
        !!
        !! Each generic resolves to a 2-D or 3-D variant, and optionally to a transient
        !! variant that accepts the strain increment `dstrain`, a time-step size `dtime`,
        !! and an internal-variable container `internal`.
        !!
        !! **Subtyping contract**
        !!
        !! A concrete subtype MUST provide implementations for the four deferred procedures:
        !! `stress_3D`, `stress_2D`, `dstress_dstrain_3D`, and `dstress_dstrain_2D`.
        !!
        !! All remaining procedures have concrete default implementations in `Base_elasticity`
        !! and need only be overridden when specialised behaviour is required (e.g., a
        !! history-dependent transient response, or a more efficient deviatoric evaluation).
        contains
            ! --- Deferred: must be implemented by every concrete subtype --------
            procedure(stress_interface_3D), deferred :: stress_3D
                !! Pure function. Computes the full 3-D Cauchy stress tensor
                !! \( \boldsymbol{\sigma} \) from a given symmetric strain tensor
                !! (`ten_3D2Osym`). Must be implemented by every concrete subtype.
            procedure(stress_interface_2D), deferred :: stress_2D
                !! Pure function. Computes the full 2-D Cauchy stress tensor
                !! \( \boldsymbol{\sigma} \) from a given symmetric strain tensor
                !! (`ten_2D2Osym`). Must be implemented by every concrete subtype.
            procedure(dstress_dstrain_interface_3D), deferred :: dstress_dstrain_3D
                !! Pure function. Computes the 3-D algorithmic tangent modulus
                !! \( \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \)
                !! (`ten_3D4O3sym`) at a given strain state.
                !! Must be implemented by every concrete subtype.
            procedure(dstress_dstrain_interface_2D), deferred :: dstress_dstrain_2D
                !! Pure function. Computes the 2-D algorithmic tangent modulus
                !! \( \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \)
                !! (`ten_2D4O3sym`) at a given strain state.
                !! Must be implemented by every concrete subtype.

            ! --- Non-deferred: deviatoric (time-independent) --------------------
            procedure :: stress_dev_3D
                !! Pure function. Returns the deviatoric part of the 3-D Cauchy stress:
                !! \( \boldsymbol{s} = \mathrm{dev}(\boldsymbol{\sigma}) \).
                !! Default implementation applies the `.dev.` operator to `stress_3D`.
                !! Override for a more efficient direct evaluation if available.
            procedure :: stress_dev_2D
                !! Pure function. Returns the deviatoric part of the 2-D Cauchy stress.
                !! Default implementation applies the `.dev.` operator to `stress_2D`.
                !! Override for a more efficient direct evaluation if available.
            procedure :: dstress_dstrain_dev_3D
                !! Pure function. Returns the deviatoric projection of the 3-D tangent:
                !! \( \mathbf{C}^{\mathrm{dev}} = \mathbb{P} : \mathbf{C} \).
                !! Default implementation applies the deviatoric projector to
                !! `dstress_dstrain_3D` using fourth-order identity tensors.
                !! Override for a more efficient direct evaluation if available.

            ! --- Non-deferred: transient full stress/tangent --------------------
            procedure :: stress_3D_transient
                !! Pure function. Transient 3-D stress evaluation.
                !! Accepts the strain increment `dstrain`, time-step size `dtime`, and
                !! a read-only internal-variable container `internal`. Default delegates
                !! to `stress_3D`, ignoring `dstrain`, `dtime`, and `internal`.
                !! Override in history-dependent subtypes.
            procedure :: stress_2D_transient
                !! Pure function. Transient 2-D stress evaluation.
                !! Default delegates to `stress_2D`. Override in history-dependent subtypes.
            procedure :: dstress_dstrain_3D_transient
                !! Pure function. Transient 3-D algorithmic tangent modulus.
                !! Default delegates to `dstress_dstrain_3D`. Override in history-dependent
                !! subtypes. Failure to override when the stress is history-dependent will
                !! cause loss of quadratic convergence in Newton–Raphson iterations.
            procedure :: dstress_dstrain_2D_transient
                !! Pure function. Transient 2-D algorithmic tangent modulus.
                !! Default delegates to `dstress_dstrain_2D`.
                !! Override in history-dependent subtypes.

            ! --- Non-deferred: transient deviatoric stress/tangent --------------
            procedure :: stress_dev_3D_transient
                !! Pure function. Transient deviatoric 3-D stress evaluation.
                !! Default applies `.dev.` to the result of `stress_3D_transient`
                !! (via `stress_3D` in the elastic default).
                !! Override in history-dependent subtypes if needed.
            procedure :: stress_dev_2D_transient
                !! Pure function. Transient deviatoric 2-D stress evaluation.
                !! Default applies `.dev.` to the result of `stress_2D_transient`
                !! (via `stress_2D` in the elastic default).
                !! Override in history-dependent subtypes if needed.
            procedure :: dstress_dstrain_dev_3D_transient
                !! Pure function. Transient deviatoric projection of the 3-D tangent.
                !! Default applies the deviatoric projector to `dstress_dstrain_3D`
                !! (time-independent fallback). Override in history-dependent subtypes
                !! to use the consistent algorithmic tangent.

            ! --- Non-deferred: state update -------------------------------------
            procedure :: stress_and_update_3D
                !! Pure subroutine. Computes the 3-D stress tensor AND advances the
                !! internal state variables `internal` to the end of the current
                !! converged increment.
                !!
                !! Must be called ONCE per converged time step, AFTER the Newton–Raphson
                !! loop has completed. Default delegates to `stress_3D` and leaves
                !! `internal` unchanged (correct for purely elastic models).
                !! Override in history-dependent subtypes to advance internal variables.

            ! --- Generic interfaces ---------------------------------------------
            generic, public :: stress => stress_3D, stress_2D, &
                                          stress_3D_transient, stress_2D_transient
                !! Generic interface for full stress evaluation (pure).
                !! Dispatches based on the type of `strain` and presence of
                !! `dstrain` / `dtime` / `internal`:
                !!   - `stress(strain)`                         → `stress_3D` or `stress_2D`
                !!   - `stress(strain, dstrain, dtime, internal)` → `stress_3D_transient`
                !!                                                 or `stress_2D_transient`
            generic, public :: stress_dev => stress_dev_3D, stress_dev_2D, &
                                              stress_dev_3D_transient, stress_dev_2D_transient
                !! Generic interface for deviatoric stress evaluation (pure).
                !! Dispatches analogously to `stress` but returns
                !! \( \boldsymbol{s} = \mathrm{dev}(\boldsymbol{\sigma}) \).
            generic, public :: dstress_dstrain => dstress_dstrain_3D, dstress_dstrain_2D, &
                                                    dstress_dstrain_3D_transient, &
                                                    dstress_dstrain_2D_transient
                !! Generic interface for the full algorithmic tangent modulus (pure).
                !! Dispatches analogously to `stress` based on dimensionality and the
                !! presence of `dstrain` / `dtime` / `internal`.
            generic, public :: dstress_dstrain_dev => dstress_dstrain_dev_3D, &
                                                       dstress_dstrain_dev_3D_transient
                !! Generic interface for the deviatoric tangent modulus (pure).
                !! Dispatches to `dstress_dstrain_dev_3D` (time-independent) or
                !! `dstress_dstrain_dev_3D_transient` (transient) based on the presence
                !! of `dstrain` / `dtime` / `internal`.
                !! Note: no 2-D deviatoric tangent variant is currently provided.
            generic, public :: stress_and_update => stress_and_update_3D
                !! Generic interface for the stress-and-state-update subroutine.
                !! Currently resolves to `stress_and_update_3D`.
                !! Additional 2-D variants may be added in future.

    end type Base_elasticity

    interface
        pure function stress_interface_3D(self, strain) result(res)
            !! Abstract interface for `stress_3D`.
            !!
            !! Computes the second-order symmetric Cauchy stress tensor (`ten_3D2Osym`)
            !! from a given second-order symmetric strain tensor (`ten_3D2Osym`).
            !! This interface must be satisfied by every concrete subtype of `Base_elasticity`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_3D2Osym
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_3D2Osym), intent(in) :: strain
                !! Input second-order symmetric strain tensor (`ten_3D2Osym`).
            type(ten_3D2Osym) :: res
                !! Output second-order symmetric Cauchy stress tensor (`ten_3D2Osym`).
        end function stress_interface_3D

        pure function stress_interface_2D(self, strain) result(res)
            !! Abstract interface for `stress_2D`.
            !!
            !! Computes the second-order symmetric Cauchy stress tensor (`ten_2D2Osym`)
            !! from a given second-order symmetric strain tensor (`ten_2D2Osym`).
            !! This interface must be satisfied by every concrete subtype of `Base_elasticity`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_2D2Osym
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_2D2Osym), intent(in) :: strain
                !! Input second-order symmetric strain tensor (`ten_2D2Osym`).
            type(ten_2D2Osym) :: res
                !! Output second-order symmetric Cauchy stress tensor (`ten_2D2Osym`).
        end function stress_interface_2D

        pure function dstress_dstrain_interface_3D(self, strain) result(res)
            !! Abstract interface for `dstress_dstrain_3D`.
            !!
            !! Computes the fourth-order major-minor-symmetric tangent modulus tensor
            !! (`ten_3D4O3sym`) representing \( \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \)
            !! evaluated at the given strain state.
            !! This interface must be satisfied by every concrete subtype of `Base_elasticity`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_3D2Osym, ten_3D4O3sym
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_3D2Osym), intent(in) :: strain
                !! Input second-order symmetric strain tensor (`ten_3D2Osym`) at which
                !! the tangent is evaluated.
            type(ten_3D4O3sym) :: res
                !! Output fourth-order tangent modulus tensor (`ten_3D4O3sym`).
        end function dstress_dstrain_interface_3D

        pure function dstress_dstrain_interface_2D(self, strain) result(res)
            !! Abstract interface for `dstress_dstrain_2D`.
            !!
            !! Computes the fourth-order major-minor-symmetric tangent modulus tensor
            !! (`ten_2D4O3sym`) representing \( \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \)
            !! evaluated at the given strain state.
            !! This interface must be satisfied by every concrete subtype of `Base_elasticity`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_2D2Osym, ten_2D4O3sym
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_2D2Osym), intent(in) :: strain
                !! Input second-order symmetric strain tensor (`ten_2D2Osym`) at which
                !! the tangent is evaluated.
            type(ten_2D4O3sym) :: res
                !! Output fourth-order tangent modulus tensor (`ten_2D4O3sym`).
        end function dstress_dstrain_interface_2D

    end interface

    contains

        ! =====================================================================
        ! Deviatoric stress — time-independent
        ! =====================================================================

        pure function stress_dev_3D(self, strain) result(res)
            !! Default implementation of the deviatoric 3-D stress.
            !!
            !! Returns the deviatoric part of the Cauchy stress tensor:
            !! \( \boldsymbol{s} = \mathrm{dev}(\boldsymbol{\sigma}) =
            !! \boldsymbol{\sigma} - \tfrac{1}{3}\,\mathrm{tr}(\boldsymbol{\sigma})\,\mathbf{I} \)
            !!
            !! Delegates to `stress_3D` and applies the `.dev.` operator.
            !! Override in a subtype if a more efficient direct evaluation is available.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_3D2Osym, operator(.dev.)
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_3D2Osym),     intent(in) :: strain
                !! Input second-order symmetric strain tensor (`ten_3D2Osym`).
            type(ten_3D2Osym) :: res
                !! Output deviatoric stress tensor (`ten_3D2Osym`).
            res = .dev. self%stress(strain)
        end function stress_dev_3D


        pure function stress_dev_2D(self, strain) result(res)
            !! Default implementation of the deviatoric 2-D stress.
            !!
            !! Returns the deviatoric part of the 2-D Cauchy stress tensor.
            !! Delegates to `stress_2D` and applies the `.dev.` operator.
            !! Override in a subtype if a more efficient direct evaluation is available.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_2D2Osym, operator(.dev.)
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_2D2Osym),     intent(in) :: strain
                !! Input second-order symmetric strain tensor (`ten_2D2Osym`).
            type(ten_2D2Osym) :: res
                !! Output deviatoric stress tensor (`ten_2D2Osym`).
            res = .dev. self%stress(strain)
        end function stress_dev_2D


        ! =====================================================================
        ! Deviatoric tangent — time-independent
        ! =====================================================================

        pure function dstress_dstrain_dev_3D(self, strain) result(res)
            !! Default implementation of the deviatoric 3-D tangent modulus.
            !!
            !! Returns the deviatoric projection of the full tangent modulus:
            !! \( \mathbf{C}^{\mathrm{dev}} = \mathbb{P} : \mathbf{C} \)
            !!
            !! where \( \mathbb{P} = \mathbb{I}^S - \tfrac{1}{3}\,\mathbf{I} \otimes \mathbf{I} \)
            !! is the fourth-order deviatoric projector, assembled from the symmetric
            !! fourth-order identity `iden_4O4T` and the trace projector `iden_4O3T`.
            !!
            !! The result is converted back to a minor-symmetric `ten_3D4O3sym` via
            !! `convert_3sym()`. Override in a subtype for a more efficient evaluation.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors
            class(Base_elasticity), intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_3D2Osym),     intent(in) :: strain
                !! Input second-order symmetric strain tensor (`ten_3D2Osym`) at which
                !! the tangent is evaluated.
            type(ten_3D4O3sym) :: res
                !! Output deviatoric fourth-order tangent modulus tensor (`ten_3D4O3sym`).
            type(ten_3D4O2sym) :: temp
            type(iden_4O3T) :: i3
            type(iden_4O4T) :: i4

            temp = (i4 - (1.0D0/3.0D0) * i3) .ddot. self%dstress_dstrain(strain)
            res  = temp%convert_3sym()
        end function dstress_dstrain_dev_3D


        ! =====================================================================
        ! Full stress — transient
        ! =====================================================================

        pure function stress_3D_transient(self, strain, dstrain, dtime, internal) result(res)
            !! Default implementation of the transient 3-D stress function.
            !!
            !! Returns the Cauchy stress tensor for the given strain state. This default
            !! delegates to `stress_3D`, ignoring `dstrain`, `dtime`, and `internal`,
            !! which is correct for purely elastic models.
            !!
            !! History-dependent subtypes (e.g., viscoelastic) must override this procedure
            !! to incorporate the viscous contributions stored in `internal`.
            !!
            !! Note: this procedure is `pure` — it does NOT modify `internal`.
            !! To advance the internal state, use `stress_and_update` instead.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_3D2Osym
            class(Base_elasticity),                    intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_3D2Osym),                        intent(in) :: strain
                !! Current strain tensor \( \boldsymbol{\varepsilon}^{n+1} \) (`ten_3D2Osym`).
            class(ten_3D2Osym),                        intent(in) :: dstrain
                !! Strain increment \( \Delta\boldsymbol{\varepsilon} \) (`ten_3D2Osym`).
                !! Ignored by this default; used by history-dependent overrides.
            real(real64),                              intent(in) :: dtime
                !! Time-step size \( \Delta t \). Ignored by this default.
            class(Base_internal_elasticity_variables), intent(in) :: internal
                !! Internal state variables (read-only). Ignored by this default.
            type(ten_3D2Osym) :: res
                !! Output Cauchy stress tensor (`ten_3D2Osym`).
            res = self%stress(strain)
        end function stress_3D_transient


        pure function stress_2D_transient(self, strain, dstrain, dtime, internal) result(res)
            !! Default implementation of the transient 2-D stress function.
            !!
            !! Returns the 2-D Cauchy stress tensor. This default delegates to `stress_2D`,
            !! ignoring `dstrain`, `dtime`, and `internal`.
            !! History-dependent subtypes must override this procedure.
            !!
            !! Note: this procedure is `pure` — it does NOT modify `internal`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_2D2Osym
            class(Base_elasticity),                    intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_2D2Osym),                        intent(in) :: strain
                !! Current strain tensor \( \boldsymbol{\varepsilon}^{n+1} \) (`ten_2D2Osym`).
            class(ten_2D2Osym),                        intent(in) :: dstrain
                !! Strain increment \( \Delta\boldsymbol{\varepsilon} \) (`ten_2D2Osym`).
                !! Ignored by this default; used by history-dependent overrides.
            real(real64),                              intent(in) :: dtime
                !! Time-step size \( \Delta t \). Ignored by this default.
            class(Base_internal_elasticity_variables), intent(in) :: internal
                !! Internal state variables (read-only). Ignored by this default.
            type(ten_2D2Osym) :: res
                !! Output Cauchy stress tensor (`ten_2D2Osym`).
            res = self%stress(strain)
        end function stress_2D_transient


        ! =====================================================================
        ! Deviatoric stress — transient
        ! =====================================================================

        pure function stress_dev_3D_transient(self, strain, dstrain, dtime, internal) result(res)
            !! Default implementation of the transient deviatoric 3-D stress function.
            !!
            !! Returns \( \boldsymbol{s} = \mathrm{dev}(\boldsymbol{\sigma}) \) for the
            !! given strain state. This default applies `.dev.` to the result of
            !! `stress_3D` (i.e., the time-independent elastic stress), ignoring
            !! `dstrain`, `dtime`, and `internal`.
            !!
            !! Override in history-dependent subtypes to return the deviatoric part of
            !! the full transient stress (e.g., by applying `.dev.` to the viscoelastic
            !! stress or by overriding `stress_3D_transient` and delegating here).
            !!
            !! Note: this procedure is `pure` — it does NOT modify `internal`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_3D2Osym, operator(.dev.)
            class(Base_elasticity),                    intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_3D2Osym),                        intent(in) :: strain
                !! Current strain tensor \( \boldsymbol{\varepsilon}^{n+1} \) (`ten_3D2Osym`).
            class(ten_3D2Osym),                        intent(in) :: dstrain
                !! Strain increment \( \Delta\boldsymbol{\varepsilon} \) (`ten_3D2Osym`).
                !! Ignored by this default; used by history-dependent overrides.
            real(real64),                              intent(in) :: dtime
                !! Time-step size \( \Delta t \). Ignored by this default.
            class(Base_internal_elasticity_variables), intent(in) :: internal
                !! Internal state variables (read-only). Ignored by this default.
            type(ten_3D2Osym) :: res
                !! Output deviatoric stress tensor (`ten_3D2Osym`).
            res = .dev. self%stress(strain)
        end function stress_dev_3D_transient


        pure function stress_dev_2D_transient(self, strain, dstrain, dtime, internal) result(res)
            !! Default implementation of the transient deviatoric 2-D stress function.
            !!
            !! Returns the deviatoric part of the 2-D Cauchy stress. This default applies
            !! `.dev.` to `stress_2D`, ignoring `dstrain`, `dtime`, and `internal`.
            !! Override in history-dependent subtypes.
            !!
            !! Note: this procedure is `pure` — it does NOT modify `internal`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_2D2Osym, operator(.dev.)
            class(Base_elasticity),                    intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_2D2Osym),                        intent(in) :: strain
                !! Current strain tensor \( \boldsymbol{\varepsilon}^{n+1} \) (`ten_2D2Osym`).
            class(ten_2D2Osym),                        intent(in) :: dstrain
                !! Strain increment \( \Delta\boldsymbol{\varepsilon} \) (`ten_2D2Osym`).
                !! Ignored by this default; used by history-dependent overrides.
            real(real64),                              intent(in) :: dtime
                !! Time-step size \( \Delta t \). Ignored by this default.
            class(Base_internal_elasticity_variables), intent(in) :: internal
                !! Internal state variables (read-only). Ignored by this default.
            type(ten_2D2Osym) :: res
                !! Output deviatoric stress tensor (`ten_2D2Osym`).
            res = .dev. self%stress(strain)
        end function stress_dev_2D_transient


        ! =====================================================================
        ! Full tangent — transient
        ! =====================================================================

        pure function dstress_dstrain_3D_transient(self, strain, dstrain, dtime, internal) result(res)
            !! Default implementation of the transient 3-D algorithmic tangent modulus.
            !!
            !! Returns \( \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} \)
            !! consistent with the time-integration scheme. This default delegates to
            !! `dstress_dstrain_3D`, ignoring `dstrain`, `dtime`, and `internal`, which
            !! is correct for purely elastic models where the tangent is time-independent.
            !!
            !! History-dependent subtypes MUST override this procedure to return the
            !! consistent algorithmic tangent matching the recurrence formula used in
            !! `stress_and_update_3D`. Failure to do so causes loss of quadratic
            !! convergence in Newton–Raphson iterations.
            !!
            !! Note: this procedure is `pure` — it does NOT modify `internal`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_3D2Osym, ten_3D4O3sym
            class(Base_elasticity),                    intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_3D2Osym),                        intent(in) :: strain
                !! Current strain tensor (`ten_3D2Osym`) at which the tangent is evaluated.
            class(ten_3D2Osym),                        intent(in) :: dstrain
                !! Strain increment \( \Delta\boldsymbol{\varepsilon} \) (`ten_3D2Osym`).
                !! Ignored by this default; used by history-dependent overrides.
            real(real64),                              intent(in) :: dtime
                !! Time-step size \( \Delta t \). Ignored by this default.
            class(Base_internal_elasticity_variables), intent(in) :: internal
                !! Internal state variables (read-only). Ignored by this default.
            type(ten_3D4O3sym) :: res
                !! Output fourth-order algorithmic tangent modulus tensor (`ten_3D4O3sym`).
            res = self%dstress_dstrain(strain)
        end function dstress_dstrain_3D_transient


        pure function dstress_dstrain_2D_transient(self, strain, dstrain, dtime, internal) result(res)
            !! Default implementation of the transient 2-D algorithmic tangent modulus.
            !!
            !! Returns the 2-D tangent modulus. This default delegates to
            !! `dstress_dstrain_2D`, ignoring `dstrain`, `dtime`, and `internal`.
            !! History-dependent subtypes must override this procedure.
            !!
            !! Note: this procedure is `pure` — it does NOT modify `internal`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_2D2Osym, ten_2D4O3sym
            class(Base_elasticity),                    intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_2D2Osym),                        intent(in) :: strain
                !! Current strain tensor (`ten_2D2Osym`) at which the tangent is evaluated.
            class(ten_2D2Osym),                        intent(in) :: dstrain
                !! Strain increment \( \Delta\boldsymbol{\varepsilon} \) (`ten_2D2Osym`).
                !! Ignored by this default; used by history-dependent overrides.
            real(real64),                              intent(in) :: dtime
                !! Time-step size \( \Delta t \). Ignored by this default.
            class(Base_internal_elasticity_variables), intent(in) :: internal
                !! Internal state variables (read-only). Ignored by this default.
            type(ten_2D4O3sym) :: res
                !! Output fourth-order algorithmic tangent modulus tensor (`ten_2D4O3sym`).
            res = self%dstress_dstrain(strain)
        end function dstress_dstrain_2D_transient


        ! =====================================================================
        ! Deviatoric tangent — transient
        ! =====================================================================

        pure function dstress_dstrain_dev_3D_transient(self, strain, dstrain, dtime, internal) result(res)
            !! Default implementation of the transient deviatoric 3-D tangent modulus.
            !!
            !! Returns the deviatoric projection of the tangent modulus (same formula
            !! as `dstress_dstrain_dev_3D`) evaluated at the time-independent elastic
            !! tangent. This default ignores `dstrain`, `dtime`, and `internal`.
            !!
            !! History-dependent subtypes that also override `dstress_dstrain_3D_transient`
            !! should override this procedure as well to ensure the deviatoric tangent is
            !! consistent with the algorithmic tangent of the transient model.
            !!
            !! Note: this procedure is `pure` — it does NOT modify `internal`.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors
            class(Base_elasticity),                    intent(in) :: self
                !! The elasticity model object (read-only).
            class(ten_3D2Osym),                        intent(in) :: strain
                !! Current strain tensor (`ten_3D2Osym`) at which the tangent is evaluated.
            class(ten_3D2Osym),                        intent(in) :: dstrain
                !! Strain increment \( \Delta\boldsymbol{\varepsilon} \) (`ten_3D2Osym`).
                !! Ignored by this default; used by history-dependent overrides.
            real(real64),                              intent(in) :: dtime
                !! Time-step size \( \Delta t \). Ignored by this default.
            class(Base_internal_elasticity_variables), intent(in) :: internal
                !! Internal state variables (read-only). Ignored by this default.
            type(ten_3D4O3sym) :: res
                !! Output deviatoric fourth-order tangent modulus tensor (`ten_3D4O3sym`).
            type(ten_3D4O2sym) :: temp
            type(iden_4O3T) :: i3
            type(iden_4O4T) :: i4

            temp = (i4 - (1.0D0/3.0D0) * i3) .ddot. self%dstress_dstrain(strain)
            res  = temp%convert_3sym()
        end function dstress_dstrain_dev_3D_transient


        ! =====================================================================
        ! State update
        ! =====================================================================

        pure subroutine stress_and_update_3D(self, stress, strain, dstrain, dtime, internal)
            !! Default implementation of the combined 3-D stress-evaluation and
            !! state-update subroutine.
            !!
            !! Computes the Cauchy stress tensor `stress` for the given converged strain
            !! state `strain` and advances the internal state variables `internal` to the
            !! end of the current time increment.
            !!
            !! **When to call:** once per converged increment, AFTER the Newton–Raphson
            !! loop has completed. Calling during iterations corrupts the solution history.
            !!
            !! This default delegates to `stress_3D` and leaves `internal` unchanged,
            !! which is correct for purely elastic models. History-dependent subtypes must
            !! override this procedure to advance their internal variables, e.g.:
            !!
            !! \[
            !!   \mathbf{q}_k^{n+1} = e^{-\Delta t/\tau_k}\,\mathbf{q}_k^{n}
            !!     + \alpha_k \left(1 - e^{-\Delta t/\tau_k}\right)
            !!       \Delta\boldsymbol{\sigma}^\infty
            !! \]
            !!
            !! Note: despite modifying `internal` via `intent(inout)`, this subroutine
            !! is `pure` because all side effects are confined to dummy arguments.
            use, intrinsic :: iso_fortran_env
            use muscle_tensors, only : ten_3D2Osym
            class(Base_elasticity),                    intent(in)    :: self
                !! The elasticity model object (read-only).
            type(ten_3D2Osym),                         intent(out)   :: stress
                !! Output second-order symmetric Cauchy stress tensor (`ten_3D2Osym`).
            class(ten_3D2Osym),                        intent(in)    :: strain
                !! Converged strain tensor \( \boldsymbol{\varepsilon}^{n+1} \) (`ten_3D2Osym`).
            class(ten_3D2Osym),                        intent(in)    :: dstrain
                !! Strain increment \( \Delta\boldsymbol{\varepsilon} \) (`ten_3D2Osym`).
                !! Ignored by this default; used by history-dependent overrides.
            real(real64),                              intent(in)    :: dtime
                !! Time-step size \( \Delta t \) [same units as relaxation times].
            class(Base_internal_elasticity_variables), intent(inout) :: internal
                !! Internal state variables. Updated in-place at the end of the converged
                !! increment. Left unchanged by this default implementation.
            stress = self%stress(strain)
        end subroutine stress_and_update_3D

end module