! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_visco_base
    !! Module muscle_visco_base
    !! ========================
    !! Abstract base type for linear viscoelastic models (generalized Maxwell, Kelvin-Voigt,
    !! Burgers, ...).
    !!
    !! A linear viscoelastic model integrated over one time step is affine in the elastic strain,
    !! sigma_n+1 = C_alg(dt) : eps_e_n+1 + sigma_hist(state_n, dt). The base splits the work by
    !! how often it has to be done:
    !!
    !! 1. `coefficients(dt)`: everything that depends only on the time step (exponentials,
    !!    algorithmic moduli). Once per step for the whole model.
    !! 2. `freeze(coefs, history)`: the algorithmic elasticity (`Elasticity_frozen`) of one
    !!    material point. Once per step and point; it is then used by the Newton iterations of
    !!    any return mapping as an ordinary elasticity.
    !! 3. `update(coefs, strain_e, history)`: advances the history tensors from state_n to
    !!    state_np1 with the elastic strain at t_n+1. Once per converged local solve.
    !!
    !! ```fortran
    !! coefs = visco%coefficients(dt)                     ! once per step (all points)
    !! elas  = visco%freeze(coefs, vhist)                 ! per point
    !! call solver%init(elasticity=elas, ...)             ! return mapping unchanged
    !! call solver%solve(strain, phist, ...)
    !! call visco%update(coefs, strain - phist%state_np1%strain_p, vhist)
    !! ! global convergence: call vhist%commit()
    !! ```
    !!
    !! Models whose step map is not affine (Schapery, stress-dependent viscosity) need a local
    !! Newton and do not fit this base.
    !!
    !! Historical OO API (host only, removed in Phase 7, D-017.1): the models are written as
    !! `pure` kernel procedures with an explicit state (`visco_step_coefficients`, `visco_freeze`,
    !! `visco_update` of muscle_visco_laws, D-017.13) and their OO types extending this base are
    !! thin adapters. `Visco_coefficients` is the kernel type `visco_coefs_dp`.
    use muscle_kinds, only : dp
    use muscle_tensors
    use muscle_visco_types, only : Visco_coefficients => visco_coefs_dp
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    implicit none
    private

    public :: Base_viscoelasticity
    public :: Visco_coefficients

    type, abstract :: Base_viscoelasticity
    contains
        procedure(coefficients_interface), deferred :: coefficients
        procedure(freeze_interface), deferred :: freeze
        procedure(update_interface), deferred :: update
        procedure(n_branches_interface), deferred :: n_branches
        procedure :: init_history
        procedure :: freeze_dt
        procedure :: update_dt
    end type Base_viscoelasticity

    abstract interface
        pure function coefficients_interface(self, dt) result(coefs)
            !! Step-dependent coefficients of the model.
            use muscle_kinds, only : dp
            import Base_viscoelasticity, Visco_coefficients
            class(Base_viscoelasticity), intent(in) :: self
            real(dp), intent(in)                    :: dt
            type(Visco_coefficients)                :: coefs
        end function coefficients_interface

        pure function freeze_interface(self, coefs, history) result(elas)
            !! Algorithmic elasticity over the step, from the converged state history%state_n.
            import Base_viscoelasticity, Visco_coefficients, Viscoelastic_history, Elasticity_frozen
            class(Base_viscoelasticity), intent(in) :: self
            type(Visco_coefficients), intent(in)    :: coefs
            type(Viscoelastic_history), intent(in)  :: history
            type(Elasticity_frozen)                 :: elas
        end function freeze_interface

        pure subroutine update_interface(self, coefs, strain_e, history)
            !! Writes history%state_np1 from history%state_n and the elastic strain at t_n+1.
            use muscle_tensors, only : ten_3D2Osym
            import Base_viscoelasticity, Visco_coefficients, Viscoelastic_history
            class(Base_viscoelasticity), intent(in)   :: self
            type(Visco_coefficients), intent(in)      :: coefs
            type(ten_3D2Osym), intent(in)             :: strain_e
            type(Viscoelastic_history), intent(inout) :: history
        end subroutine update_interface

        pure function n_branches_interface(self) result(res)
            !! Number of history tensors used by the model.
            import Base_viscoelasticity
            class(Base_viscoelasticity), intent(in) :: self
            integer                                 :: res
        end function n_branches_interface
    end interface

contains

    pure subroutine init_history(self, history, strain_e)
        !! Initializes a history with the number of branches of the model (zero history).
        class(Base_viscoelasticity), intent(in)   :: self
        type(Viscoelastic_history), intent(inout) :: history
        type(ten_3D2Osym), intent(in), optional   :: strain_e
        call history%init(self%n_branches(), strain_e)
    end subroutine init_history

    pure function freeze_dt(self, dt, history) result(elas)
        !! Convenience: `freeze` computing the coefficients (use `freeze` in hot loops).
        class(Base_viscoelasticity), intent(in) :: self
        real(dp), intent(in)                    :: dt
        type(Viscoelastic_history), intent(in)  :: history
        type(Elasticity_frozen)                 :: elas
        elas = self%freeze(self%coefficients(dt), history)
    end function freeze_dt

    pure subroutine update_dt(self, dt, strain_e, history)
        !! Convenience: `update` computing the coefficients (use `update` in hot loops).
        class(Base_viscoelasticity), intent(in)   :: self
        real(dp), intent(in)                      :: dt
        type(ten_3D2Osym), intent(in)             :: strain_e
        type(Viscoelastic_history), intent(inout) :: history
        call self%update(self%coefficients(dt), strain_e, history)
    end subroutine update_dt

end module muscle_visco_base
