! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_solver_closest_point
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_base, only : Base_hardening_laws
    use muscle_yield_base, only : Base_yield_critera
    use muscle_elasticity_base, only : Base_elasticity
    use muscle_plastic_history, only : Plastic_material_history, Plastic_state

    private

    integer, parameter, public :: STATUS_UNDEFINED=-1
    integer, parameter, public :: STATUS_ELASTIC_CASE=0
    integer, parameter, public :: STATUS_CONVERGED=1
    integer, parameter, public :: STATUS_NONCONVERGED=2
    integer, parameter, public :: STATUS_ITER_CONVERGED    = 10
    integer, parameter, public :: STATUS_ITER_NONCONVERGED = 11


    type, public :: Closest_point
        class(Base_elasticity), allocatable     :: elasticity
        class(Base_hardening_laws), allocatable :: hardening
        class(Base_yield_critera), allocatable  :: yield
        integer, private                        :: iter_nw = 200
        
        contains
        procedure, public :: init => closest_point_init
        procedure, public :: solve => closest_point_solve
        procedure, public :: iter => closest_point_iter
        procedure, public :: tangent => closest_point_tangent
        procedure, public :: tangent_numerical => closest_point_tangent_numerical
    end type Closest_point

    public :: closest_point2

    contains

        pure subroutine closest_point_init(self, elasticity, hardening, yield, iter_nw)
            implicit none
            class(Closest_point), intent(inout) :: self
            class(Base_elasticity), intent(in) :: elasticity
            class(Base_hardening_laws), intent(in) :: hardening
            class(Base_yield_critera), intent(in) :: yield
            integer, optional, intent(in) :: iter_nw
            self%elasticity = elasticity
            self%hardening = hardening
            self%yield = yield

            self%iter_nw = 200
            if (present(iter_nw)) self%iter_nw = iter_nw
        end subroutine closest_point_init

        pure subroutine closest_point_iter(self, strain, history, dgamma, omega, ddgamma, status)
            !! Performs a single Newton-Raphson iteration step with adaptive under-relaxation omega.
            implicit none
            class(Closest_point), intent(in)              :: self
            type(ten_3D2Osym), intent(in)                 :: strain   !! Input total strain tensor
            type(Plastic_material_history), intent(inout) :: history  !! Material history container
            real(real64), intent(inout)                   :: dgamma   !! Accumulated plastic multiplier increment
            real(real64), intent(inout)                   :: omega    !! Relaxation factor (0 < omega <= 1)
            real(real64), intent(out)                     :: ddgamma  !! Newton-Raphson update step (delta^2 gamma)
            integer, intent(out)                          :: status   !! Iteration status

            real(real64), parameter :: TOL = 1D-5

            ! Local variables for the current iteration
            type(ten_3D2Osym)  :: strain_p_n, strain_p
            real(real64)       :: strain_pf_n, strain_pf
            type(ten_3D2Osym)  :: df, residual1, dstrain_p, stress
            type(ten_3D4O3sym) :: elas_tan, hess, ddf
            real(real64)       :: hard, dhard, f, norm_res

            ! 1. Read state t_n
            strain_p_n  = history%state_n%strain_p
            strain_pf_n = history%state_n%strain_pf

            ! 2. Read state t_n+1 (iter k)
            strain_p  = history%state_np1%strain_p
            strain_pf = history%state_np1%strain_pf

            ! 3. Elastic Tangent & Current Stress
            elas_tan = self%elasticity%dstress_dstrain(strain - strain_p)
            stress   = self%elasticity%stress(strain - strain_p)
            hard     = self%hardening%stress(strain_pf)
            f        = self%yield%stress_eq(stress) - hard

            ! 4. Derivatives & Plastic Strain Residual
            dhard     = self%hardening%dstress_dep(strain_pf)
            df        = self%yield%dstressEq_dstress(stress)
            residual1 = strain_p_n - strain_p + (dgamma * df)
            norm_res  = sqrt(sum(residual1%vals**2))

            ! 5. Check Convergence
            if (abs(f) < TOL .and. norm_res < TOL) then
                status = STATUS_ITER_CONVERGED
                history%state_np1%stress = stress
                return
            end if

            ! 6. Compute Hessian & Raw Newton-Raphson correction (ddgamma)
            ddf  = self%yield%ddstressEq_ddstress(stress)
            hess = .inv. ((.inv. elas_tan) + dgamma * ddf)

            ddgamma = (f - (df .ddot. hess .ddot. residual1)) / ((df .ddot. hess .ddot. df) + dhard)

            ! 7. Apply Under-Relaxed Step: dgamma = dgamma + omega * ddgamma
            dgamma = dgamma + omega * ddgamma

            ! Enforce non-negativity constraint on plastic multiplier
            if (dgamma < 0.0D0) then
                dgamma = 0.0D0
                omega  = omega * 0.75D0  ! Reduce relaxation factor if hitting non-negative boundary
            end if

            strain_pf = strain_pf_n + dgamma

            dstrain_p = ((.inv. elas_tan) .ddot. hess) .ddot. (residual1 + ddgamma * df)
            strain_p  = strain_p + dstrain_p

            ! 8. Update candidate state t_n+1 (iter k+1)
            history%state_np1%strain_p  = strain_p
            history%state_np1%strain_pf = strain_pf
            history%state_np1%stress    = stress
            status = STATUS_ITER_NONCONVERGED

        end subroutine closest_point_iter

        pure subroutine closest_point_solve(self, strain, history, status, iters)
            implicit none
            class(Closest_point), intent(in)              :: self
            type(ten_3D2Osym), intent(in)                 :: strain
            type(Plastic_material_history), intent(inout) :: history
            integer, intent(out), optional                :: status
            integer, intent(out), optional                :: iters

            real(real64), parameter :: TOL2 = 1D-8
            type(ten_3D2Osym)       :: stress_trial
            real(real64)            :: hard_n, f_trial
            real(real64)            :: dgamma, ddgamma, omega
            integer                 :: i, iter_status, local_status


           ! *** Step 1: Elastic Trial Check ***
            stress_trial = self%elasticity%stress(strain - history%state_n%strain_p)
            hard_n       = self%hardening%stress(history%state_n%strain_pf)
            f_trial      = self%yield%stress_eq(stress_trial) - hard_n

            if (f_trial / max(hard_n, 1.0D-10) <= TOL2) then  ! Elastic Case
                history%state_np1%stress    = stress_trial
                history%state_np1%strain_p  = history%state_n%strain_p
                history%state_np1%strain_pf = history%state_n%strain_pf

                if (present(status)) status = STATUS_ELASTIC_CASE
                if (present(iters))  iters  = 0
                return
            end if

            ! *** Step 2: Initialize Candidate State and Relaxation Factor ***
            history%state_np1 = history%state_n
            dgamma       = 0.0D0
            ddgamma      = 0.0D0
            omega        = 1.0D0  ! Initial unrelaxed factor
            local_status = STATUS_NONCONVERGED

            ! *** Step 3: Newton-Raphson Loop with Adaptive Relaxation ***
            do i = 1, self%iter_nw
                ! Stagnation protection: reduce omega if taking too many iterations
                if (i == 10 .or. i == 20 .or. i == 30 .or. i == 40) then
                    omega = omega * 0.75D0
                end if

                call self%iter(strain, history, dgamma, omega, ddgamma, iter_status)

                if (iter_status == STATUS_ITER_CONVERGED) then
                    local_status = STATUS_CONVERGED
                    exit
                end if
            end do

            if (present(status)) status = local_status
            if (present(iters))  iters  = i
        end subroutine closest_point_solve



        pure subroutine closest_point_tangent(self, strain, history, tangent)
            !! Computes the consistent algorithmic tangent stiffness tensor C_mat = dS/dE.
            !! Pure read-only function that evaluates the tangent at the converged state t_n+1.
            implicit none
            class(Closest_point), intent(in)           :: self
            type(ten_3D2Osym), intent(in)              :: strain   !! Input total strain tensor
            type(Plastic_material_history), intent(in) :: history  !! Material history container (Read-Only)
            type(ten_3D4O2sym), intent(out)            :: tangent  !! Algorithmic tangent tensor C_mat

            type(ten_3D2Osym)  :: df, stress, strain_p
            type(ten_3D4O3sym) :: elas_tan, hess, ddf
            real(real64)       :: dhard, dgamma, strain_pf

            ! 1. Read state variables at t_n+1 from history
            strain_p  = history%state_np1%strain_p
            strain_pf = history%state_np1%strain_pf
            stress    = history%state_np1%stress

            ! 2. Compute step plastic multiplier increment Delta gamma = bar_eps_p_{n+1} - bar_eps_p_n
            dgamma = history%state_np1%strain_pf - history%state_n%strain_pf

            ! 3. Compute elastic tangent at current elastic strain
            elas_tan = self%elasticity%dstress_dstrain(strain - strain_p)

            ! 4. Check for elastic step (Delta gamma = 0)
            if (dgamma <= 1.0D-12) then
                tangent = elas_tan
                return
            end if

            ! 5. Compute derivatives at converged stress state
            dhard = self%hardening%dstress_dep(strain_pf)
            df    = self%yield%dstressEq_dstress(stress)
            ddf   = self%yield%ddstressEq_ddstress(stress)

            ! 6. Algorithmic Hessian matrix: H_alg = [ C^-1 + dgamma * d2f/dsigma2 ]^-1
            hess = .inv. ((.inv. elas_tan) + dgamma * ddf)

            ! 7. Consistent Elastoplastic Tangent Modulus:
            ! C_ep = H_alg - (H_alg : N x N : H_alg) / (N : H_alg : N + H')
            tangent = hess - (.tdotsym. (df .ddot. hess)) / ((df .ddot. hess .ddot. df) + dhard)

        end subroutine closest_point_tangent

        subroutine closest_point_tangent_numerical(self, strain, history, tangent)
            use muscle_math_derivatives
            implicit none
            class(Closest_point), intent(in)           :: self
            type(ten_3D2Osym), intent(in)              :: strain
            type(Plastic_material_history), intent(in) :: history
            type(ten_3D4O2sym), intent(out)            :: tangent

            tangent = derivative(wrapper, strain)

            contains
            pure function wrapper(x) result(stress_out)
                implicit none
                type(ten_3D2Osym), intent(in) :: x
                type(ten_3D2Osym)              :: stress_out
                type(Plastic_material_history) :: local_history

                local_history = history
                call self%solve(x, local_history)
                stress_out = local_history%state_np1%stress
            end function wrapper
        end subroutine closest_point_tangent_numerical

end module