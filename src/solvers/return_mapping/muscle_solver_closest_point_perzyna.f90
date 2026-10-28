! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_solver_closest_point_perzyna
    !! Closest point return mapping for overstress viscoplasticity (Perzyna and similar).
    !!
    !! Unlike `muscle_solver_closest_point` the plastic multiplier is not fixed by the
    !! consistency condition f = 0 but by the flow law dlambda = dt*g(sigma_eq, ep), so f > 0
    !! remains at the end of the step. The residuals are
    !!
    !!     R_ep     = ep_n - ep + dlambda*df/dsigma
    !!     R_lambda = dlambda - dt*g(sigma_eq, ep_n + dlambda)
    !!
    !! See Theory/closest_point/closest_point_perzyna.tex. The time step `dt` is an argument of
    !! `solve`, `iter` and `tangent` (it is not part of the material history).
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_overstress_base, only : Base_overstress_law
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

    integer, parameter :: MAX_BACKTRACK = 10

    type, public :: Closest_point
        class(Base_elasticity), allocatable      :: elasticity
        class(Base_overstress_law), allocatable  :: viscoplastic
        class(Base_yield_critera), allocatable   :: yield
        integer, private                         :: iter_nw = 200
        real(real64), private                    :: tol = 1D-8

        contains
        procedure, public :: init => closest_point_init
        procedure, public :: solve => closest_point_solve
        procedure, public :: iter => closest_point_iter
        procedure, public :: tangent => closest_point_tangent
        procedure, public :: tangent_numerical => closest_point_tangent_numerical
    end type Closest_point

    contains

        pure subroutine closest_point_init(self, elasticity, viscoplastic, yield, iter_nw, tol)
            implicit none
            class(Closest_point), intent(inout) :: self
            class(Base_elasticity), intent(in) :: elasticity
            class(Base_overstress_law), intent(in) :: viscoplastic
            class(Base_yield_critera), intent(in) :: yield
            integer, optional, intent(in) :: iter_nw
            real(real64), optional, intent(in) :: tol
            ! Reallocate only when the dynamic type changes: gfortran 13/14 do not do it on assignment
            ! Nested ifs: .and. does not short-circuit, and same_type_as on an unallocated component
            ! reads its _vptr, which is garbage in an OpenMP PRIVATE copy (segfault)
            if (allocated(self%elasticity)) then
                if (.not. same_type_as(self%elasticity, elasticity)) deallocate(self%elasticity)
            end if
            self%elasticity = elasticity
            if (allocated(self%viscoplastic)) then
                if (.not. same_type_as(self%viscoplastic, viscoplastic)) deallocate(self%viscoplastic)
            end if
            self%viscoplastic = viscoplastic
            if (allocated(self%yield)) then
                if (.not. same_type_as(self%yield, yield)) deallocate(self%yield)
            end if
            self%yield = yield

            self%iter_nw = 200
            self%tol = 1D-8
            if (present(iter_nw)) self%iter_nw = iter_nw
            if (present(tol)) self%tol = tol
        end subroutine closest_point_init

        pure function merit(self, strain, strain_p, strain_p_n, dgamma, strain_pf_n, dt) result(res)
            !! Squared norm of the residuals (R_ep, R_lambda) at a trial state
            class(Closest_point), intent(in) :: self
            type(ten_3D2Osym), intent(in) :: strain, strain_p, strain_p_n
            real(real64), intent(in) :: dgamma, strain_pf_n, dt
            real(real64) :: res

            type(ten_3D2Osym) :: stress, residual1
            real(real64) :: stress_eq

            stress = self%elasticity%stress(strain - strain_p)
            stress_eq = self%yield%stress_eq(stress)
            residual1 = strain_p_n - strain_p + (dgamma*self%yield%dstressEq_dstress(stress))
            res = sum(residual1%vals**2) &
                  + (dgamma - dt*self%viscoplastic%flow_rate(stress_eq, strain_pf_n + dgamma))**2
        end function merit

        pure subroutine closest_point_iter(self, strain, history, dt, ddgamma, status)
            !! Performs a single Newton-Raphson iteration on (R_ep, R_lambda) with a backtracking
            !! line search. The plastic multiplier increment is
            !! dgamma = history%state_np1%strain_pf - history%state_n%strain_pf.
            implicit none
            class(Closest_point), intent(in)              :: self
            type(ten_3D2Osym), intent(in)                 :: strain   !! Input total strain tensor
            type(Plastic_material_history), intent(inout) :: history  !! Material history container
            real(real64), intent(in)                      :: dt       !! Time step
            real(real64), intent(out)                     :: ddgamma  !! Newton-Raphson update step (delta^2 gamma)
            integer, intent(out)                          :: status   !! Iteration status

            type(ten_3D2Osym)  :: strain_p_n, strain_p, strain_p_new
            real(real64)       :: strain_pf_n, strain_pf, dgamma, dgamma_new
            type(ten_3D2Osym)  :: df, residual1, dstrain_p, stress
            type(ten_3D4O3sym) :: elas_tan, elas_inv, hess, ddf
            real(real64)       :: stress_eq, g, c, gp
            real(real64)       :: res_lambda, norm_res, merit0, alpha
            logical            :: joint  ! The criterion evaluates f, df and ddf together
            integer            :: i

            ! 1. Read state t_n and t_n+1 (iter k)
            strain_p_n  = history%state_n%strain_p
            strain_pf_n = history%state_n%strain_pf
            strain_p    = history%state_np1%strain_p
            strain_pf   = history%state_np1%strain_pf
            dgamma      = strain_pf - strain_pf_n

            ! 2. Elastic Tangent, Current Stress, Equivalent Stress and its Gradient (and the Hessian
            ! too when the criterion evaluates the three together)
            elas_tan = self%elasticity%dstress_dstrain(strain - strain_p)
            stress   = self%elasticity%stress(strain - strain_p)
            joint    = self%yield%joint_evaluation()
            call evaluate_yield(self, stress, stress_eq, df, ddf, joint, .false.)

            ! 3. Flow law and its derivatives
            g  = self%viscoplastic%flow_rate(stress_eq, strain_pf)
            c  = self%viscoplastic%dflow_dsigma(stress_eq, strain_pf)
            gp = self%viscoplastic%dflow_dep(stress_eq, strain_pf)

            ! 4. Residuals
            residual1  = strain_p_n - strain_p + (dgamma * df)
            res_lambda = dgamma - dt*g
            norm_res   = sqrt(sum(residual1%vals**2))

            ! 5. Check Convergence
            if (abs(res_lambda) < self%tol .and. norm_res < self%tol) then
                status = STATUS_ITER_CONVERGED
                ddgamma = 0D0
                history%state_np1%stress = stress
                return
            end if

            ! 6. Hessian & Newton-Raphson correction (finite for c = 0: then ddgamma = -dgamma)
            if (.not. joint) ddf = self%yield%ddstressEq_ddstress(stress)
            elas_inv = .inv. elas_tan
            hess = .inv. (elas_inv + dgamma * ddf)

            ddgamma = -(res_lambda + dt*c*(df .ddot. hess .ddot. residual1)) &
                      /(1D0 - dt*gp + dt*c*(df .ddot. hess .ddot. df))
            dstrain_p = (elas_inv .ddot. hess) .ddot. (residual1 + ddgamma * df)

            ! 7. Backtracking line search on the residual norm (Newton can cycle when the
            ! hardening slope is very steep, e.g. Swift near zero plastic strain). dgamma >= 0.
            merit0 = res_lambda**2 + norm_res**2
            alpha = 1D0
            do i = 1, MAX_BACKTRACK
                dgamma_new   = max(dgamma + alpha*ddgamma, 0D0)
                strain_p_new = strain_p + alpha*dstrain_p
                if (merit(self, strain, strain_p_new, strain_p_n, dgamma_new, strain_pf_n, dt) &
                    < merit0) exit
                alpha = 0.5D0*alpha
            end do
            ddgamma = dgamma_new - dgamma

            ! 8. Update candidate state t_n+1 (iter k+1)
            history%state_np1%strain_p  = strain_p_new
            history%state_np1%strain_pf = strain_pf_n + dgamma_new
            history%state_np1%stress    = stress
            status = STATUS_ITER_NONCONVERGED
        end subroutine closest_point_iter

        pure subroutine closest_point_solve(self, strain, history, dt, status, iters)
            implicit none
            class(Closest_point), intent(in)              :: self
            type(ten_3D2Osym), intent(in)                 :: strain
            type(Plastic_material_history), intent(inout) :: history
            real(real64), intent(in)                      :: dt
            integer, intent(out), optional                :: status
            integer, intent(out), optional                :: iters

            type(ten_3D2Osym) :: stress_trial
            real(real64)      :: ddgamma
            integer           :: i, iter_status, local_status

            ! *** Step 1: Elastic Trial Check (no flow if the flow rate vanishes) ***
            stress_trial = self%elasticity%stress(strain - history%state_n%strain_p)
            if (self%viscoplastic%flow_rate(self%yield%stress_eq(stress_trial), &
                                            history%state_n%strain_pf) .le. 0D0) then  ! Elastic Case
                history%state_np1%stress    = stress_trial
                history%state_np1%strain_p  = history%state_n%strain_p
                history%state_np1%strain_pf = history%state_n%strain_pf

                if (present(status)) status = STATUS_ELASTIC_CASE
                if (present(iters))  iters  = 0
                return
            end if

            ! *** Step 2: Initialize Candidate State ***
            history%state_np1 = history%state_n
            local_status = STATUS_NONCONVERGED

            ! *** Step 3: Newton-Raphson Loop ***
            do i = 1, self%iter_nw
                call self%iter(strain, history, dt, ddgamma, iter_status)
                if (iter_status == STATUS_ITER_CONVERGED) then
                    local_status = STATUS_CONVERGED
                    exit
                end if
            end do

            if (present(status)) status = local_status
            if (present(iters))  iters  = min(i, self%iter_nw)
        end subroutine closest_point_solve

        pure subroutine closest_point_tangent(self, strain, history, dt, tangent)
            !! Consistent algorithmic tangent at the converged state t_n+1.
            implicit none
            class(Closest_point), intent(in)           :: self
            type(ten_3D2Osym), intent(in)              :: strain   !! Input total strain tensor
            type(Plastic_material_history), intent(in) :: history  !! Material history container (Read-Only)
            real(real64), intent(in)                   :: dt       !! Time step
            type(ten_3D4O2sym), intent(out)            :: tangent  !! Algorithmic tangent tensor C_mat

            type(ten_3D2Osym)  :: df, stress, strain_p
            type(ten_3D4O3sym) :: elas_tan, hess, ddf
            real(real64)       :: dgamma, strain_pf, stress_eq, c, gp, denom

            strain_p  = history%state_np1%strain_p
            strain_pf = history%state_np1%strain_pf
            stress    = history%state_np1%stress
            dgamma    = strain_pf - history%state_n%strain_pf

            elas_tan = self%elasticity%dstress_dstrain(strain - strain_p)

            ! Elastic step (Delta gamma = 0)
            if (dgamma <= 1.0D-12) then
                tangent = elas_tan
                return
            end if

            call evaluate_yield(self, stress, stress_eq, df, ddf, self%yield%joint_evaluation(), .true.)
            c  = self%viscoplastic%dflow_dsigma(stress_eq, strain_pf)
            gp = self%viscoplastic%dflow_dep(stress_eq, strain_pf)
            hess = .inv. ((.inv. elas_tan) + dgamma * ddf)

            if (c .le. 0D0) then  ! No overstress sensitivity: no plastic correction
                tangent = hess
                return
            end if

            ! denom = (1 - dt*gp + dt*c*n:Xi:n)/(dt*c),  gp = dg/dep = -c*htilde
            denom = (1D0 - dt*gp)/(dt*c) + (df .ddot. hess .ddot. df)
            tangent = hess - (.tdotsym. (df .ddot. hess)) / denom
        end subroutine closest_point_tangent

        pure subroutine evaluate_yield(self, stress, stress_eq, df, ddf, joint, need_hessian)
            ! Evaluate the criterion through its declared route for iter and tangent.
            ! The separate Hessian is deferred until after iter's convergence check.
            implicit none
            class(Closest_point), intent(in) :: self
            type(ten_3D2Osym), intent(in) :: stress
            real(real64), intent(out) :: stress_eq
            type(ten_3D2Osym), intent(out) :: df
            type(ten_3D4O3sym), intent(out) :: ddf
            logical, intent(in) :: joint, need_hessian

            if (joint) then
                call self%yield%stressEq_and_derivatives(stress, stress_eq, df, ddf)
            else
                stress_eq = self%yield%stress_eq(stress)
                df = self%yield%dstressEq_dstress(stress)
                if (need_hessian) ddf = self%yield%ddstressEq_ddstress(stress)
            end if
        end subroutine evaluate_yield

        subroutine closest_point_tangent_numerical(self, strain, history, dt, tangent)
            use muscle_math_derivatives
            implicit none
            class(Closest_point), intent(in)           :: self
            type(ten_3D2Osym), intent(in)              :: strain
            type(Plastic_material_history), intent(in) :: history
            real(real64), intent(in)                   :: dt
            type(ten_3D4O2sym), intent(out)            :: tangent

            tangent = derivative(wrapper, strain)

            contains
            pure function wrapper(x) result(stress_out)
                implicit none
                type(ten_3D2Osym), intent(in) :: x
                type(ten_3D2Osym)              :: stress_out
                type(Plastic_material_history) :: local_history

                local_history = history
                call self%solve(x, local_history, dt)
                stress_out = local_history%state_np1%stress
            end function wrapper
        end subroutine closest_point_tangent_numerical

end module
