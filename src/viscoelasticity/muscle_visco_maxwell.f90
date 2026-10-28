! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_visco_maxwell
    !! Module muscle_visco_maxwell
    !! ===========================
    !! Isotropic generalized Maxwell model (Prony series) with deviatoric relaxation:
    !!
    !! \[
    !!   G(t) = G_\infty + \sum_k G_k e^{-t/\tau_k}, \qquad K(t) = K
    !! \]
    !!
    !! The single Maxwell element (G_inf = 0) and the standard linear solid (one branch) are
    !! particular cases. The history tensors are the deviatoric branch stresses s_k, integrated
    !! exactly for an elastic strain that varies linearly within the step:
    !!
    !! \[
    !!   a_k = e^{-\Delta t/\tau_k}, \quad \gamma_k = \frac{1 - a_k}{\Delta t/\tau_k}, \quad
    !!   s_k^{n+1} = a_k s_k^n + 2 G_k \gamma_k (e^{n+1} - e^n)
    !! \]
    !!
    !! which gives G_alg = G_inf + sum(G_k gamma_k) and
    !! sigma_hist = sum(a_k s_k^n - 2 G_k gamma_k e^n). See
    !! Theory/viscoelasticity/generalized_maxwell.tex.
    !!
    !! Coefficients: `decay(k) = a_k`, `weight(k) = 2 G_k gamma_k`.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_visco_state, only : Viscoelastic_history, MAX_VISCO_BRANCHES
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Base_viscoelasticity, Visco_coefficients
    implicit none
    private

    public :: Generalized_maxwell

    type, extends(Base_viscoelasticity) :: Generalized_maxwell
        real(real64), private :: k = 0D0                                !! Bulk modulus
        real(real64), private :: g_inf = 0D0                            !! Long-term shear modulus
        integer, private      :: n = 0                                  !! Number of Maxwell branches
        real(real64), private :: g(MAX_VISCO_BRANCHES) = 0D0            !! Branch shear moduli
        real(real64), private :: tau(MAX_VISCO_BRANCHES) = 1D0          !! Branch relaxation times
    contains
        procedure :: set_parameters
        procedure :: set_parameters_normalized
        procedure :: relaxation_shear
        procedure :: shear_instantaneous
        procedure :: n_branches => n_branches_maxwell
        procedure :: coefficients => coefficients_maxwell
        procedure :: freeze => freeze_maxwell
        procedure :: update => update_maxwell
    end type Generalized_maxwell

contains

    pure subroutine set_parameters(self, bulk, shear_inf, shear, tau)
        !! Sets the Prony series with absolute moduli.
        class(Generalized_maxwell), intent(inout) :: self
        real(real64), intent(in) :: bulk       !! Bulk modulus K
        real(real64), intent(in) :: shear_inf  !! Long-term shear modulus G_inf
        real(real64), intent(in) :: shear(:)   !! Branch shear moduli G_k
        real(real64), intent(in) :: tau(:)     !! Branch relaxation times tau_k

        if (size(shear) /= size(tau)) then
            error stop "ERROR FATAL [Generalized_maxwell%set_parameters]: " // &
                       "'shear' y 'tau' deben tener el mismo tamano."
        end if
        if (size(shear) > MAX_VISCO_BRANCHES) then
            error stop "ERROR FATAL [Generalized_maxwell%set_parameters]: " // &
                       "Numero de ramas mayor que MAX_VISCO_BRANCHES."
        end if
        if (any(tau <= 0D0)) then
            error stop "ERROR FATAL [Generalized_maxwell%set_parameters]: " // &
                       "Los tiempos de relajacion 'tau' deben ser positivos."
        end if

        self%k = bulk
        self%g_inf = shear_inf
        self%n = size(shear)
        self%g = 0D0
        self%tau = 1D0
        self%g(1:self%n) = shear
        self%tau(1:self%n) = tau
    end subroutine set_parameters

    pure subroutine set_parameters_normalized(self, young, poisson, g, tau)
        !! Sets the Prony series as in Abaqus: instantaneous E and nu, and normalized shear
        !! moduli g_k = G_k/G_0 (sum(g) <= 1, sum(g) = 1 is a Maxwell fluid).
        class(Generalized_maxwell), intent(inout) :: self
        real(real64), intent(in) :: young    !! Instantaneous Young modulus
        real(real64), intent(in) :: poisson  !! Instantaneous Poisson ratio
        real(real64), intent(in) :: g(:)     !! Normalized branch shear moduli g_k
        real(real64), intent(in) :: tau(:)   !! Branch relaxation times tau_k
        real(real64) :: g0

        ! sum(g) = 1 (Maxwell fluid) is a valid input, but the sum of a few normalized terms
        ! carries a rounding error below n*epsilon/2: [0.2, 0.4, 0.3, 0.1] gives 1 + 2.2e-16.
        if (sum(g) > 1D0 + 10D0*epsilon(1D0)) then
            error stop "ERROR FATAL [Generalized_maxwell%set_parameters_normalized]: " // &
                       "La suma de 'g' no puede ser mayor que 1."
        end if

        g0 = 0.5D0*young/(1D0 + poisson)
        call self%set_parameters(bulk=young/(3D0*(1D0 - 2D0*poisson)), &
                                 shear_inf=g0*max(0D0, 1D0 - sum(g)), shear=g0*g, tau=tau)
    end subroutine set_parameters_normalized

    pure function relaxation_shear(self, t) result(res)
        !! Shear relaxation modulus G(t).
        class(Generalized_maxwell), intent(in) :: self
        real(real64), intent(in) :: t
        real(real64) :: res
        res = self%g_inf + sum(self%g(1:self%n)*exp(-t/self%tau(1:self%n)))
    end function relaxation_shear

    pure function shear_instantaneous(self) result(res)
        !! Instantaneous shear modulus G_0 = G(0).
        class(Generalized_maxwell), intent(in) :: self
        real(real64) :: res
        res = self%g_inf + sum(self%g(1:self%n))
    end function shear_instantaneous

    pure function n_branches_maxwell(self) result(res)
        class(Generalized_maxwell), intent(in) :: self
        integer :: res
        res = self%n
    end function n_branches_maxwell

    pure function coefficients_maxwell(self, dt) result(coefs)
        class(Generalized_maxwell), intent(in) :: self
        real(real64), intent(in) :: dt
        type(Visco_coefficients) :: coefs
        real(real64) :: x, gamma
        integer :: i

        if (dt < 0D0) then
            error stop "ERROR FATAL [Generalized_maxwell%coefficients]: " // &
                       "El paso de tiempo 'dt' no puede ser negativo (dt = 0 es valido)."
        end if
        coefs%dt = dt
        coefs%k_alg = self%k
        coefs%g_alg = self%g_inf
        do i = 1, self%n
            x = dt/self%tau(i)
            coefs%decay(i) = exp(-x)
            if (x < 1D-5) then  ! (1 - exp(-x))/x without cancellation
                gamma = 1D0 - x/2D0 + x*x/6D0
            else
                gamma = (1D0 - coefs%decay(i))/x
            end if
            coefs%weight(i) = 2D0*self%g(i)*gamma
            coefs%g_alg = coefs%g_alg + self%g(i)*gamma
        end do
    end function coefficients_maxwell

    pure function freeze_maxwell(self, coefs, history) result(elas)
        class(Generalized_maxwell), intent(in) :: self
        type(Visco_coefficients), intent(in)   :: coefs
        type(Viscoelastic_history), intent(in) :: history
        type(Elasticity_frozen)                :: elas
        type(ten_3D2Osym) :: e_n, hist
        integer :: i

        e_n = .dev. history%state_n%strain_e
        hist%vals = 0D0
        do i = 1, self%n
            hist%vals = hist%vals + coefs%decay(i)*history%state_n%h(:, i) - coefs%weight(i)*e_n%vals
        end do
        call elas%set(coefs%k_alg, coefs%g_alg, hist)
    end function freeze_maxwell

    pure subroutine update_maxwell(self, coefs, strain_e, history)
        class(Generalized_maxwell), intent(in)    :: self
        type(Visco_coefficients), intent(in)      :: coefs
        type(ten_3D2Osym), intent(in)             :: strain_e
        type(Viscoelastic_history), intent(inout) :: history
        type(ten_3D2Osym) :: de
        integer :: i

        de = (.dev. strain_e) - (.dev. history%state_n%strain_e)
        history%state_np1%n_branches = self%n
        do i = 1, self%n
            history%state_np1%h(:, i) = coefs%decay(i)*history%state_n%h(:, i) + coefs%weight(i)*de%vals
        end do
        history%state_np1%strain_e = strain_e
    end subroutine update_maxwell

end module muscle_visco_maxwell
