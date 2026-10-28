! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_vp_peric
    !! Direct Peric overstress: de Souza Neto, Peric and Owen (2008), Eq. (11.86).
    !! g = ((sigma_eq/Y(ep))**(1/eps) - 1)/mu above yield; zero otherwise.
    !! g and its derivatives are 0 for sigma_eq <= Y by convention; the right-hand
    !! derivative at the threshold is not 0.
    use, intrinsic :: iso_fortran_env, only: real64
    use muscle_vp_overstress_base, only: Base_overstress_law
    use muscle_hard_base, only: Base_hardening_laws
    implicit none
    private
    public :: Peric_overstress

    type, extends(Base_overstress_law) :: Peric_overstress
        class(Base_hardening_laws), allocatable :: hard_law
        real(real64) :: mu = 1D0   !! Positive relaxation time mu
        real(real64) :: eps = 1D0  !! Positive rate sensitivity epsilon
    contains
        procedure :: flow_rate => flow_rate_peric
        procedure :: dflow_dsigma => dflow_dsigma_peric
        procedure :: dflow_dep => dflow_dep_peric
    end type Peric_overstress

contains
    pure function dflow_dr(self, r) result(res)
        ! Slope in r=sigma_eq/Y; de Souza Neto et al. (2008), Eq. (11.86), chain rule.
        ! Called by the two analytical derivative bindings.
        implicit none
        class(Peric_overstress), intent(in) :: self
        real(real64), intent(in) :: r
        real(real64) :: res
        if (r <= 1D0) then
            res = 0D0
        else
            res = r**(1D0/self%eps - 1D0)/(self%mu*self%eps)
        end if
    end function dflow_dr

    pure function flow_rate_peric(self, sigma_eq, ep) result(res)
        !! Flow rate for Closest_point; de Souza Neto et al. (2008), Eq. (11.86).
        implicit none
        class(Peric_overstress), intent(in) :: self
        real(real64), intent(in) :: sigma_eq, ep
        real(real64) :: res, r  ! r: stress ratio, assumes Y > 0
        r = sigma_eq/self%hard_law%stress(ep)
        if (r <= 1D0) then
            res = 0D0
        else
            res = (r**(1D0/self%eps) - 1D0)/self%mu
        end if
    end function flow_rate_peric

    pure function dflow_dsigma_peric(self, sigma_eq, ep) result(res)
        !! Stress sensitivity for Closest_point; de Souza Neto et al. (2008), Eq. (11.86), chain rule.
        implicit none
        class(Peric_overstress), intent(in) :: self
        real(real64), intent(in) :: sigma_eq, ep
        real(real64) :: res, y  ! y: current static yield stress Y(ep)
        y = self%hard_law%stress(ep)
        res = dflow_dr(self, sigma_eq/y)/y
    end function dflow_dsigma_peric

    pure function dflow_dep_peric(self, sigma_eq, ep) result(res)
        !! Hardening sensitivity for Closest_point; de Souza Neto et al. (2008), Eq. (11.86), chain rule.
        implicit none
        class(Peric_overstress), intent(in) :: self
        real(real64), intent(in) :: sigma_eq, ep
        real(real64) :: res, y, r  ! y=Y(ep); r=sigma_eq/y
        y = self%hard_law%stress(ep)
        r = sigma_eq/y
        res = -dflow_dr(self, r)*r*self%hard_law%dstress_dep(ep)/y
    end function dflow_dep_peric
end module muscle_vp_peric
