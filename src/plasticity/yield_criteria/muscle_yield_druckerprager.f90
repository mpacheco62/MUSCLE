! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_druckerprager
    !! Module muscle_yield_druckerprager
    !! =================================
    !! Implements the Extended Linear Drucker-Prager Yield Criterion matching Abaqus formulation.
    !!
    !! Mathematical Formulation
    !! ------------------------
    !! The yield function is defined as:
    !! \[ F = t + \tan(\beta) \cdot f \cdot p - Y = 0 \]
    !! where:
    !! - \( p = \frac{1}{3} \mathrm{tr}(\boldsymbol{\sigma}) \) is the hydrostatic stress.
    !! - \( q = \sqrt{\frac{3}{2} \mathbf{S}:\mathbf{S}} \) is the Mises equivalent stress.
    !! - \( t = \frac{1}{2} q \left[ 1 + \frac{1}{K} - \left(1 - \frac{1}{K}\right) \left(\frac{r}{q}\right)^3 \right] \)
    !! - \( K \) is the ratio of yield stress in triaxial tension to triaxial compression (\(0.778 \le K \le 1.0\)).

    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_base
    implicit none
    private

    integer, parameter, public :: DP_HARDENING_TENSION     = 1
    integer, parameter, public :: DP_HARDENING_COMPRESSION = 2
    integer, parameter, public :: DP_HARDENING_SHEAR       = 3

    public :: DruckerPrager

    type, extends(Base_yield_critera) :: DruckerPrager
        !! Extended Linear Drucker-Prager Yield Criterion (Abaqus Compatible)
        private
        real(real64) :: f_factor = 0.0D0
        real(real64) :: I1_factor = 0.0D0
        real(real64) :: K = 0.0D0
        integer      :: hardening_mode = DP_HARDENING_TENSION
    contains
        procedure, public :: init => init_dp
        procedure, public :: stress_eq => stress_eq_dp
        procedure, public :: dstressEq_dstress => dstressEq_dstress_dp
    end type DruckerPrager

contains

    pure subroutine init_dp(self, beta_deg, K, hardening_mode)
        !! Initializes Drucker-Prager parameters matching Abaqus material options.
        implicit none
        class(DruckerPrager), intent(inout) :: self
        real(real64), intent(in)          :: beta_deg       !! Friction angle beta in degrees
        real(real64), intent(in)          :: K              !! Yield stress ratio K (0.778 <= K <= 1.0)
        integer, intent(in)               :: hardening_mode !! Hardening calibration mode

        real(real64), parameter :: PI = acos(-1.0D0)
        real(real64) :: tanbeta

        tanbeta = tan(beta_deg * PI / 180.0D0)
        select case (hardening_mode)
        case(DP_HARDENING_TENSION)
            self%f_factor = (1.0D0 / K + tanbeta / 3.0D0) ** (-1)
        case(DP_HARDENING_COMPRESSION)
            self%f_factor = (1.0D0 - tanbeta / 3.0D0) ** (-1)
        case(DP_HARDENING_SHEAR)
            self%f_factor = 1.0D0
        end select

        self%I1_factor = tanbeta * self%f_factor
        self%K = K
        self%hardening_mode = hardening_mode

    end subroutine init_dp


    pure function stress_eq_dp(self, stress) result(res)
        !! Computes the Drucker-Prager equivalent stress.
        implicit none
        class(DruckerPrager), intent(in) :: self
        type(ten_3D2Osym), intent(in)   :: stress
        real(real64)                     :: res

        real(real64) :: q, r3, t, p, stress_norm
        type(ten_3D2Osym) :: s
        real(real64), parameter :: TOL_REL = 1.0D-9
        real(real64), parameter :: TOL_ABS = 1.0D-40

        stress_norm = sqrt(stress .ddot. stress)

        s = .dev. stress
        q = sqrt(1.5D0 * (s .ddot. s))
        p = (1.0D0 / 3.0D0) * (stress%xx() + stress%yy() + stress%zz())

        ! Zero stress check
        if (stress_norm < TOL_ABS) then
            res = 0.0D0
            return
        end if

        ! Hydrostatic singularity check (q -> 0)
        if ((q / stress_norm) < TOL_REL) then
            res = self%I1_factor * p
            return
        end if

        ! Third invariant term r^3 = (9/2) * tr(S^3)
        r3 = (9.0D0 / 2.0D0) * ((s * s) .ddot. s)

        ! Abaqus equivalent deviatoric stress t
        t = 0.5D0 * q * (1.0D0 + 1.0D0 / self%K - (1.0D0 - 1.0D0 / self%K) * (r3 / (q**3)))
        t = t * self%f_factor

        res = t + self%I1_factor * p

    end function stress_eq_dp

    pure function dstressEq_dstress_dp(self, stress) result(N_tensor)
        !! Analytical First Derivative (Gradient / Flow Direction): N = d(sigma_eq) / d(sigma)
        implicit none
        class(DruckerPrager), intent(in) :: self
        type(ten_3D2Osym), intent(in)    :: stress
        type(ten_3D2Osym)                :: N_tensor

        type(ten_3D2Osym) :: s, s2
        real(real64)      :: q, r3, c1, c2, f, d
        real(real64)      :: alpha1, alpha2, alpha3, stress_norm
        type(iden_2O)     :: I2O

        real(real64), parameter :: TOL_REL = 1.0D-9
        real(real64), parameter :: TOL_ABS = 1.0D-40

        f = self%f_factor
        d = self%I1_factor

        stress_norm = sqrt(stress .ddot. stress)

        ! Handle Zero Stress State
        if (stress_norm < TOL_ABS) then
            N_tensor = (d / 3.0D0) * I2O
            return
        end if

        s = .dev. stress
        q = sqrt(1.5D0 * (s .ddot. s))

        ! Handle Hydrostatic State (q -> 0)
        if ((q / stress_norm) < TOL_REL) then
            N_tensor = (d / 3.0D0) * I2O
            return
        end if

        ! Invariant constants
        c1 = 0.5D0 * (1.0D0 + 1.0D0 / self%K)
        c2 = 0.5D0 * (1.0D0 - 1.0D0 / self%K)

        s2 = s%square()
        r3 = 4.5D0 * (s2 .ddot. s)

        ! Scalar coefficients for gradient
        alpha1 = (1.5D0 / q) * (c1 + 2.0D0 * c2 * (r3 / (q**3)))
        alpha2 = -13.5D0 * c2 / (q**2)
        alpha3 = (d / 3.0D0) + 3.0D0 * f * c2

        ! N = f * (alpha1 * S + alpha2 * S^2) + alpha3 * I
        N_tensor = f * (alpha1 * s + alpha2 * s2) + alpha3 * I2O

    end function dstressEq_dstress_dp

end module muscle_yield_druckerprager