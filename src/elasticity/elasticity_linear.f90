module mod_elasticity_linear
    !! Module mod_elasticity_linear
    !! ============================
    !!
    !! Defines a concrete implementation for linear isotropic elasticity.
    !!
    !! This module provides the `Elasticity_linear` derived type, which represents a
    !! standard linear isotropic elastic material model. It extends the abstract
    !! `Base_elasticity` type defined in [[mod_base_elasticity]].
    !!
    !! The model calculates stress based on Hooke's law for isotropic materials using
    !! Young's modulus (E) and Poisson's ratio (nu). It also provides the constant
    !! fourth-order elasticity tensor (tangent modulus) corresponding to these parameters.
    !!
    !! This type overrides the deviatoric procedures `stress_dev_3D` and
    !! `dstress_dstrain_dev_3D` inherited from `Base_elasticity` with direct analytic
    !! implementations that avoid the general projector computation.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `Elasticity_linear`: Concrete type for linear isotropic elasticity.
    !!     - Extends: `Base_elasticity`.
    !!     - Components: Lamé parameters (`lambda`, `mu`), the diagonal constant
    !!       `e1 = lambda + 2*mu`, pre-cached stiffness tensors (`tan_3D`, `tan_2D`),
    !!       and flags indicating whether each cache is populated.
    !!     - `set_parameters(young, poisson, [calc_tan_3D], [calc_tan_2D])` — initialises
    !!       the material and optionally pre-calculates the stiffness tensors (default: yes).
    !!     - `stress_3D(strain)` — full 3-D stress via Hooke's law.
    !!     - `stress_2D(strain)` — full 2-D stress via Hooke's law.
    !!     - `stress_dev_3D(strain)` — deviatoric 3-D stress, analytically computed
    !!       as \( s_{ij} = 2\mu\,\varepsilon_{ij} - \tfrac{2}{3}\mu\,\varepsilon_{kk}\delta_{ij} \).
    !!       **Overrides** the default in `Base_elasticity`.
    !!     - `dstress_dstrain_3D(strain)` — full 3-D tangent modulus (constant).
    !!     - `dstress_dstrain_2D(strain)` — full 2-D tangent modulus (constant).
    !!     - `dstress_dstrain_dev_3D(strain)` — deviatoric 3-D tangent modulus, analytically
    !!       computed as \( C^{\mathrm{dev}}_{ijkl} = 2\mu\left(\mathbb{I}^S_{ijkl} -
    !!       \tfrac{1}{3}\delta_{ij}\delta_{kl}\right) \).
    !!       **Overrides** the default in `Base_elasticity`.
    !!
    !! Mathematical Background
    !! -----------------------
    !!
    !! Hooke's law for isotropic media:
    !! \[ \sigma_{ij} = \lambda\,\delta_{ij}\,\varepsilon_{kk} + 2\mu\,\varepsilon_{ij} \]
    !!
    !! Lamé parameters from engineering constants:
    !! \[ \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}, \qquad \mu = \frac{E}{2(1+\nu)} \]
    !!
    !! The deviatoric stress is:
    !! \[ s_{ij} = 2\mu\,\varepsilon_{ij} - \tfrac{2}{3}\mu\,\varepsilon_{kk}\,\delta_{ij} \]
    !!
    !! The full tangent modulus:
    !! \[ C_{ijkl} = \lambda\,\delta_{ij}\delta_{kl} + \mu(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}) \]
    !!
    !! The deviatoric tangent modulus:
    !! \[ C^{\mathrm{dev}}_{ijkl} = 2\mu\left(\mathbb{I}^S_{ijkl} - \tfrac{1}{3}\delta_{ij}\delta_{kl}\right) \]
    !!
    !! In Voigt notation this gives diagonal normal components \( 4\mu/3 \),
    !! off-diagonal normal components \( -2\mu/3 \), and shear components \( \mu \).
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_linear_elasticity_usage
    !!   use mod_elasticity_linear
    !!   use tensors_types
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(Elasticity_linear) :: steel_material
    !!   type(ten_3D2Osym) :: strain_tensor, stress_tensor, stress_dev_tensor
    !!   type(ten_3D4O3sym) :: stiffness_tensor, stiffness_dev_tensor
    !!   real(real64) :: E, nu
    !!
    !!   E  = 200000.0D0   ! Young's Modulus [MPa]
    !!   nu = 0.3D0        ! Poisson's Ratio
    !!
    !!   call steel_material%set_parameters(young=E, poisson=nu)
    !!
    !!   call strain_tensor%init(xx=0.001D0, yy=-0.0003D0, zz=-0.0003D0, &
    !!                           xy=0.0005D0, yz=0.0D0, xz=0.0D0)
    !!
    !!   stress_tensor     = steel_material%stress(strain_tensor)
    !!   stress_dev_tensor = steel_material%stress_dev(strain_tensor)
    !!   stiffness_tensor  = steel_material%dstress_dstrain(strain_tensor)
    !!   stiffness_dev_tensor = steel_material%dstress_dstrain_dev(strain_tensor)
    !!
    !! end program example_linear_elasticity_usage
    !! ```
    !!
    !! For the base class definition see [[mod_base_elasticity]].
    !! For tensor type definitions see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use tensors_types
    use mod_base_elasticity
    implicit none
    PRIVATE

    PUBLIC :: Elasticity_linear
    type, extends(Base_elasticity) :: Elasticity_linear
        !! Concrete type for Linear Isotropic Elasticity.
        !!
        !! Extends `Base_elasticity` and provides analytic implementations for all
        !! stress and tangent modulus evaluations, including direct deviatoric variants
        !! that override the general projector-based defaults of `Base_elasticity`.
        real(real64), private :: lambda = 0.0D0
            !! Lamé's first parameter \( \lambda = E\nu/[(1+\nu)(1-2\nu)] \).
        real(real64), private :: mu = 0.0D0
            !! Shear modulus \( \mu = E/[2(1+\nu)] \).
        real(real64), private :: e1 = 0.0D0
            !! Diagonal constant \( e_1 = \lambda + 2\mu \) used in Hooke's law.
        type(ten_3D4O3sym), private :: tan_3D
            !! Pre-cached 3-D full stiffness tensor \( \mathbf{C} \).
            !! Populated by `set_parameters` when `calc_tan_3D = .TRUE.` (default).
        type(ten_2D4O3sym), private :: tan_2D
            !! Pre-cached 2-D full stiffness tensor \( \mathbf{C} \).
            !! Populated by `set_parameters` when `calc_tan_2D = .TRUE.` (default).
        logical, private :: tan_init_3D = .FALSE.
            !! `.TRUE.` if `tan_3D` has been populated and can be returned directly.
        logical, private :: tan_init_2D = .FALSE.
            !! `.TRUE.` if `tan_2D` has been populated and can be returned directly.
    contains
        procedure :: set_parameters
            !! Initialises \( \lambda \), \( \mu \), and optionally pre-caches the
            !! stiffness tensors. Must be called before any stress evaluation.
        procedure :: stress_3D          => stress_linear_3D
            !! Pure function. Full 3-D Cauchy stress via Hooke's law.
        procedure :: stress_2D          => stress_linear_2D
            !! Pure function. Full 2-D Cauchy stress via Hooke's law.
        procedure :: stress_dev_3D      => stress_linear_dev_3D
            !! Pure function. Deviatoric 3-D stress computed analytically.
            !! Overrides the `.dev.`-based default of `Base_elasticity`.
        procedure :: dstress_dstrain_3D  => dstress_dstrain_linear_3D
            !! Pure function. Full 3-D tangent modulus (constant for linear elasticity).
        procedure :: dstress_dstrain_2D  => dstress_dstrain_linear_2D
            !! Pure function. Full 2-D tangent modulus (constant for linear elasticity).
        procedure :: dstress_dstrain_dev_3D => dstress_dstrain_dev_linear_3D
            !! Pure function. Deviatoric 3-D tangent modulus computed analytically.
            !! Overrides the projector-based default of `Base_elasticity`.
    end type Elasticity_linear

contains

    pure subroutine set_parameters(self, young, poisson, calc_tan_3D, calc_tan_2D)
        !! Initialises the linear elastic material with Young's modulus and Poisson's ratio.
        !!
        !! Computes the Lamé parameters \( \lambda \) and \( \mu \), the derived constant
        !! \( e_1 = \lambda + 2\mu \), and optionally pre-calculates and caches the full
        !! isotropic stiffness tensors for 3-D and 2-D problems.
        !!
        !! Pre-calculation is enabled by default (`calc_tan_3D = .TRUE.`,
        !! `calc_tan_2D = .TRUE.`) and is recommended when `dstress_dstrain` or
        !! `dstress_dstrain_dev` will be called repeatedly, since it avoids redundant
        !! assembly of the stiffness tensor at each call.
        implicit none
        class(Elasticity_linear), intent(inout) :: self
            !! The linear elasticity model object to initialise.
        real(real64), intent(in) :: young
            !! Young's modulus \( E \).
        real(real64), intent(in) :: poisson
            !! Poisson's ratio \( \nu \).
        logical, optional, intent(in) :: calc_tan_3D
            !! If `.TRUE.` (default), pre-calculates and caches the 3-D stiffness tensor.
        logical, optional, intent(in) :: calc_tan_2D
            !! If `.TRUE.` (default), pre-calculates and caches the 2-D stiffness tensor.
        real(real64) :: e1, lambda, mu
        logical :: ccalc_tan_3D, ccalc_tan_2D

        ccalc_tan_3D = .TRUE.
        if (present(calc_tan_3D)) ccalc_tan_3D = calc_tan_3D

        ccalc_tan_2D = .TRUE.
        if (present(calc_tan_2D)) ccalc_tan_2D = calc_tan_2D

        ! Compute Lamé parameters
        lambda = young * poisson / ((1.0D0 + poisson) * (1.0D0 - 2.0D0 * poisson))
        mu     = 0.5D0 * young / (1.0D0 + poisson)

        self%lambda = lambda
        self%mu     = mu
        self%e1     = lambda + 2.0D0 * mu   ! = lambda + 2*mu

        ! Pre-cache 3-D full stiffness tensor C
        self%tan_init_3D = ccalc_tan_3D
        if (ccalc_tan_3D) then
            e1 = self%e1
            call self%tan_3D%init( xxxx=e1,     yyyy=e1,     zzzz=e1,      &
                                   xxyy=lambda, yyzz=lambda, xxzz=lambda,  &
                                   xxxy=0.0D0,  xxyz=0.0D0,  xxxz=0.0D0,  &
                                   yyxy=0.0D0,  yyyz=0.0D0,  yyxz=0.0D0,  &
                                   zzxy=0.0D0,  zzyz=0.0D0,  zzxz=0.0D0,  &
                                   xyxy=mu,     yzyz=mu,     xzxz=mu,      &
                                   xyyz=0.0D0,  yzxz=0.0D0,  xyxz=0.0D0   )
        end if

        ! Pre-cache 2-D full stiffness tensor C
        self%tan_init_2D = ccalc_tan_2D
        if (ccalc_tan_2D) then
            e1 = self%e1
            call self%tan_2D%init( xxxx=e1,     yyyy=e1,     zzzz=e1,      &
                                   xxyy=lambda, yyzz=lambda, xxzz=lambda,  &
                                   xxxy=0.0D0,  yyxy=0.0D0,  zzxy=0.0D0,  &
                                   xyxy=mu                                  )
        end if
    end subroutine set_parameters


    ! =========================================================================
    ! Full stress
    ! =========================================================================

    pure function stress_linear_3D(self, strain) result(res)
        !! Computes the full 3-D Cauchy stress tensor using Hooke's law:
        !! \[ \sigma_{ij} = \lambda\,\varepsilon_{kk}\,\delta_{ij} + 2\mu\,\varepsilon_{ij} \]
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object (read-only).
        class(ten_3D2Osym), intent(in) :: strain
            !! Input second-order symmetric strain tensor (`ten_3D2Osym`).
        type(ten_3D2Osym) :: res
            !! Output second-order symmetric Cauchy stress tensor (`ten_3D2Osym`).
        real(real64) :: e1, lam, mu, tr

        e1  = self%e1       ! lambda + 2*mu
        lam = self%lambda
        mu  = self%mu

        tr = strain%xx() + strain%yy() + strain%zz()   ! volumetric strain

        call res%init( xx = e1*strain%xx() + lam*(strain%yy() + strain%zz()), &
                       yy = e1*strain%yy() + lam*(strain%zz() + strain%xx()), &
                       zz = e1*strain%zz() + lam*(strain%xx() + strain%yy()), &
                       xy = 2.0D0*mu*strain%xy(), &
                       yz = 2.0D0*mu*strain%yz(), &
                       xz = 2.0D0*mu*strain%xz()  )
    end function stress_linear_3D


    pure function stress_linear_2D(self, strain) result(res)
        !! Computes the full 2-D Cauchy stress tensor using Hooke's law.
        !! See `stress_linear_3D` for the mathematical description.
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object (read-only).
        class(ten_2D2Osym), intent(in) :: strain
            !! Input second-order symmetric strain tensor (`ten_2D2Osym`).
        type(ten_2D2Osym) :: res
            !! Output second-order symmetric Cauchy stress tensor (`ten_2D2Osym`).
        real(real64) :: e1, lam, mu

        e1  = self%e1
        lam = self%lambda
        mu  = self%mu

        call res%init( xx = e1*strain%xx() + lam*(strain%yy() + strain%zz()), &
                       yy = e1*strain%yy() + lam*(strain%zz() + strain%xx()), &
                       zz = e1*strain%zz() + lam*(strain%xx() + strain%yy()), &
                       xy = 2.0D0*mu*strain%xy() )
    end function stress_linear_2D


    ! =========================================================================
    ! Deviatoric stress — analytic override
    ! =========================================================================

    pure function stress_linear_dev_3D(self, strain) result(res)
        !! Computes the deviatoric part of the 3-D Cauchy stress analytically:
        !!
        !! \[ s_{ij} = 2\mu\,\varepsilon_{ij}
        !!            - \tfrac{2}{3}\mu\,\varepsilon_{kk}\,\delta_{ij} \]
        !!
        !! This is equivalent to \( \mathrm{dev}(\boldsymbol{\sigma}) \) for an
        !! isotropic linear elastic material, since the volumetric part of \( \boldsymbol{\sigma} \)
        !! is \( \lambda\,\varepsilon_{kk}\,\mathbf{I} \), which is purely spherical.
        !!
        !! Overrides the default `.dev.`-based implementation in `Base_elasticity`
        !! with a direct evaluation that avoids computing the full stress first.
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object (read-only).
        class(ten_3D2Osym), intent(in) :: strain
            !! Input second-order symmetric strain tensor (`ten_3D2Osym`).
        type(ten_3D2Osym) :: res
            !! Output deviatoric stress tensor (`ten_3D2Osym`).
        real(real64) :: mu, two_mu, tr_comp

        mu      = self%mu
        two_mu  = 2.0D0 * mu
        tr_comp = two_mu / 3.0D0 * (strain%xx() + strain%yy() + strain%zz())
            ! = (2mu/3) * tr(eps) — subtracted from each normal component

        call res%init( xx = two_mu*strain%xx() - tr_comp, &
                       yy = two_mu*strain%yy() - tr_comp, &
                       zz = two_mu*strain%zz() - tr_comp, &
                       xy = two_mu*strain%xy(),            &
                       yz = two_mu*strain%yz(),            &
                       xz = two_mu*strain%xz()             )
    end function stress_linear_dev_3D


    ! =========================================================================
    ! Full tangent modulus
    ! =========================================================================

    pure function dstress_dstrain_linear_3D(self, strain) result(res)
        !! Returns the full 3-D isotropic tangent modulus (constant for linear elasticity):
        !! \[ C_{ijkl} = \lambda\,\delta_{ij}\delta_{kl}
        !!              + \mu(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}) \]
        !!
        !! Returns the pre-cached tensor `tan_3D` if available (populated by
        !! `set_parameters`); otherwise assembles it on the fly.
        !!
        !! The input `strain` is ignored since the tangent is strain-independent.
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object (read-only).
        class(ten_3D2Osym), intent(in) :: strain
            !! Input strain tensor. Ignored; retained for interface consistency.
        type(ten_3D4O3sym) :: res
            !! Output fourth-order tangent modulus tensor (`ten_3D4O3sym`).
        real(real64) :: e1, lam, mu

        if (self%tan_init_3D) then
            res = self%tan_3D
            return
        end if

        e1  = self%e1
        lam = self%lambda
        mu  = self%mu
        call res%init( xxxx=e1,   yyyy=e1,   zzzz=e1,    &
                       xxyy=lam,  yyzz=lam,  xxzz=lam,   &
                       xxxy=0.0D0, xxyz=0.0D0, xxxz=0.0D0, &
                       yyxy=0.0D0, yyyz=0.0D0, yyxz=0.0D0, &
                       zzxy=0.0D0, zzyz=0.0D0, zzxz=0.0D0, &
                       xyxy=mu,   yzyz=mu,   xzxz=mu,    &
                       xyyz=0.0D0, yzxz=0.0D0, xyxz=0.0D0  )
    end function dstress_dstrain_linear_3D


    pure function dstress_dstrain_linear_2D(self, strain) result(res)
        !! Returns the full 2-D isotropic tangent modulus.
        !! See `dstress_dstrain_linear_3D` for the mathematical description.
        !!
        !! Returns the pre-cached tensor `tan_2D` if available; otherwise assembles
        !! it on the fly. The input `strain` is ignored.
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object (read-only).
        class(ten_2D2Osym), intent(in) :: strain
            !! Input strain tensor. Ignored; retained for interface consistency.
        type(ten_2D4O3sym) :: res
            !! Output fourth-order tangent modulus tensor (`ten_2D4O3sym`).
        real(real64) :: e1, lam, mu

        if (self%tan_init_2D) then
            res = self%tan_2D
            return
        end if

        e1  = self%e1
        lam = self%lambda
        mu  = self%mu
        call res%init( xxxx=e1,   yyyy=e1,   zzzz=e1,   &
                       xxyy=lam,  yyzz=lam,  xxzz=lam,  &
                       xxxy=0.0D0, yyxy=0.0D0, zzxy=0.0D0, &
                       xyxy=mu                             )
    end function dstress_dstrain_linear_2D


    ! =========================================================================
    ! Deviatoric tangent modulus — analytic override
    ! =========================================================================

    pure function dstress_dstrain_dev_linear_3D(self, strain) result(res)
        !! Returns the deviatoric part of the 3-D isotropic tangent modulus analytically:
        !!
        !! \[ C^{\mathrm{dev}}_{ijkl} = 2\mu\left(\mathbb{I}^S_{ijkl}
        !!    - \tfrac{1}{3}\,\delta_{ij}\delta_{kl}\right) \]
        !!
        !! In Voigt notation this gives:
        !! - Normal diagonal components \( (C^{\mathrm{dev}}_{iiii}) \): \( 4\mu/3 \)
        !! - Normal off-diagonal components \( (C^{\mathrm{dev}}_{iijj},\,i\neq j) \): \( -2\mu/3 \)
        !! - Shear components \( (C^{\mathrm{dev}}_{ijij},\,i\neq j) \): \( \mu \)
        !!
        !! Overrides the projector-based default of `Base_elasticity` with a direct
        !! evaluation. The pre-cached `tan_3D` is NOT used here because it stores the
        !! full (non-deviatoric) stiffness.
        !!
        !! The input `strain` is ignored since the tangent is strain-independent.
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object (read-only).
        class(ten_3D2Osym), intent(in) :: strain
            !! Input strain tensor. Ignored; retained for interface consistency.
        type(ten_3D4O3sym) :: res
            !! Output deviatoric fourth-order tangent modulus tensor (`ten_3D4O3sym`).
        real(real64) :: mu, c_diag, c_off

        mu     = self%mu
        c_diag =  4.0D0/3.0D0 * mu   ! normal diagonal:     4mu/3
        c_off  = -2.0D0/3.0D0 * mu   ! normal off-diagonal: -2mu/3

        call res%init( xxxx=c_diag, yyyy=c_diag, zzzz=c_diag,  &
                       xxyy=c_off,  yyzz=c_off,  xxzz=c_off,   &
                       xxxy=0.0D0,  xxyz=0.0D0,  xxxz=0.0D0,   &
                       yyxy=0.0D0,  yyyz=0.0D0,  yyxz=0.0D0,   &
                       zzxy=0.0D0,  zzyz=0.0D0,  zzxz=0.0D0,   &
                       xyxy=mu,     yzyz=mu,     xzxz=mu,       &
                       xyyz=0.0D0,  yzxz=0.0D0,  xyxz=0.0D0    )
    end function dstress_dstrain_dev_linear_3D

end module