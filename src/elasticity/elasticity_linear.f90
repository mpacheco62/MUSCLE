module mod_elasticity_linear
    !! Module mod_elasticity_linear
    !! ============================
    !!
    !! Defines a concrete implementation for linear isotropic elasticity.
    !!
    !! This module provides the `Elasticity_linear` derived type, which represents a
    !! standard linear isotropic elastic material model. It extends the abstract
    !! `Base_elasticity` type defined in `mod_base_elasticity`.
    !!
    !! The model calculates stress based on Hooke's law for isotropic materials using
    !! Young's modulus (E) and Poisson's ratio (nu). It also provides the constant
    !! fourth-order elasticity tensor (tangent modulus) corresponding to these parameters.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `Elasticity_linear`: Concrete type representing a linear isotropic elastic material.
    !!     - Extends: `Base_elasticity`.
    !!     - Private Components: Stores Young's modulus (`young`), Poisson's ratio (`poisson`),
    !!       derived elastic constants (`e1`, `e2`, `e3`), the pre-calculated tangent
    !!       modulus (`tan`), and a flag indicating if the tangent is initialized (`tan_init`).
    !!     - Procedure: `set_parameters(young, poisson, [calc_tan])` - Initializes the material
    !!       properties (Young's modulus and Poisson's ratio). Optionally pre-calculates
    !!       and stores the tangent modulus tensor if `calc_tan` is present and `.true.` (default).
    !!     - Procedure: `stress(strain)` - Calculates the stress tensor (`ten_3D2Osym`)
    !!       corresponding to the input `strain` tensor (`ten_3D2Osym`) using Hooke's law.
    !!     - Procedure: `dstress_dstrain(strain)` - Returns the constant linear elastic
    !!       tangent modulus tensor (`ten_3D4O3sym`) for the defined material properties.
    !!       The input `strain` is ignored as the tangent is constant for linear elasticity.
    !!
    !! Mathematical Background
    !! -----------------------
    !!
    !! The stress (\(\sigma\)) is calculated from strain (\(\epsilon\)) using Hooke's Law:
    !! \[ \sigma_{ij} = \lambda \delta_{ij} \epsilon_{kk} + 2 \mu \epsilon_{ij} \]
    !! where \(\lambda\) and \(\mu\) are Lamé's parameters, derived from Young's modulus \(E\)
    !! and Poisson's ratio \(\nu\):
    !! \[ \lambda = \frac{E \nu}{(1+\nu)(1-2\nu)} \]
    !! \[ \mu = G = \frac{E}{2(1+\nu)} \]
    !! The implementation uses equivalent constants `e1`, `e2`, `e3` derived from E and nu.
    !!
    !! The fourth-order elasticity tensor \(C_{ijkl} = \frac{\partial \sigma_{ij}}{\partial \epsilon_{kl}}\) is constant:
    !! \[ C_{ijkl} = \lambda \delta_{ij} \delta_{kl} + \mu (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}) \]
    !! This tensor is stored internally in Voigt notation (`ten_3D4O3sym`).
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
    !!   type(ten_3D2Osym) :: strain_tensor, stress_tensor
    !!   type(ten_3D4O3sym) :: stiffness_tensor
    !!   real(real64) :: E, nu
    !!
    !!   ! Define material properties
    !!   E = 200000.0D0 ! Young's Modulus (e.g., MPa)
    !!   nu = 0.3       ! Poisson's Ratio
    !!
    !!   ! Initialize the material model (tangent is calculated by default)
    !!   call steel_material%set_parameters(young=E, poisson=nu)
    !!
    !!   ! Define a strain state
    !!   call strain_tensor%init(xx=0.001, yy=-0.0003, zz=-0.0003, xy=0.0005, yz=0.0, xz=0.0)
    !!
    !!   ! Calculate stress
    !!   stress_tensor = steel_material%stress(strain_tensor)
    !!
    !!   ! Get the tangent modulus (stiffness)
    !!   stiffness_tensor = steel_material%dstress_dstrain(strain_tensor) ! strain is ignored here
    !!
    !!   print *, "Steel Properties:"
    !!   print *, "  Young's Modulus (E) =", E
    !!   print *, "  Poisson's Ratio (nu) =", nu
    !!   print *, "Calculated Stress:"
    !!   print *, "  Stress XX =", stress_tensor%xx()
    !!   print *, "  Stress XY =", stress_tensor%xy()
    !!   print *, "Stiffness Tensor:"
    !!   print *, "  C_1111 =", stiffness_tensor%vals(1)
    !!   print *, "  C_1122 =", stiffness_tensor%vals(7)
    !!   print *, "  C_1212 =", stiffness_tensor%vals(4)
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
        !! Extends `Base_elasticity` and implements the stress and tangent modulus calculations
        !! based on Young's modulus and Poisson's ratio.
        real(real64), private :: young        !! Young's Modulus (E)
        real(real64), private :: poisson      !! Poisson's Ratio (nu)
        real(real64), private :: e1, e2, e3   !! Derived elastic constants related to Lame parameters
        type(ten_3D4O3sym), private :: tan_3D    !! Pre-calculated tangent modulus tensor (stiffness)
        type(ten_2D4O3sym), private :: tan_2D    !! Pre-calculated tangent modulus tensor (stiffness)
        logical, private :: tan_init_3D=.FALSE.  !! Flag indicating if 'tan' has been calculated
        logical, private :: tan_init_2D=.FALSE.  !! Flag indicating if 'tan' has been calculated
    contains
        procedure :: dstress_dstrain_3D => dstress_dstrain_linear_3D
            !! Calculates the constant tangent modulus tensor.
        procedure :: dstress_dstrain_2D => dstress_dstrain_linear_2D
            !! Calculates the constant tangent modulus tensor.
        procedure :: stress_3D => stress_linear_3D
            !! Calculates stress using Hooke's law.
        procedure :: stress_2D => stress_linear_2D
            !! Calculates stress using Hooke's law.
        procedure :: set_parameters
            !! Sets the material parameters (E, nu) and optionally pre-calculates the tangent.
    end type Elasticity_linear

    contains
    pure subroutine set_parameters(self, young, poisson, calc_tan_3D, calc_tan_2D)
        !! Sets the material parameters (Young's modulus, Poisson's ratio) for the linear elastic model.
        !! Optionally pre-calculates the tangent modulus tensor.
        implicit none
        logical, optional, intent(in) :: calc_tan_3D, calc_tan_2D
            !! If present and .TRUE. (default), pre-calculates and stores the tangent modulus.
        class(Elasticity_linear), intent(inout) :: self
            !! The linear elasticity model object
        real(real64), intent(in) :: young
            !! Young's Modulus (E).
        real(real64), intent(in) :: poisson
            !! Poisson's Ratio (nu).
        real(real64) :: e1, e2, e3
        logical :: ccalc_tan_3D, ccalc_tan_2D
        
        ccalc_tan_3D = .True.
        if (present(calc_tan_3D)) ccalc_tan_3D = calc_tan_3D

        ccalc_tan_2D = .True.
        if (present(calc_tan_2D)) ccalc_tan_2D = calc_tan_2D

        self%young = young
        self%poisson = poisson

        ! Calculate derived constants (related to Lame parameters lambda and mu)
        ! e1 = lambda + 2*mu
        ! e2 = lambda
        ! e3 = mu (Shear Modulus G)
        e1 = young*(1D0-poisson)/((1D0+poisson)*(1D0-2D0*poisson))
        e2 = young*poisson/((1D0+poisson)*(1D0-2D0*poisson))
        e3 = 0.5D0*young/(1D0+poisson)
        self%e1 = e1
        self%e2 = e2
        self%e3 = e3

        self%tan_init_3D = ccalc_tan_3D
        if (ccalc_tan_3D) then
            call self%tan_3D%init( xxxx= e1, yyyy= e1, zzzz= e1,  &
                                   xxyy= e2, yyzz= e2, xxzz= e2,  &
                                   xxxy=0D0, xxyz=0D0, xxxz=0D0,  &
                                   yyxy=0D0, yyyz=0D0, yyxz=0D0,  &
                                   zzxy=0D0, zzyz=0D0, zzxz=0D0,  &
                                   xyxy= e3, yzyz= e3, xzxz= e3,  &
                                   xyyz=0D0, yzxz=0D0, xyxz=0D0   &
                                   )
        end if
        
        if (ccalc_tan_2D) then
            call self%tan_2D%init( xxxx= e1, yyyy= e1, zzzz= e1,  &
                                   xxyy= e2, yyzz= e2, xxzz= e2,  &
                                   xxxy=0D0, yyxy=0D0, zzxy=0D0,  &
                                   xyxy= e3                       &
                                   )
        end if

    end subroutine

    pure function stress_linear_3D(self, strain) result(res)
        !! Calculates the stress tensor using Hooke's law for linear isotropic elasticity.
        !! sigma = C : epsilon
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object containing material parameters.
        class(ten_3D2Osym), intent(in) :: strain
            !! Input strain tensor (`ten_3D2Osym`).
        type(ten_3D2Osym) :: res
            !! Output stress tensor (`ten_3D2Osym`).
        real(real64) :: e1, e2, e3, xx, yy, zz, xy, yz, xz

        ! Retrieve derived elastic constants
        e1 = self%e1 ! lambda + 2*mu
        e2 = self%e2 ! lambda
        e3 = self%e3 ! mu

        ! Calculate stress components using Hooke's law
        xx = e1*strain%xx() + e2*(strain%yy()+strain%zz())
        yy = e1*strain%yy() + e2*(strain%zz()+strain%xx())
        zz = e1*strain%zz() + e2*(strain%xx()+strain%yy())
        xy = 2*e3*strain%xy()  ! Note: 2*mu*epsilon_xy
        yz = 2*e3*strain%yz()  ! Note: 2*mu*epsilon_xy
        xz = 2*e3*strain%xz()  ! Note: 2*mu*epsilon_xy

        ! Initialize the result tensor
        call res%init(xx=xx, yy=yy, zz=zz, xy=xy, yz=yz, xz=xz)
        return 
    end function stress_linear_3D

    pure function stress_linear_2D(self, strain) result(res)
        !! Calculates the stress tensor using Hooke's law for linear isotropic elasticity.
        !! sigma = C : epsilon
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object containing material parameters.
        class(ten_2D2Osym), intent(in) :: strain
            !! Input strain tensor (`ten_2D2Osym`).
        type(ten_2D2Osym) :: res
            !! Output stress tensor (`ten_2D2Osym`).
        real(real64) :: e1, e2, e3, xx, yy, zz, xy

        ! Retrieve derived elastic constants
        e1 = self%e1 ! lambda + 2*mu
        e2 = self%e2 ! lambda
        e3 = self%e3 ! mu

        ! Calculate stress components using Hooke's law
        xx = e1*strain%xx() + e2*(strain%yy()+strain%zz())
        yy = e1*strain%yy() + e2*(strain%zz()+strain%xx())
        zz = e1*strain%zz() + e2*(strain%xx()+strain%yy())
        xy = 2*e3*strain%xy()  ! Note: 2*mu*epsilon_xy

        ! Initialize the result tensor
        call res%init(xx=xx, yy=yy, zz=zz, xy=xy)
        return 
    end function stress_linear_2D

    pure function dstress_dstrain_linear_3D(self, strain) result(res)
        !! Returns the constant tangent modulus (stiffness) tensor for linear isotropic elasticity.
        !! C_ijkl = d(sigma_ij) / d(epsilon_kl)
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object.
        class(ten_3D2Osym), intent(in) :: strain
            !! Input strain tensor (`ten_3D2Osym`). Ignored for linear elasticity as the tangent is constant.
        type(ten_3D4O3sym) :: res
            !! Output tangent modulus tensor (`ten_3D4O3sym`).
        real(real64) :: e1, e2, e3

        if (self%tan_init_3D) then
            ! Return the pre-calculated tangent if available
            res = self%tan_3D
            return
        else
            ! Calculate the tangent on the fly if not pre-calculated
            e1 = self%e1
            e2 = self%e2
            e3 = self%e3
            call res%init( xxxx= e1, yyyy= e1, zzzz= e1,  &
                           xxyy= e2, yyzz= e2, xxzz= e2,  &
                           xxxy=0D0, xxyz=0D0, xxxz=0D0,  &
                           yyxy=0D0, yyyz=0D0, yyxz=0D0,  &
                           zzxy=0D0, zzyz=0D0, zzxz=0D0,  &
                           xyxy= e3, yzyz= e3, xzxz= e3,  &
                           xyyz=0D0, yzxz=0D0, xyxz=0D0   &
                          )
            return
        end if
    end function dstress_dstrain_linear_3D

    pure function dstress_dstrain_linear_2D(self, strain) result(res)
        !! Returns the constant tangent modulus (stiffness) tensor for linear isotropic elasticity.
        !! C_ijkl = d(sigma_ij) / d(epsilon_kl)
        implicit none
        class(Elasticity_linear), intent(in) :: self
            !! The linear elasticity model object.
        class(ten_2D2Osym), intent(in) :: strain
            !! Input strain tensor (`ten_3D2Osym`). Ignored for linear elasticity as the tangent is constant.
        type(ten_2D4O3sym) :: res
            !! Output tangent modulus tensor (`ten_3D4O3sym`).
        real(real64) :: e1, e2, e3

        if (self%tan_init_2D) then
            ! Return the pre-calculated tangent if available
            res = self%tan_2D
            return
        else
            ! Calculate the tangent on the fly if not pre-calculated
            e1 = self%e1
            e2 = self%e2
            e3 = self%e3
            call res%init( xxxx= e1, yyyy= e1, zzzz= e1,  &
                           xxyy= e2, yyzz= e2, xxzz= e2,  &
                           xxxy=0D0, yyxy=0D0, zzxy=0D0,  &
                           xyxy= e3                       &
                          )
            return
        end if
    end function dstress_dstrain_linear_2D
end module