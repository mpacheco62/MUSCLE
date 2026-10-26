! Standard linear solid material added to MUSCLE by the XFEM integration.
module muscle_viscoelastic_sls
  use iso_fortran_env, only: real64
  use muscle_elasticity_linear, only: Elasticity_linear
  use muscle_elasticity_base, only: Base_internal_elasticity_variables
  use muscle_tensors
  implicit none
  private
  public :: Viscoelastic_sls, Sls_internal, sls_factor
  type, extends(Base_internal_elasticity_variables) :: Sls_internal
    real(real64) :: q(6)=0.0_real64
  end type
  type, extends(Elasticity_linear) :: Viscoelastic_sls
    real(real64) :: ratio=1.0_real64, tau=0.5_real64
  contains
    procedure :: to_string => sls_description
    procedure :: stress_3D_transient => sls_stress
    procedure :: dstress_dstrain_3D_transient => sls_tangent
    procedure :: stress_and_update_3D => sls_update
    procedure :: stress_dev_3D_transient => sls_dev
    procedure :: dstress_dstrain_dev_3D_transient => sls_dev_tangent
    procedure :: stress_2D_transient => unsupported_stress_2d
    procedure :: dstress_dstrain_2D_transient => unsupported_tangent_2d
    procedure :: stress_dev_2D_transient => unsupported_stress_2d
  end type
contains
  pure function sls_description(self) result(str)
    class(Viscoelastic_sls), intent(in) :: self
    character(len=275) :: str
    write(str,'(A,ES14.6,A,ES14.6)') 'Viscoelastic_sls (3D): ratio=',self%ratio,', tau=',self%tau
  end function
  pure real(real64) function sls_factor(dt,tau) result(b)
    real(real64), intent(in) :: dt,tau
    real(real64) :: x
    if(dt<=0 .or. tau<=0) error stop 'SLS requires dt,tau > 0'
    x=dt/tau
    if(x<1e-5_real64) then
      b=1-x/2+x*x/6-x*x*x/24+x**4/120
    else
      b=(1-exp(-x))/x
    end if
  end function
  pure function sls_stress(self,strain,dstrain,dtime,internal) result(res)
    class(Viscoelastic_sls), intent(in) :: self
    class(ten_3D2Osym), intent(in) :: strain,dstrain
    real(real64), intent(in) :: dtime
    class(Base_internal_elasticity_variables), intent(in) :: internal
    type(ten_3D2Osym) :: res,ds
    if(self%ratio<0) error stop 'SLS requires ratio >= 0'
    res=self%stress_3D(strain)
    ds=self%stress_3D(dstrain)
    select type(internal)
    type is(Sls_internal)
      res%vals=res%vals+exp(-dtime/self%tau)*internal%q+self%ratio*sls_factor(dtime,self%tau)*ds%vals
    class default
      error stop 'SLS requires Sls_internal'
    end select
  end function
  pure function sls_tangent(self,strain,dstrain,dtime,internal) result(res)
    class(Viscoelastic_sls), intent(in) :: self
    class(ten_3D2Osym), intent(in) :: strain,dstrain
    real(real64), intent(in) :: dtime
    class(Base_internal_elasticity_variables), intent(in) :: internal
    type(ten_3D4O3sym) :: res
    res=self%dstress_dstrain_3D(strain)
    res%vals=(1+self%ratio*sls_factor(dtime,self%tau))*res%vals
  end function
  pure subroutine sls_update(self,stress,strain,dstrain,dtime,internal)
    class(Viscoelastic_sls), intent(in) :: self
    type(ten_3D2Osym), intent(out) :: stress
    class(ten_3D2Osym), intent(in) :: strain,dstrain
    real(real64), intent(in) :: dtime
    class(Base_internal_elasticity_variables), intent(inout) :: internal
    type(ten_3D2Osym) :: eq
    stress=self%stress(strain,dstrain,dtime,internal)
    eq=self%stress_3D(strain)
    select type(internal)
    type is(Sls_internal)
      internal%q=stress%vals-eq%vals
    class default
      error stop 'SLS requires Sls_internal'
    end select
  end subroutine
  pure function sls_dev(self,strain,dstrain,dtime,internal) result(res)
    class(Viscoelastic_sls), intent(in) :: self
    class(ten_3D2Osym), intent(in) :: strain,dstrain
    real(real64), intent(in) :: dtime
    class(Base_internal_elasticity_variables), intent(in) :: internal
    type(ten_3D2Osym) :: res
    res=.dev. self%stress(strain,dstrain,dtime,internal)
  end function
  pure function sls_dev_tangent(self,strain,dstrain,dtime,internal) result(res)
    class(Viscoelastic_sls), intent(in) :: self
    class(ten_3D2Osym), intent(in) :: strain,dstrain
    real(real64), intent(in) :: dtime
    class(Base_internal_elasticity_variables), intent(in) :: internal
    type(ten_3D4O3sym) :: res
    type(ten_3D4O2sym) :: tmp
    type(iden_4O3T) :: i3
    type(iden_4O4T) :: i4
    tmp=(i4-(1.0_real64/3)*i3) .ddot. self%dstress_dstrain(strain,dstrain,dtime,internal)
    res=tmp%convert_3sym()
  end function
  pure function unsupported_stress_2d(self,strain,dstrain,dtime,internal) result(res)
    class(Viscoelastic_sls), intent(in) :: self
    class(ten_2D2Osym), intent(in) :: strain,dstrain
    real(real64), intent(in) :: dtime
    class(Base_internal_elasticity_variables), intent(in) :: internal
    type(ten_2D2Osym) :: res
    error stop 'SLS: use 3D material and external plane-stress condensation'
  end function
  pure function unsupported_tangent_2d(self,strain,dstrain,dtime,internal) result(res)
    class(Viscoelastic_sls), intent(in) :: self
    class(ten_2D2Osym), intent(in) :: strain,dstrain
    real(real64), intent(in) :: dtime
    class(Base_internal_elasticity_variables), intent(in) :: internal
    type(ten_2D4O3sym) :: res
    error stop 'SLS: use 3D material and external plane-stress condensation'
  end function
end module
