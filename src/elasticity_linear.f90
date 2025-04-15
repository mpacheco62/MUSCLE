module mod_elasticity_linear
    use, intrinsic :: iso_fortran_env
    use tensors_types
    use mod_base_elasticity
    implicit none
    PRIVATE

    PUBLIC :: Elasticity_linear
    type, extends(Base_elasticity) :: Elasticity_linear
        real(real64), private :: young, poisson
        real(real64), private :: e1, e2, e3
        type(ten_3D4O3sym), private :: tan
        logical, private :: tan_init=.FALSE.
    contains
        procedure :: dstress_dstrain => dstress_dstrain_linear
        procedure :: stress => stress_linear
        procedure :: set_parameters
    end type Elasticity_linear

    contains
    pure subroutine set_parameters(self, young, poisson, calc_tan)
        implicit none
        logical, optional, intent(in) :: calc_tan
        class(Elasticity_linear), intent(inout) :: self
        real(real64), intent(in) :: young, poisson
        real(real64) :: e1, e2, e3
        logical :: ccalc_tan
        
        ccalc_tan = .True.
        if (present(calc_tan)) ccalc_tan = calc_tan

        self%young = young
        self%poisson = poisson
        e1 = young*(1D0-poisson)/((1D0+poisson)*(1D0-2D0*poisson))
        e2 = young*poisson/((1D0+poisson)*(1D0-2D0*poisson))
        e3 = 0.5D0*young/(1D0+poisson)
        self%e1 = e1
        self%e2 = e2
        self%e3 = e3

        self%tan_init = ccalc_tan
        call self%tan%init( xxxx=e1,  yyyy=e1,  zzzz=e1,  &
                            xxyy=e2,  yyzz=e2,  xxzz=e2,  &
                           xxxy=0D0, xxyz=0D0, xxxz=0D0,  &
                           yyxy=0D0, yyyz=0D0, yyxz=0D0,  &
                           zzxy=0D0, zzyz=0D0, zzxz=0D0,  &
                            xyxy=e3,  yzyz=e3,  xzxz=e3,  &
                           xyyz=0D0, yzxz=0D0,  xyxz=0D0  &
                           )

    end subroutine

    pure function stress_linear(self, strain) result(res)
        implicit none
        class(Elasticity_linear), intent(in) :: self
        class(ten_3D2Osym), intent(in) :: strain
        type(ten_3D2Osym) :: res
        real(real64) :: e1, e2, e3, xx, yy, zz, xy, yz, xz
        e1 = self%e1
        e2 = self%e2
        e3 = self%e3
        xx = e1*strain%xx() + e2*(strain%yy()+strain%zz())
        yy = e1*strain%yy() + e2*(strain%zz()+strain%xx())
        zz = e1*strain%zz() + e2*(strain%xx()+strain%yy())
        xy = 2*e3*strain%xy()
        yz = 2*e3*strain%yz()
        xz = 2*e3*strain%xz()
        call res%init(xx=xx, yy=yy, zz=zz, xy=xy, yz=yz, xz=xz)
        return 
    end function stress_linear

    pure function dstress_dstrain_linear(self, strain) result(res)
        implicit none
        class(Elasticity_linear), intent(in) :: self
        class(ten_3D2Osym), intent(in) :: strain
        type(ten_3D4O3sym) :: res
        real(real64) :: e1, e2, e3
        if (.not. self%tan_init) then
            e1 = self%e1
            e2 = self%e2
            e3 = self%e3
            call res%init( xxxx=e1,  yyyy=e1,  zzzz=e1,  &
                           xxyy=e2,  yyzz=e2,  xxzz=e2,  &
                          xxxy=0D0, xxyz=0D0, xxxz=0D0,  &
                          yyxy=0D0, yyyz=0D0, yyxz=0D0,  &
                          zzxy=0D0, zzyz=0D0, zzxz=0D0,  &
                           xyxy=e3,  yzyz=e3,  xzxz=e3,  &
                          xyyz=0D0, yzxz=0D0,  xyxz=0D0  &
                          )
            return
        else 
            res = self%tan
            return
        end if
    end function dstress_dstrain_linear

end module