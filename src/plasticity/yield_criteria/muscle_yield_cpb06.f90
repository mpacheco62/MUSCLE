! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_cpb06
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_base
    implicit none
    PRIVATE

    PUBLIC :: CPB06
    type, extends(Base_yield_critera) :: CPB06
      real(real64), dimension(3,3) :: C1
      real(real64), dimension(3) :: C2
      real(real64) :: k
      real(real64) :: a

    contains
        procedure :: stress_eq
        procedure :: init
        ! procedure :: dstressEq_dstress => dstressEq_dstress_vm
        ! procedure :: ddstressEq_ddstress => ddstressEq_ddstress_vm
    end type CPB06

    contains
    subroutine init(self,           &
                    c11, c12, c13,  &
                    c21, c22, c23,  &
                    c31, c32, c33,  &
                    c44, c55, c66,  &
                    k, a)
        implicit none
        class(CPB06), intent(inout) :: self
        real(real64), intent(in) :: c11, c12, c13, c21, c22, c23, c31, c32, c33, c44, c55, c66
        real(real64), intent(in) :: k, a

        self%C1 = reshape(  &
                          [c11, c21, c31, &
                           c12, c22, c32, &
                           c13, c23, c33], &
                          [3,3] )

        self%C2 = [c44, c55, c66]
        self%k = k
        self%a = a

    end subroutine init


    pure function stress_eq(self, stress) result(res)
        use muscle_tensors
        use muscle_math_operations, only : eigenvals
        implicit None
        class(CPB06), intent(in) :: self
        class(ten_3D2Osym), intent(in) :: stress

        real(real64) :: k, a

        type(ten_3D2Osym) :: sigma
        real(real64) :: res
        type(ten_3D2Osym) :: dev
        real(real64), dimension(3) :: sigma_eig

        real(real64) :: sxx, syy, szz, sxy, syz, sxz
        real(real64), dimension(3):: gamma
        real(real64) :: b

        k = self%k
        a = self%a

        dev = .dev. stress
        sxx = dev%xx()*self%C1(1,1) + dev%yy()*self%C1(1,2) + dev%zz()*self%C1(1,3)
        syy = dev%xx()*self%C1(2,1) + dev%yy()*self%C1(2,2) + dev%zz()*self%C1(2,3)
        szz = dev%xx()*self%C1(3,1) + dev%yy()*self%C1(3,2) + dev%zz()*self%C1(3,3)
        sxy = dev%xy()*self%C2(1)
        syz = dev%yz()*self%C2(2)
        sxz = dev%xz()*self%C2(3)

        call sigma%init(xx=sxx, yy=syy, zz=szz, xy=sxy, yz=syz, xz=sxz)
        sigma_eig = eigenvals(sigma)

        res = sum((abs(sigma_eig) - k*sigma_eig)**a)**(1D0/a)

        gamma(1) = (2D0*self%C1(1,1) - self%C1(1,2) - self%C1(1,3))/3D0
        gamma(2) = (2D0*self%C1(2,1) - self%C1(2,2) - self%C1(2,3))/3D0
        gamma(3) = (2D0*self%C1(3,1) - self%C1(3,2) - self%C1(3,3))/3D0

        b = sum((abs(gamma) - k*gamma)**a)**(-1D0/a)

        res = res*b

    end function stress_eq

  !   pure function dstressEq_dstress_vm(self, stress) result(res)
  !     ! use muscle_math_operations
  !     use muscle_tensors
  !     implicit None
  !     class(VonMises), intent(in) :: self
  !     class(ten_3D2Osym), intent(in) :: stress
  !     type(ten_3D2Osym) :: res
     
  !     type(ten_3D2Osym) :: dev
  !     ! real(real64) :: stress_eq

  !     dev = .dev. stress
  !     res = 1.5D0*dev/self%stress_eq(stress)
  !     ! res = 1.5D0*dev
  !     ! stress_eq = self%stress_eq(stress)
  !     ! res = res/stress_eq
  !     ! res = 1.5D0*dev/self%stress_eq(stress)
  !     ! res = 1.5D0*dev/self%stress_eq(stress)

  !   end function dstressEq_dstress_vm

  !   pure function ddstressEq_ddstress_vm(self, stress) result(res)
  !     use muscle_tensors
  !     implicit None
  !     class(VonMises), intent(in) :: self
  !     class(ten_3D2Osym), intent(in) :: stress
  !     type(ten_3D4O3sym) :: res
    
  !     type(ten_3D2Osym) :: dev
  !     real(real64) :: st_eq, st_eq_3

  !     dev = .dev. stress
  !     st_eq = self%stress_eq(stress)
  !     st_eq_3 = st_eq**3

  !     !
  !     !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
  !     !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
  !     !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
  !     !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
  !     !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
  !     !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |

  !     ! ! res = (3D0/(2D0*st_eq))*IDEN4_3sym() - (9D0/(4D0*st_eq_3))* (dev .tdot. dev) 
  !     ! res = (3D0/(2D0*st_eq))*iden_4O4T()
  !     ! res = (9D0/(4D0*st_eq_3))* (dev .tdot. dev)
  !     ! ! res = (1D0/(2D0*st_eq))*iden_4O3T()
  !     ! res = (3D0/(2D0*st_eq))*iden_4O4T() - (9D0/(4D0*st_eq_3))* (dev .tdot. dev)
  !     ! res = (9D0/(4D0*st_eq_3))* (dev .tdot. dev) - (1D0/(2D0*st_eq))*iden_4O3T()
  !     ! res = ((3D0/(2D0*st_eq))*iden_4O4T() - (9D0/(4D0*st_eq_3))* (dev .tdot. dev)) - (1D0/(2D0*st_eq))*iden_4O3T()
  !     res = (3D0/(2D0*st_eq))*iden_43DO4T() - (9D0/(4D0*st_eq_3))* (dev .tdot. dev) - (1D0/(2D0*st_eq))*iden_4O3T()
  ! end function ddstressEq_ddstress_vm
end module