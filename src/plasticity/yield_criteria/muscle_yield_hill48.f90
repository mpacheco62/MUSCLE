! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_hill48
    use, intrinsic :: iso_fortran_env
    use muscle_yield_base
    implicit none
    PRIVATE

    PUBLIC :: Hill48
    type, extends(Base_yield_critera) :: Hill48
        real(real64) :: f,g,h,l,m,n
    contains
        procedure :: stress_eq
        ! procedure :: dstressEq_dstress => dstressEq_dstress_hill48
        ! procedure :: ddstressEq_ddstress => ddstressEq_ddstress_hill48
    end type Hill48

    contains
    pure function stress_eq(self, stress) result(res)
        use muscle_tensors
        implicit None
        class(Hill48), intent(in) :: self
        class(ten_3D2Osym), intent(in) :: stress
        real(real64) :: res

        res = (self%f*(stress%yy()-stress%zz())**2 + &
               self%g*(stress%zz()-stress%xx())**2 + &
               self%h*(stress%xx()-stress%yy())**2 + &
               2*self%l*stress%yz()**2 + &
               2*self%m*stress%xz()**2 + &
               2*self%n*stress%xy()**2   &
               )**0.5D0



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
  !     res = (3D0/(2D0*st_eq))*iden_4O4T() - (9D0/(4D0*st_eq_3))* (dev .tdot. dev) - (1D0/(2D0*st_eq))*iden_4O3T()
  ! end function ddstressEq_ddstress_vm
end module