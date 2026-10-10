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
      ! Normalization factor B (Cazacu et al., 2006, Eq. 12). It depends only on C1, k and a,
      ! so it is computed once in init instead of in every stress_eq call
      real(real64) :: B
      ! Normal block of L = C1 P_dev, the map sigma -> Sigma of Cazacu et al. (2006, Eq. 8) on
      ! the normal components (on the shears it is C2). Computed in init for the derivatives
      real(real64) :: L(3,3)

    contains
        procedure :: stress_eq
        procedure :: init
        procedure :: dstressEq_dstress => dstressEq_dstress_cpb06
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
        ! Principal values of C1 : dev(sigma) for a unit uniaxial tension along x;
        ! only used to compute B
        real(real64), dimension(3):: gamma

        self%C1 = reshape(  &
                          [c11, c21, c31, &
                           c12, c22, c32, &
                           c13, c23, c33], &
                          [3,3] )

        self%C2 = [c44, c55, c66]
        self%k = k
        self%a = a

        gamma(1) = (2D0*self%C1(1,1) - self%C1(1,2) - self%C1(1,3))/3D0
        gamma(2) = (2D0*self%C1(2,1) - self%C1(2,2) - self%C1(2,3))/3D0
        gamma(3) = (2D0*self%C1(3,1) - self%C1(3,2) - self%C1(3,3))/3D0

        self%B = sum((abs(gamma) - k*gamma)**a)**(-1D0/a)
        self%L = deviatoric_block(self%C1)

    end subroutine init


    pure function stress_eq(self, stress) result(res)
        use muscle_tensors
        use muscle_math_operations, only : eigenvals
        implicit None
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress

        real(real64) :: k, a

        type(ten_3D2Osym) :: sigma
        real(real64) :: res
        type(ten_3D2Osym) :: dev
        real(real64), dimension(3) :: sigma_eig

        real(real64) :: sxx, syy, szz, sxy, syz, sxz

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

        res = res*self%B

    end function stress_eq

    pure function dstressEq_dstress_cpb06(self, stress) result(res)
        !! Analytical gradient of the CPB06 equivalent stress (Cazacu et al., 2006, Eqs. 8, 9
        !! and 12): df/dsigma = L^T : sum_i (df/dlam_i) E_i, with lam_i and E_i the principal
        !! values and eigenprojections of Sigma = L : sigma = C : dev(sigma).
        !! Called by the return mapping through the dstressEq_dstress binding.
        use muscle_math_spectral_derivs, only : spectral_decomposition, spectral_gradient
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res

        real(real64) :: lam(3), V(3,3), df_dlam(3)

        call spectral_decomposition(linear_transform(self, stress), lam, V)
        call principal_derivatives(self, lam, df_dlam)
        res = pull_back(self, spectral_gradient(V, df_dlam))
    end function dstressEq_dstress_cpb06

    pure subroutine principal_derivatives(self, lam, df_dlam)
        ! Derivatives of f = B Phi**(1/a) (Cazacu et al., 2006, Eqs. 9 and 12) with respect to
        ! the principal values lam_i of Sigma, with Phi = sum_i psi_i**a, psi_i = |lam_i| - k lam_i:
        ! df/dlam_i = (f/Phi) g_i, with g_i = psi_i**(a-1) (sgn(lam_i) - k).
        ! Called by the analytical gradient.
        implicit none
        class(CPB06), intent(in) :: self
        real(real64), intent(in) :: lam(3)       ! Principal values lam_i of Sigma
        real(real64), intent(out) :: df_dlam(3)  ! df/dlam_i

        real(real64) :: s(3), psi(3), p(3), phi, f_phi

        ! s_i = sgn(lam_i) - k, psi_i = |lam_i| - k lam_i and p_i = psi_i**(a-1)
        s(1) = sign(1D0, lam(1)) - self%k
        s(2) = sign(1D0, lam(2)) - self%k
        s(3) = sign(1D0, lam(3)) - self%k
        psi(1) = lam(1)*s(1)
        psi(2) = lam(2)*s(2)
        psi(3) = lam(3)*s(3)
        p(1) = psi(1)**(self%a - 1D0)
        p(2) = psi(2)**(self%a - 1D0)
        p(3) = psi(3)**(self%a - 1D0)
        phi = psi(1)*p(1) + psi(2)*p(2) + psi(3)*p(3)
        f_phi = self%B*phi**(1D0/self%a)/phi  ! f/Phi, with f the equivalent stress of stress_eq
        df_dlam(1) = f_phi*p(1)*s(1)
        df_dlam(2) = f_phi*p(2)*s(2)
        df_dlam(3) = f_phi*p(3)*s(3)
    end subroutine principal_derivatives

    pure function deviatoric_block(C) result(L)
        ! Normal block of L = C P_dev, so that L : sigma = C : dev(sigma) (Cazacu et al., 2006,
        ! Eq. 8), with C the normal block C1 (general, not only symmetric) and P_dev the
        ! deviatoric projection: L_ij = C_ij - (1/3) sum_k C_ik. Each row of L sums to zero
        ! (L annihilates a hydrostatic stress). Called by init.
        implicit none
        real(real64), intent(in) :: C(3,3)
        real(real64) :: L(3,3)

        real(real64) :: m1, m2, m3

        ! m_i = (1/3) sum_k C_ik. The third column is taken from the row sum, so that each row
        ! sums to zero also in floating point
        m1 = (C(1,1) + C(1,2) + C(1,3))/3D0
        m2 = (C(2,1) + C(2,2) + C(2,3))/3D0
        m3 = (C(3,1) + C(3,2) + C(3,3))/3D0
        L(1,1) = C(1,1) - m1
        L(1,2) = C(1,2) - m1
        L(1,3) = -L(1,1) - L(1,2)
        L(2,1) = C(2,1) - m2
        L(2,2) = C(2,2) - m2
        L(2,3) = -L(2,1) - L(2,2)
        L(3,1) = C(3,1) - m3
        L(3,2) = C(3,2) - m3
        L(3,3) = -L(3,1) - L(3,2)
    end function deviatoric_block

    pure function linear_transform(self, stress) result(res)
        ! Sigma = L : stress = C : dev(stress) (Cazacu et al., 2006, Eq. 8): L acts on the
        ! normal components and C2 on the tensorial shears (xy, yz, xz). Called by the
        ! analytical gradient.
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res

        real(real64) :: dxx, dyy

        ! The rows of L sum to zero, so L acts on the differences of the normal stresses: the
        ! pressure cancels before the products and not by rounding after them
        dxx = stress%vals(1) - stress%vals(3)
        dyy = stress%vals(2) - stress%vals(3)
        res%vals(1) = self%L(1,1)*dxx + self%L(1,2)*dyy
        res%vals(2) = self%L(2,1)*dxx + self%L(2,2)*dyy
        res%vals(3) = self%L(3,1)*dxx + self%L(3,2)*dyy
        res%vals(4) = self%C2(1)*stress%vals(4)
        res%vals(5) = self%C2(2)*stress%vals(5)
        res%vals(6) = self%C2(3)*stress%vals(6)
    end function linear_transform

    pure function pull_back(self, x) result(res)
        ! L^T : x, the adjoint of linear_transform: by the chain rule through Sigma = L : stress
        ! (Cazacu et al., 2006, Eq. 8) it takes a derivative with respect to Sigma to one with
        ! respect to the stress. Called by the analytical gradient.
        implicit none
        class(CPB06), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res

        real(real64) :: n(3)

        n = normal_pull_back(self, x%vals(1), x%vals(2), x%vals(3))
        call res%init(xx=n(1), yy=n(2), zz=n(3), xy=self%C2(1)*x%vals(4), &
                      yz=self%C2(2)*x%vals(5), xz=self%C2(3)*x%vals(6))
    end function pull_back

    pure function normal_pull_back(self, x, y, z) result(res)
        ! L^T (x, y, z) on the normal components, from the chain rule through Sigma = L : stress
        ! (Cazacu et al., 2006, Eq. 8). Called by pull_back.
        implicit none
        class(CPB06), intent(in) :: self
        real(real64), intent(in) :: x, y, z
        real(real64) :: res(3)

        res(1) = self%L(1,1)*x + self%L(2,1)*y + self%L(3,1)*z
        res(2) = self%L(1,2)*x + self%L(2,2)*y + self%L(3,2)*z
        res(3) = self%L(1,3)*x + self%L(2,3)*y + self%L(3,3)*z
    end function normal_pull_back

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