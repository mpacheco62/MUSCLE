! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_cpb06ex2
    !! CPB06ex2 criterion, with two transformations of the stress deviator (Plunkett et al., 2008, Eq. 8).
    !! Library (c44, c55, c66) on (xy, yz, xz) = paper (C66, C44, C55): xy from Eq. (13),
    !! yz and xz from the out-of-plane values of Tables 9 and 11.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_base
    implicit none
    private

    public :: CPB06ex2
    type, extends(Base_yield_critera) :: CPB06ex2
        real(real64) :: L(3,3)
            !! Normal block of C P_dev, which maps the normal stresses to those of Sigma (Eq. 9)
        real(real64) :: Lp(3,3)
            !! Normal block of Cp P_dev, for Sigmap (Eq. 9; suffix p = primed in the paper)
        real(real64) :: c_shear(3)
            !! Shear factors of C on the tensorial shears (xy, yz, xz)
        real(real64) :: cp_shear(3)
            !! Shear factors of Cp on the tensorial shears (xy, yz, xz), for Sigmap
        real(real64) :: k
            !! Strength-differential parameter k of Sigma (Eq. 8)
        real(real64) :: kp
            !! Strength-differential parameter kp of Sigmap (Eq. 8)
        real(real64) :: a
            !! Degree of homogeneity a (Eq. 8)
        real(real64) :: B
            !! Normalization: f = s in uniaxial tension s along RD (Eq. 11 with theta = 0)
    contains
        procedure :: init
        procedure :: stress_eq
        procedure :: dstressEq_dstress => dstressEq_dstress_cpb06ex2
        procedure :: ddstressEq_ddstress => ddstressEq_ddstress_cpb06ex2
    end type CPB06ex2

contains

    subroutine init(self, c11, c12, c13, c21, c22, c23, c31, c32, c33, c44, c55, c66, k, &
                    cp11, cp12, cp13, cp21, cp22, cp23, cp31, cp32, cp33, cp44, cp55, cp66, kp, a)
        !! Stores the maps of Plunkett et al. (2008), Eqs. (9)-(10). Normalizes Eq. (8)
        !! to unit uniaxial tension in RD, as in Eq. (11). Called when setting the material.
        implicit none
        class(CPB06ex2), intent(inout) :: self
        real(real64), intent(in) :: c11, c12, c13, c21, c22, c23, c31, c32, c33, c44, c55, c66, k
        real(real64), intent(in) :: cp11, cp12, cp13, cp21, cp22, cp23, cp31, cp32, cp33, &
                                    cp44, cp55, cp66, kp, a
        type(ten_3D2Osym) :: rd

        self%L = deviatoric_block(c11, c12, c13, c21, c22, c23, c31, c32, c33)
        self%Lp = deviatoric_block(cp11, cp12, cp13, cp21, cp22, cp23, cp31, cp32, cp33)
        self%c_shear = [c44, c55, c66]
        self%cp_shear = [cp44, cp55, cp66]
        self%k = k
        self%kp = kp
        self%a = a
        call rd%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
        self%B = 1D0
        self%B = 1D0/self%stress_eq(rd)
    end subroutine init

    pure function stress_eq(self, stress) result(res)
        !! Equivalent stress B F**(1/a), Plunkett et al. (2008), Eq. (8).
        !! Called by the return mapping; principal values come from eigenvals.
        use muscle_math_operations, only : eigenvals
        implicit none
        class(CPB06ex2), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        real(real64) :: res
        real(real64) :: lam(3), lamp(3)

        lam = eigenvals(linear_transform(self%L, self%c_shear, stress))
        lamp = eigenvals(linear_transform(self%Lp, self%cp_shear, stress))
        res = self%B*((abs(lam(1)) - self%k*lam(1))**self%a &
                      + (abs(lam(2)) - self%k*lam(2))**self%a &
                      + (abs(lam(3)) - self%k*lam(3))**self%a &
                      + (abs(lamp(1)) - self%kp*lamp(1))**self%a &
                      + (abs(lamp(2)) - self%kp*lamp(2))**self%a &
                      + (abs(lamp(3)) - self%kp*lamp(3))**self%a)**(1D0/self%a)
    end function stress_eq

    pure function dstressEq_dstress_cpb06ex2(self, stress) result(res)
        !! Analytical gradient of Plunkett et al. (2008), Eqs. (8)-(9).
        !! Called by the return mapping through the binding; requires nonzero f and a >= 2.
        use muscle_math_spectral_derivs, only : spectral_decomposition, spectral_gradient
        implicit none
        class(CPB06ex2), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res
        real(real64) :: lam(3), lamp(3), V(3,3), Vp(3,3), df(3), dfp(3)

        call spectral_decomposition(linear_transform(self%L, self%c_shear, stress), lam, V)
        call spectral_decomposition(linear_transform(self%Lp, self%cp_shear, stress), lamp, Vp)
        call principal_derivatives(self, lam, lamp, df, dfp)
        res = pull_back(self%L, self%c_shear, spectral_gradient(V, df)) &
              + pull_back(self%Lp, self%cp_shear, spectral_gradient(Vp, dfp))
    end function dstressEq_dstress_cpb06ex2

    pure function ddstressEq_ddstress_cpb06ex2(self, stress) result(res)
        !! Hessian of Plunkett et al. (2008), Eqs. (8)-(9), called by the return mapping.
        !! The mixed principal block is (1-a)/f times the outer product of the two gradients.
        use muscle_math_spectral_derivs, only : spectral_decomposition, spectral_gradient, &
                                                spectral_hessian
        implicit none
        class(CPB06ex2), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D4O3sym) :: res
        type(ten_3D2Osym) :: g, gp
        real(real64) :: lam(3), lamp(3), V(3,3), Vp(3,3), df(3), dfp(3), d2f(3,3), d2fp(3,3), cross

        call spectral_decomposition(linear_transform(self%L, self%c_shear, stress), lam, V)
        call spectral_decomposition(linear_transform(self%Lp, self%cp_shear, stress), lamp, Vp)
        call principal_derivatives(self, lam, lamp, df, dfp, d2f, d2fp, cross)
        g = pull_back(self%L, self%c_shear, spectral_gradient(V, df))
        gp = pull_back(self%Lp, self%cp_shear, spectral_gradient(Vp, dfp))
        res = pull_back_hessian(self%L, self%c_shear, spectral_hessian(lam, V, df, d2f)) &
              + pull_back_hessian(self%Lp, self%cp_shear, spectral_hessian(lamp, Vp, dfp, d2fp)) &
              + (2D0*cross)*(g .tdotsym. gp)
    end function ddstressEq_ddstress_cpb06ex2

    pure subroutine principal_derivatives(self, lam, lamp, df, dfp, d2f, d2fp, cross)
        ! Principal derivatives of f = B Psi**(1/a), Plunkett et al. (2008), Eq. (8):
        ! df/dlam_i = (f/Psi) psi_i**(a-1) (sgn(lam_i)-k), Psi = Phi + Phip adds Sigma and Sigmap.
        ! Called by both analytical derivatives. Each real power is evaluated once.
        implicit none
        class(CPB06ex2), intent(in) :: self
        real(real64), intent(in) :: lam(3), lamp(3)
        real(real64), intent(out) :: df(3), dfp(3)
        real(real64), intent(out), optional :: d2f(3,3), d2fp(3,3), cross
            ! The two diagonal blocks and the coefficient of the mixed block, supplied together.
        real(real64) :: psi(3), psip(3), p(3), pp(3), s(3), sp(3), phi, f, f_phi, q

        call principal_terms(lam, self%k, self%a, psi, p, s)
        call principal_terms(lamp, self%kp, self%a, psip, pp, sp)
        phi = psi(1)*p(1) + psi(2)*p(2) + psi(3)*p(3) &
              + psip(1)*pp(1) + psip(2)*pp(2) + psip(3)*pp(3)
        f = self%B*phi**(1D0/self%a)
        f_phi = f/phi
        df(1) = f_phi*p(1)*s(1)
        df(2) = f_phi*p(2)*s(2)
        df(3) = f_phi*p(3)*s(3)
        dfp(1) = f_phi*pp(1)*sp(1)
        dfp(2) = f_phi*pp(2)*sp(2)
        dfp(3) = f_phi*pp(3)*sp(3)
        if (.not. present(d2f)) return
        q = (self%a - 1D0)/phi
        cross = (1D0 - self%a)/f
        d2f = principal_hessian(lam, self%k, self%a, psi, p, s, df, f_phi, q)
        d2fp = principal_hessian(lamp, self%kp, self%a, psip, pp, sp, dfp, f_phi, q)
    end subroutine principal_derivatives

    pure subroutine principal_terms(lam, k, a, psi, p, s)
        ! psi_i = |lam_i|-k lam_i and p_i = psi_i**(a-1), Plunkett et al. (2008), Eq. (8).
        ! Called by principal_derivatives for either transformation without repeating its work.
        implicit none
        real(real64), intent(in) :: lam(3), k, a
        real(real64), intent(out) :: psi(3), p(3), s(3)

        s(1) = sign(1D0, lam(1)) - k
        s(2) = sign(1D0, lam(2)) - k
        s(3) = sign(1D0, lam(3)) - k
        psi(1) = lam(1)*s(1)
        psi(2) = lam(2)*s(2)
        psi(3) = lam(3)*s(3)
        p(1) = psi(1)**(a - 1D0)
        p(2) = psi(2)**(a - 1D0)
        p(3) = psi(3)**(a - 1D0)
    end subroutine principal_terms

    pure function principal_hessian(lam, k, a, psi, p, s, df, f_phi, q) result(res)
        ! One diagonal block of d2f/dlam2, Plunkett et al. (2008), Eq. (8), by the chain rule.
        ! Called twice by principal_derivatives, using the total Psi in f_phi and q.
        implicit none
        real(real64), intent(in) :: lam(3), k, a, psi(3), p(3), s(3), df(3), f_phi, q
        real(real64) :: res(3,3)

        res(1,1) = (a - 1D0)*f_phi*curvature(lam(1), k, a, psi(1), p(1), s(1)) - q*df(1)*p(1)*s(1)
        res(2,2) = (a - 1D0)*f_phi*curvature(lam(2), k, a, psi(2), p(2), s(2)) - q*df(2)*p(2)*s(2)
        res(3,3) = (a - 1D0)*f_phi*curvature(lam(3), k, a, psi(3), p(3), s(3)) - q*df(3)*p(3)*s(3)
        res(1,2) = -q*df(1)*p(2)*s(2)
        res(1,3) = -q*df(1)*p(3)*s(3)
        res(2,3) = -q*df(2)*p(3)*s(3)
        res(2,1) = res(1,2)
        res(3,1) = res(1,3)
        res(3,2) = res(2,3)
    end function principal_hessian

    pure function curvature(lam, k, a, psi, p, s) result(res)
        ! psi**(a-2)*(sgn(lam)-k)**2, Plunkett et al. (2008), Eq. (8).
        ! Called by principal_hessian, reusing p = psi**(a-1) instead of another power.
        implicit none
        real(real64), intent(in) :: lam, k, a, psi, p, s
        real(real64) :: res

        if (psi > 0D0) then
            res = p*s**2/psi
        else if (abs(lam) > 0D0) then
            res = 0D0  ! |k|=1: the contribution vanishes on this side of zero.
        else
            ! At lam=0, a=2 uses the mean of the lateral curvatures; zero for a>2.
            ! For k/=0, a=2 this is a convention at a point without a classical Hessian.
            res = merge(1D0 + k**2, 0D0, a <= 2D0)
        end if
    end function curvature

    pure function deviatoric_block(c11, c12, c13, c21, c22, c23, c31, c32, c33) result(L)
        ! Normal block L = C P_dev (Plunkett et al., 2008, Eq. 9). Called by init.
        ! The third column comes from the zero row sum, also in floating point.
        implicit none
        real(real64), intent(in) :: c11, c12, c13, c21, c22, c23, c31, c32, c33
        real(real64) :: L(3,3), m1, m2, m3

        m1 = (c11 + c12 + c13)/3D0
        m2 = (c21 + c22 + c23)/3D0
        m3 = (c31 + c32 + c33)/3D0
        L(1,1) = c11 - m1
        L(1,2) = c12 - m1
        L(1,3) = -L(1,1) - L(1,2)
        L(2,1) = c21 - m2
        L(2,2) = c22 - m2
        L(2,3) = -L(2,1) - L(2,2)
        L(3,1) = c31 - m3
        L(3,2) = c32 - m3
        L(3,3) = -L(3,1) - L(3,2)
    end function deviatoric_block

    pure function linear_transform(L, c_shear, stress) result(res)
        ! Sigma = L : stress (Plunkett et al., 2008, Eq. 9), for either transformation.
        ! Called by stress_eq and the derivatives. Zero row sums cancel pressure before products.
        implicit none
        real(real64), intent(in) :: L(3,3), c_shear(3)
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res
        real(real64) :: dxx, dyy

        dxx = stress%vals(1) - stress%vals(3)
        dyy = stress%vals(2) - stress%vals(3)
        res%vals(1) = L(1,1)*dxx + L(1,2)*dyy
        res%vals(2) = L(2,1)*dxx + L(2,2)*dyy
        res%vals(3) = L(3,1)*dxx + L(3,2)*dyy
        res%vals(4) = c_shear(1)*stress%vals(4)
        res%vals(5) = c_shear(2)*stress%vals(5)
        res%vals(6) = c_shear(3)*stress%vals(6)
    end function linear_transform

    pure function pull_back(L, c_shear, x) result(res)
        ! Adjoint L^T : x through Sigma = L : stress (Plunkett et al., 2008, Eq. 9).
        ! Called by the analytical derivatives; tensorial shear needs no extra factor two.
        implicit none
        real(real64), intent(in) :: L(3,3), c_shear(3)
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res
        real(real64) :: n(3)

        n = normal_pull_back(L, x%vals(1), x%vals(2), x%vals(3))
        call res%init(xx=n(1), yy=n(2), zz=n(3), xy=c_shear(1)*x%vals(4), &
                      yz=c_shear(2)*x%vals(5), xz=c_shear(3)*x%vals(6))
    end function pull_back

    pure function pull_back_hessian(L, c_shear, H) result(res)
        ! L^T : H : L through Sigma = L : stress (Plunkett et al., 2008, Eq. 9).
        ! Called by the analytical Hessian, for either transformation, in tensorial Voigt order.
        implicit none
        real(real64), intent(in) :: L(3,3), c_shear(3)
        type(ten_3D4O3sym), intent(in) :: H
        type(ten_3D4O3sym) :: res
        real(real64) :: cxx(3), cyy(3), czz(3), rxx(3), ryy(3), rzz(3), nxy(3), nyz(3), nxz(3)

        cxx = normal_pull_back(L, H%vals(1), H%vals(7), H%vals(12))
        cyy = normal_pull_back(L, H%vals(7), H%vals(2), H%vals(8))
        czz = normal_pull_back(L, H%vals(12), H%vals(8), H%vals(3))
        rxx = normal_pull_back(L, cxx(1), cyy(1), czz(1))
        ryy = normal_pull_back(L, cxx(2), cyy(2), czz(2))
        rzz = normal_pull_back(L, cxx(3), cyy(3), czz(3))
        nxy = c_shear(1)*normal_pull_back(L, H%vals(16), H%vals(13), H%vals(9))
        nyz = c_shear(2)*normal_pull_back(L, H%vals(19), H%vals(17), H%vals(14))
        nxz = c_shear(3)*normal_pull_back(L, H%vals(21), H%vals(20), H%vals(18))
        call res%init(xxxx=rxx(1), yyyy=ryy(2), zzzz=rzz(3), &
                      xyxy=c_shear(1)**2*H%vals(4), yzyz=c_shear(2)**2*H%vals(5), &
                      xzxz=c_shear(3)**2*H%vals(6), &
                      xxyy=rxx(2), yyzz=ryy(3), zzxy=nxy(3), &
                      xyyz=c_shear(1)*c_shear(2)*H%vals(10), &
                      yzxz=c_shear(2)*c_shear(3)*H%vals(11), &
                      xxzz=rxx(3), yyxy=nxy(2), zzyz=nyz(3), &
                      xyxz=c_shear(1)*c_shear(3)*H%vals(15), &
                      xxxy=nxy(1), yyyz=nyz(2), zzxz=nxz(3), &
                      xxyz=nyz(1), yyxz=nxz(2), xxxz=nxz(1))
    end function pull_back_hessian

    pure function normal_pull_back(L, x, y, z) result(res)
        ! Normal part of L^T : x (Plunkett et al., 2008, Eq. 9), used by both pull backs.
        implicit none
        real(real64), intent(in) :: L(3,3), x, y, z
        real(real64) :: res(3)

        res(1) = L(1,1)*x + L(2,1)*y + L(3,1)*z
        res(2) = L(1,2)*x + L(2,2)*y + L(3,2)*z
        res(3) = L(1,3)*x + L(2,3)*y + L(3,3)*z
    end function normal_pull_back

end module muscle_yield_cpb06ex2
