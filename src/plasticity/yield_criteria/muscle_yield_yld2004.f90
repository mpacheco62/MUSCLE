! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_yield_yld2004
    !! Module muscle_yield_yld2004
    !! ===========================
    !! Yld2004-18p anisotropic yield criterion (Barlat et al., 2005).
    !!
    !! Two linear transformations of the deviator s, S' = C':s and S'' = C'':s, with
    !! principal values lam'_i and lam''_j, define
    !! \[ \phi = \sum_i \sum_j |\lambda'_i - \lambda''_j|^a = 4 \bar{\sigma}^a \]
    !! (Barlat et al., 2005, Eq. 14) and the equivalent stress \( \bar{\sigma} = (\phi/4)^{1/a} \).
    !! C' and C'' take the 9 coefficients each of Barlat et al. (2005, Eq. 15): c12, c13, c21, c23,
    !! c31, c32 on the normal components and c44 (xy), c55 (yz), c66 (xz) on the shears, in the
    !! library order; the paper orders the shears yz, zx, xy (Eq. A.15), so
    !! (c44, c55, c66) here = (c66, c44, c55) in the paper.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_base
    implicit none
    private

    public :: Yld2004
    type, extends(Base_yield_critera) :: Yld2004
        !! Yld2004-18p yield criterion
        real(real64) :: Lp(3,3)
            !! Normal block of L' = C' P_dev, which maps (sxx, syy, szz) to (S'xx, S'yy, S'zz)
        real(real64) :: cp_shear(3)
            !! Shear factors of C' on (xy, yz, xz): (c'44, c'55, c'66)
        real(real64) :: Lpp(3,3)
            !! Normal block of L'' = C'' P_dev
        real(real64) :: cpp_shear(3)
            !! Shear factors of C'' on (xy, yz, xz): (c''44, c''55, c''66)
        real(real64) :: a
            !! Exponent a >= 2 for the analytical derivative bindings
    contains
        procedure :: init
        procedure :: stress_eq
        procedure :: dstressEq_dstress => dstressEq_dstress_yld2004
        procedure :: ddstressEq_ddstress => ddstressEq_ddstress_yld2004
    end type Yld2004

contains

    subroutine init(self, cp12, cp13, cp21, cp23, cp31, cp32, cp44, cp55, cp66, &
                    cpp12, cpp13, cpp21, cpp23, cpp31, cpp32, cpp44, cpp55, cpp66, a)
        !! Stores L' = C' P_dev and L'' = C'' P_dev (Barlat et al., 2005, Eqs. 11 and 15), so
        !! that each call applies them to the stress without forming the deviator. Analytical
        !! derivatives require a >= 2 and a nonzero equivalent stress. cpIJ = c'_IJ
        !! and cppIJ = c''_IJ; the shears follow the library order: 44 -> xy, 55 -> yz, 66 -> xz.
        implicit none
        class(Yld2004), intent(inout) :: self
        real(real64), intent(in) :: cp12, cp13, cp21, cp23, cp31, cp32, cp44, cp55, cp66
        real(real64), intent(in) :: cpp12, cpp13, cpp21, cpp23, cpp31, cpp32, cpp44, cpp55, cpp66
        real(real64), intent(in) :: a

        self%Lp = deviatoric_block(cp12, cp13, cp21, cp23, cp31, cp32)
        self%cp_shear(1) = cp44
        self%cp_shear(2) = cp55
        self%cp_shear(3) = cp66
        self%Lpp = deviatoric_block(cpp12, cpp13, cpp21, cpp23, cpp31, cpp32)
        self%cpp_shear(1) = cpp44
        self%cpp_shear(2) = cpp55
        self%cpp_shear(3) = cpp66
        self%a = a
    end subroutine init

    pure function stress_eq(self, stress) result(res)
        !! Equivalent stress (phi/4)**(1/a) (Barlat et al., 2005, Eq. 14), with the
        !! eigenvalues of S' and S'' from eigenvals.
        use muscle_math_operations, only : eigenvals
        implicit none
        class(Yld2004), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        real(real64) :: res

        real(real64) :: lp(3), lpp(3)

        lp = eigenvals(linear_transform(self%Lp, self%cp_shear, stress))     ! lam'
        lpp = eigenvals(linear_transform(self%Lpp, self%cpp_shear, stress))  ! lam''
        res = (0.25D0*(abs(lp(1) - lpp(1))**self%a + abs(lp(1) - lpp(2))**self%a &
                       + abs(lp(1) - lpp(3))**self%a                             &
                       + abs(lp(2) - lpp(1))**self%a + abs(lp(2) - lpp(2))**self%a &
                       + abs(lp(2) - lpp(3))**self%a                             &
                       + abs(lp(3) - lpp(1))**self%a + abs(lp(3) - lpp(2))**self%a &
                       + abs(lp(3) - lpp(3))**self%a))**(1D0/self%a)
    end function stress_eq

    pure function dstressEq_dstress_yld2004(self, stress) result(res)
        !! Analytical gradient df/dstress = L'^T : G' + L''^T : G'' (Barlat et al., 2005, Eqs. 11
        !! and 14), with G' = sum_i df/dlam'_i v'_i (x) v'_i from spectral_gradient, and the same
        !! for G''. Called by the return mapping through the dstressEq_dstress binding.
        use muscle_math_spectral_derivs, only : spectral_decomposition, spectral_gradient
        implicit none
        class(Yld2004), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res

        real(real64) :: lp(3), lpp(3), Vp(3,3), Vpp(3,3), dfp(3), dfpp(3)

        call spectral_decomposition(linear_transform(self%Lp, self%cp_shear, stress), lp, Vp)
        call spectral_decomposition(linear_transform(self%Lpp, self%cpp_shear, stress), lpp, Vpp)
        call principal_derivatives(self, lp, lpp, dfp, dfpp)
        res = pull_back(self%Lp, self%cp_shear, spectral_gradient(Vp, dfp)) &
              + pull_back(self%Lpp, self%cpp_shear, spectral_gradient(Vpp, dfpp))
    end function dstressEq_dstress_yld2004

    pure function ddstressEq_ddstress_yld2004(self, stress) result(res)
        !! Analytical Hessian Lp^T Hp Lp + Lpp^T Hpp Lpp + Lp^T M Lpp + (Lp^T M Lpp)^T
        !! (Barlat et al., 2005, Eqs. 11 and 14), where the suffixes p and pp mark the single- and
        !! double-primed quantities of the paper. Hp and Hpp come from spectral_hessian, each with
        !! the other tensor fixed. The cross term M = sum_ij d2f/dlp_i dlpp_j Ep_i (x) Epp_j,
        !! with Ep_i = vp_i (x) vp_i, has no eigenvector-rotation part: Ep_i does not depend on Spp.
        !! Called by the return mapping through the ddstressEq_ddstress binding.
        use muscle_math_spectral_derivs, only : spectral_decomposition, spectral_gradient, &
                                                spectral_hessian
        implicit none
        class(Yld2004), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D4O3sym) :: res

        real(real64) :: lp(3), lpp(3), Vp(3,3), Vpp(3,3), dfp(3), dfpp(3)
        real(real64) :: d2fp(3,3), d2fpp(3,3), d2fx(3,3)
        type(ten_3D2Osym) :: a1, a2, a3, b1, b2, b3

        call spectral_decomposition(linear_transform(self%Lp, self%cp_shear, stress), lp, Vp)
        call spectral_decomposition(linear_transform(self%Lpp, self%cpp_shear, stress), lpp, Vpp)
        call principal_derivatives(self, lp, lpp, dfp, dfpp, d2fp, d2fpp, d2fx)
        ! Cross term: Lp^T M Lpp + (Lp^T M Lpp)^T = sum_i 2 sym(A_i (x) B_i), with
        ! A_i = Lp^T : Ep_i and B_i = Lpp^T : sum_j d2f/dlp_i dlpp_j Epp_j
        a1 = pull_back(self%Lp, self%cp_shear, spectral_gradient(Vp, [1D0, 0D0, 0D0]))
        a2 = pull_back(self%Lp, self%cp_shear, spectral_gradient(Vp, [0D0, 1D0, 0D0]))
        a3 = pull_back(self%Lp, self%cp_shear, spectral_gradient(Vp, [0D0, 0D0, 1D0]))
        b1 = pull_back(self%Lpp, self%cpp_shear, spectral_gradient(Vpp, d2fx(1,:)))
        b2 = pull_back(self%Lpp, self%cpp_shear, spectral_gradient(Vpp, d2fx(2,:)))
        b3 = pull_back(self%Lpp, self%cpp_shear, spectral_gradient(Vpp, d2fx(3,:)))
        res = pull_back_hessian(self%Lp, self%cp_shear, spectral_hessian(lp, Vp, dfp, d2fp))     &
              + pull_back_hessian(self%Lpp, self%cpp_shear, spectral_hessian(lpp, Vpp, dfpp, d2fpp)) &
              + 2D0*((a1 .tdotsym. b1) + (a2 .tdotsym. b2) + (a3 .tdotsym. b3))
    end function ddstressEq_ddstress_yld2004

    pure subroutine principal_derivatives(self, lp, lpp, dfp, dfpp, d2fp, d2fpp, d2fx)
        ! Derivatives of f = (phi/4)**(1/a) with respect to lam'_i and lam''_j (Barlat et al.,
        ! 2005, Eq. 14). With d_ij = lam'_i - lam''_j and h(d) = |d|**a: df/dlam'_i =
        ! c sum_j h'(d_ij) and df/dlam''_j = -c sum_i h'(d_ij), c = f/(a phi).
        ! |d_ij|**(a-2) is the only real power per pair: h = p d**2, h' = a p d.
        ! Second derivatives, over x, y in (lp, lpp) (lam with one and two primes):
        ! d2f/dx dy = c d2phi/dx dy + (1 - a) (df/dx) (df/dy)/f, with d2phi/dlp_i dlpp_j = -d2h(d_ij),
        ! d2h = a (a-1) p the second derivative of h, and the sums of d2h on the diagonal blocks.
        ! Called by the analytical derivatives.
        implicit none
        class(Yld2004), intent(in) :: self
        real(real64), intent(in) :: lp(3)     ! lam'
        real(real64), intent(in) :: lpp(3)    ! lam''
        real(real64), intent(out) :: dfp(3)   ! df/dlam'_i
        real(real64), intent(out) :: dfpp(3)  ! df/dlam''_j
        real(real64), intent(out), optional :: d2fp(3,3)   ! d2f/dlp_i dlp_k
        real(real64), intent(out), optional :: d2fpp(3,3)  ! d2f/dlpp_j dlpp_l
        real(real64), intent(out), optional :: d2fx(3,3)   ! d2f/dlp_i dlpp_j

        real(real64) :: d(3,3), p(3,3), dh(3,3), phi, f, c, g, q

        ! d_ij = lam'_i - lam''_j
        d(1,1) = lp(1) - lpp(1)
        d(1,2) = lp(1) - lpp(2)
        d(1,3) = lp(1) - lpp(3)
        d(2,1) = lp(2) - lpp(1)
        d(2,2) = lp(2) - lpp(2)
        d(2,3) = lp(2) - lpp(3)
        d(3,1) = lp(3) - lpp(1)
        d(3,2) = lp(3) - lpp(2)
        d(3,3) = lp(3) - lpp(3)
        ! p_ij = |d_ij|**(a-2). Assumes a >= 2 (Barlat et al., 2005, use a = 6 or 8): with a < 2,
        ! p is infinite at d_ij = 0 and h' = a p d gives NaN instead of 0.
        p(1,1) = abs(d(1,1))**(self%a - 2D0)
        p(1,2) = abs(d(1,2))**(self%a - 2D0)
        p(1,3) = abs(d(1,3))**(self%a - 2D0)
        p(2,1) = abs(d(2,1))**(self%a - 2D0)
        p(2,2) = abs(d(2,2))**(self%a - 2D0)
        p(2,3) = abs(d(2,3))**(self%a - 2D0)
        p(3,1) = abs(d(3,1))**(self%a - 2D0)
        p(3,2) = abs(d(3,2))**(self%a - 2D0)
        p(3,3) = abs(d(3,3))**(self%a - 2D0)
        phi = p(1,1)*d(1,1)**2 + p(2,1)*d(2,1)**2 + p(3,1)*d(3,1)**2 &
              + p(1,2)*d(1,2)**2 + p(2,2)*d(2,2)**2 + p(3,2)*d(3,2)**2 &
              + p(1,3)*d(1,3)**2 + p(2,3)*d(2,3)**2 + p(3,3)*d(3,3)**2
        f = (0.25D0*phi)**(1D0/self%a)      ! f, as in stress_eq
        c = f/(self%a*phi)
        ! h'(d_ij) = a p_ij d_ij
        dh(1,1) = self%a*p(1,1)*d(1,1)
        dh(1,2) = self%a*p(1,2)*d(1,2)
        dh(1,3) = self%a*p(1,3)*d(1,3)
        dh(2,1) = self%a*p(2,1)*d(2,1)
        dh(2,2) = self%a*p(2,2)*d(2,2)
        dh(2,3) = self%a*p(2,3)*d(2,3)
        dh(3,1) = self%a*p(3,1)*d(3,1)
        dh(3,2) = self%a*p(3,2)*d(3,2)
        dh(3,3) = self%a*p(3,3)*d(3,3)
        dfp(1) = c*(dh(1,1) + dh(1,2) + dh(1,3))
        dfp(2) = c*(dh(2,1) + dh(2,2) + dh(2,3))
        dfp(3) = c*(dh(3,1) + dh(3,2) + dh(3,3))
        dfpp(1) = -c*(dh(1,1) + dh(2,1) + dh(3,1))
        dfpp(2) = -c*(dh(1,2) + dh(2,2) + dh(3,2))
        dfpp(3) = -c*(dh(1,3) + dh(2,3) + dh(3,3))
        if (.not. present(d2fp)) return

        ! c d2h(d_ij) = q p_ij, stored in p
        q = c*self%a*(self%a - 1D0)
        p(1,1) = q*p(1,1)
        p(1,2) = q*p(1,2)
        p(1,3) = q*p(1,3)
        p(2,1) = q*p(2,1)
        p(2,2) = q*p(2,2)
        p(2,3) = q*p(2,3)
        p(3,1) = q*p(3,1)
        p(3,2) = q*p(3,2)
        p(3,3) = q*p(3,3)
        g = (1D0 - self%a)/f
        d2fp(1,1) = g*dfp(1)*dfp(1) + p(1,1) + p(1,2) + p(1,3)
        d2fp(2,2) = g*dfp(2)*dfp(2) + p(2,1) + p(2,2) + p(2,3)
        d2fp(3,3) = g*dfp(3)*dfp(3) + p(3,1) + p(3,2) + p(3,3)
        d2fp(1,2) = g*dfp(1)*dfp(2)
        d2fp(1,3) = g*dfp(1)*dfp(3)
        d2fp(2,3) = g*dfp(2)*dfp(3)
        d2fp(2,1) = d2fp(1,2)
        d2fp(3,1) = d2fp(1,3)
        d2fp(3,2) = d2fp(2,3)
        d2fpp(1,1) = g*dfpp(1)*dfpp(1) + p(1,1) + p(2,1) + p(3,1)
        d2fpp(2,2) = g*dfpp(2)*dfpp(2) + p(1,2) + p(2,2) + p(3,2)
        d2fpp(3,3) = g*dfpp(3)*dfpp(3) + p(1,3) + p(2,3) + p(3,3)
        d2fpp(1,2) = g*dfpp(1)*dfpp(2)
        d2fpp(1,3) = g*dfpp(1)*dfpp(3)
        d2fpp(2,3) = g*dfpp(2)*dfpp(3)
        d2fpp(2,1) = d2fpp(1,2)
        d2fpp(3,1) = d2fpp(1,3)
        d2fpp(3,2) = d2fpp(2,3)
        d2fx(1,1) = g*dfp(1)*dfpp(1) - p(1,1)
        d2fx(1,2) = g*dfp(1)*dfpp(2) - p(1,2)
        d2fx(1,3) = g*dfp(1)*dfpp(3) - p(1,3)
        d2fx(2,1) = g*dfp(2)*dfpp(1) - p(2,1)
        d2fx(2,2) = g*dfp(2)*dfpp(2) - p(2,2)
        d2fx(2,3) = g*dfp(2)*dfpp(3) - p(2,3)
        d2fx(3,1) = g*dfp(3)*dfpp(1) - p(3,1)
        d2fx(3,2) = g*dfp(3)*dfpp(2) - p(3,2)
        d2fx(3,3) = g*dfp(3)*dfpp(3) - p(3,3)
    end subroutine principal_derivatives

    pure function deviatoric_block(c12, c13, c21, c23, c31, c32) result(L)
        ! Normal block of L = C P_dev, with C the normal block of C' or C'' (Barlat et al., 2005,
        ! Eq. 15: zero diagonal, -cIJ off it): L_ij = C_ij - (1/3) sum_k C_ik. Called by init.
        implicit none
        real(real64), intent(in) :: c12, c13, c21, c23, c31, c32
        real(real64) :: L(3,3)

        L(1,1) = (c12 + c13)/3D0
        L(1,2) = (c13 - 2D0*c12)/3D0
        L(1,3) = (c12 - 2D0*c13)/3D0
        L(2,1) = (c23 - 2D0*c21)/3D0
        L(2,2) = (c21 + c23)/3D0
        L(2,3) = (c21 - 2D0*c23)/3D0
        L(3,1) = (c32 - 2D0*c31)/3D0
        L(3,2) = (c31 - 2D0*c32)/3D0
        L(3,3) = (c31 + c32)/3D0
    end function deviatoric_block

    pure function linear_transform(L, c_shear, stress) result(res)
        ! S = L:stress with L = C P_dev (Barlat et al., 2005, Eq. 11): S' from (Lp, cp_shear),
        ! S'' from (Lpp, cpp_shear). The 3x3 block acts on the normal components and each factor
        ! on its tensorial shear. Called by stress_eq and the analytical derivatives.
        implicit none
        real(real64), intent(in) :: L(3,3)
        real(real64), intent(in) :: c_shear(3)
        type(ten_3D2Osym), intent(in) :: stress
        type(ten_3D2Osym) :: res

        res%vals(1) = L(1,1)*stress%vals(1) + L(1,2)*stress%vals(2) + L(1,3)*stress%vals(3)
        res%vals(2) = L(2,1)*stress%vals(1) + L(2,2)*stress%vals(2) + L(2,3)*stress%vals(3)
        res%vals(3) = L(3,1)*stress%vals(1) + L(3,2)*stress%vals(2) + L(3,3)*stress%vals(3)
        res%vals(4) = c_shear(1)*stress%vals(4)
        res%vals(5) = c_shear(2)*stress%vals(5)
        res%vals(6) = c_shear(3)*stress%vals(6)
    end function linear_transform

    pure function pull_back(L, c_shear, x) result(res)
        ! L^T : x, the adjoint of linear_transform: by the chain rule through S = L : stress
        ! (Barlat et al., 2005, Eq. 11) it takes a derivative with respect to S' or S'' to one
        ! with respect to the stress. Called by the analytical gradient and Hessian.
        implicit none
        real(real64), intent(in) :: L(3,3)
        real(real64), intent(in) :: c_shear(3)
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res

        real(real64) :: n(3)

        n = normal_pull_back(L, x%vals(1), x%vals(2), x%vals(3))
        call res%init(xx=n(1), yy=n(2), zz=n(3), xy=c_shear(1)*x%vals(4), &
                      yz=c_shear(2)*x%vals(5), xz=c_shear(3)*x%vals(6))
    end function pull_back

    pure function pull_back_hessian(L, c_shear, H) result(res)
        ! L^T : H : L, the second derivative by the chain rule through S = L : stress (Barlat
        ! et al., 2005, Eq. 11), for a fourth-order H with major and minor symmetries, by blocks
        ! of its 6x6 components: L^T Hnn L on the normal block, L^T Hns c on the normal-shear
        ! block and c Hss c on the shear block. Called by the analytical Hessian.
        implicit none
        real(real64), intent(in) :: L(3,3)
        real(real64), intent(in) :: c_shear(3)
        type(ten_3D4O3sym), intent(in) :: H
        type(ten_3D4O3sym) :: res

        real(real64), dimension(3) :: cxx, cyy, czz, rxx, ryy, rzz, nxy, nyz, nxz

        ! Normal block: L^T on the columns of Hnn, then on the rows of the result
        cxx = normal_pull_back(L, H%vals(1), H%vals(7), H%vals(12))
        cyy = normal_pull_back(L, H%vals(7), H%vals(2), H%vals(8))
        czz = normal_pull_back(L, H%vals(12), H%vals(8), H%vals(3))
        rxx = normal_pull_back(L, cxx(1), cyy(1), czz(1))
        ryy = normal_pull_back(L, cxx(2), cyy(2), czz(2))
        rzz = normal_pull_back(L, cxx(3), cyy(3), czz(3))
        ! Normal-shear block: L^T on each shear column of Hns, times its shear factor
        nxy = c_shear(1)*normal_pull_back(L, H%vals(16), H%vals(13), H%vals(9))
        nyz = c_shear(2)*normal_pull_back(L, H%vals(19), H%vals(17), H%vals(14))
        nxz = c_shear(3)*normal_pull_back(L, H%vals(21), H%vals(20), H%vals(18))
        ! Shear block, written in place
        call res%init(xxxx=rxx(1), yyyy=ryy(2), zzzz=rzz(3),                           &
                      xyxy=c_shear(1)**2*H%vals(4), yzyz=c_shear(2)**2*H%vals(5),      &
                      xzxz=c_shear(3)**2*H%vals(6),                                     &
                      xxyy=rxx(2), yyzz=ryy(3), zzxy=nxy(3),                           &
                      xyyz=c_shear(1)*c_shear(2)*H%vals(10),                            &
                      yzxz=c_shear(2)*c_shear(3)*H%vals(11),                            &
                      xxzz=rxx(3), yyxy=nxy(2), zzyz=nyz(3),                           &
                      xyxz=c_shear(1)*c_shear(3)*H%vals(15),                            &
                      xxxy=nxy(1), yyyz=nyz(2), zzxz=nxz(3),                           &
                      xxyz=nyz(1), yyxz=nxz(2),                                         &
                      xxxz=nxz(1))
    end function pull_back_hessian

    pure function normal_pull_back(L, x, y, z) result(res)
        ! L^T (x, y, z) on the normal components, from the chain rule through S = L : stress
        ! (Barlat et al., 2005, Eq. 11). Called by pull_back and, on the columns and rows of the
        ! 6x6 components, by pull_back_hessian.
        implicit none
        real(real64), intent(in) :: L(3,3)
        real(real64), intent(in) :: x, y, z
        real(real64) :: res(3)

        res(1) = L(1,1)*x + L(2,1)*y + L(3,1)*z
        res(2) = L(1,2)*x + L(2,2)*y + L(3,2)*z
        res(3) = L(1,3)*x + L(2,3)*y + L(3,3)*z
    end function normal_pull_back

end module muscle_yield_yld2004
