! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_math_spectral_derivs
    !! # Module mod_muscle_math_spectral_derivs
    !!
    !! This module provides exact, analytical, eigenvector-free first and second derivatives 
    !! of eigenvalues with respect to their parent symmetric second-order tensor.
    !!
    !! ## Mathematical Formulation
    !!
    !! The eigenvalues \(\Sigma_a\) (\(a=1,2,3\)) are roots of the characteristic polynomial:
    !!
    !! \[ P(\Sigma_a) = \Sigma_a^3 - I_1 \Sigma_a^2 + I_2 \Sigma_a - I_3 = 0 \]
    !!
    !! By applying implicit differentiation and the Cayley-Hamilton theorem, the derivatives 
    !! are calculated without computing eigenvectors or performing matrix inversions, ensuring 
    !! absolute numerical stability even for singular tensors (e.g., when \(\det(\Sigma) = 0\)).
    !!
    !! See the technical documentation for further algebraic details.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    implicit none
    private

    public :: dEigenvalues_dTensor
    public :: d2Eigenvalues_dTensor2
    public :: spectral_decomposition
    public :: spectral_gradient
    public :: spectral_hessian

    real(real64), parameter :: EPS_TOL = 1.0D-11
    real(real64), parameter :: EPS_TOL_SECOND = 1.0D-6
    !! Safe tolerance threshold to prevent division by zero at repeated roots (singularities)
    real(real64), parameter :: REPEATED_TOL = 1.0D-9
    !! Relative eigenvalue gap below which theta_ab of spectral_hessian takes its limit value.
    !! The limit assumes phi is twice differentiable between the two eigenvalues: 1e-9 keeps
    !! that band narrow where it is not (CPB06, a = 2, k /= 0: h'' jumps at lam = 0) and the
    !! quotient error below 1e-7.

contains
    pure function dEigenvalues_dTensor(Sigma, eigenvalues) result(res)
        !!# First Derivative of Eigenvalues
        !!
        !! Computes the first analytical derivative of all three eigenvalues of a symmetric 
        !! second-order tensor \(\Sigma_{ij}\) with respect to itself:
        !!
        !! \[ H_{mn}^{(a)} = \frac{\partial \Sigma_a}{\partial \Sigma_{mn}} \]
        !!
        !! ## Mathematical Formulation
        !! The computation is based on the implicit differentiation of the characteristic 
        !! polynomial and the Cayley-Hamilton theorem. This approach yields a closed-form 
        !! expression free of eigenvectors:
        !!
        !! \[ H_{mn}^{(a)} = \frac{\Sigma^2_{mn} - (I_1 - \Sigma_a)\Sigma_{mn} + (\Sigma_a^2 - I_1\Sigma_a + I_2)\delta_{mn}}{3\Sigma_a^2 - 2I_1\Sigma_a + I_2} \]
        !!
        !! ## Singularity Handling (Repeated Roots)
        !! When eigenvalues coincide (e.g., uniaxial tension or hydrostatic states), the denominator 
        !! approaches zero, leading to an indeterminate form \(0/0\). This routine robustly handles 
        !! such states by exploiting the identity property \(\sum_{a=1}^3 \frac{\partial \Sigma_a}{\partial \pmb{\Sigma}} = \mathbf{I}\).
        !! 
        !! - **Two identical eigenvalues:** The remaining subspace is isotropically divided among the repeated roots.
        !! - **Three identical eigenvalues:** The identity tensor is isotropically divided by 3.
        !!
        !! @note This function is pure and returns an array of three symmetric second-order tensors.
        type(ten_3D2Osym), intent(in)  :: Sigma          !! Symmetric second-order tensor \(\Sigma_{ij}\)
        real(real64), intent(in)       :: eigenvalues(3) !! Evaluated eigenvalues \(\Sigma_1, \Sigma_2, \Sigma_3\)
        type(ten_3D2Osym)              :: res(3)         !! Resulting derivatives \(H_{ij}^{(a)}\) for \(a=1,2,3\)
        
        real(real64)      :: I1, I2, D_scal(3), lam
        type(ten_3D2Osym) :: F_tensor
        type(iden_2O)     :: I2O
        integer           :: a, num_sing
        integer           :: good_idx, sing_idx(3)

        ! 1. Compute Invariants
        ! First invariant: trace(\Sigma)
        I1 = Sigma%vals(1) + Sigma%vals(2) + Sigma%vals(3)
        ! Second invariant: 0.5*(I1^2 - tr(\Sigma^2))
        I2 = 0.5D0 * (I1**2 - (Sigma .ddot. Sigma))
        ! Second-order identity tensor (\delta_ij)

        num_sing = 0
        good_idx = 1 ! Default initialization

        ! 2. First Pass: Compute denominators and identify singular roots
        do a = 1, 3
            lam = eigenvalues(a)
            ! Characteristic polynomial derivative (scalar denominator)
            D_scal(a) = (3.0D0*lam - 2.0D0*I1)*lam + I2
            
            ! Check for distinct eigenvalues (well-conditioned denominator)
            if (abs(D_scal(a)) > EPS_TOL) then
                ! Numerator tensor using Cayley-Hamilton expression
                F_tensor = Sigma%square() - (I1 - lam)*Sigma + (lam*(lam - I1))*I2O + I2*I2O 
                res(a) = F_tensor / D_scal(a)
                good_idx = a
            else
                ! Mark as repeated/singular root for post-processing
                num_sing = num_sing + 1
                sing_idx(num_sing) = a
            end if
        end do

        ! 3. Second Pass: Repair singular roots using orthogonal projections
        if (num_sing == 2) then
            ! Two repeated eigenvalues (e.g., pure uniaxial tension/compression)
            ! The degenerate subspace is isotropically averaged from the remaining identity space.
            F_tensor = 0.5D0 * (I2O - res(good_idx))
            res(sing_idx(1)) = F_tensor
            res(sing_idx(2)) = F_tensor
            
        else if (num_sing == 3) then
            ! Three repeated eigenvalues (e.g., pure hydrostatic stress)
            ! All principal directions are isotropic.
            F_tensor = (1.0D0 / 3.0D0) * I2O
            res(1) = F_tensor
            res(2) = F_tensor
            res(3) = F_tensor
        end if

    end function dEigenvalues_dTensor



    pure function d2Eigenvalues_dTensor2(Sigma, eigenvalues, H_tensors) result(res)
        !!# Second Derivative (Hessian) of Eigenvalues
        !!
        !! Computes the second analytical derivative of all three eigenvalues of a symmetric 
        !! second-order tensor \(\Sigma_{ij}\) with respect to itself:
        !!
        !! \[ \mathcal{H}_{mnop}^{(a)} = \frac{\partial^2 \Sigma_a}{\partial \Sigma_{mn} \partial \Sigma_{op}} \]
        !!
        !! This function is pure and returns an array of three fully symmetric fourth-order tensors.
        implicit none
        type(ten_3D2Osym), intent(in)   :: Sigma          !! Symmetric second-order tensor \(\Sigma_{ij}\)
        real(real64), intent(in)        :: eigenvalues(3) !! Evaluated eigenvalues \(\Sigma_1, \Sigma_2, \Sigma_3\)
        type(ten_3D2Osym), intent(in)   :: H_tensors(3)   !! Pre-computed first derivatives \(H_{ij}^{(a)}\)
        type(ten_3D4O3sym)              :: res(3)         !! Resulting Hessians for \(a=1,2,3\)
        
        type(ten_3D4O3sym) :: dF_dSigma, M_tensor, test
        type(ten_3D4O3sym) :: N_tensor
        type(ten_3D2Osym) :: H, K
        type(iden_2O) :: I2O
        type(iden_4O4T) :: I4O4T
        type(iden_4O3T) :: I4O3T
        real(real64) :: I1, I2, D_scal, lam
        integer :: a, num_sing, good_idx, sing_idx(3)

        ! First invariant: trace(\Sigma)
        I1 = Sigma%vals(1) + Sigma%vals(2) + Sigma%vals(3)
        ! Second invariant
        I2 = 0.5D0 * (I1**2 - (Sigma .ddot. Sigma))
       
        ! Fourth-order tensor M_ijkl (analytical derivative of \Sigma^2) 
        M_tensor = Sigma%dsquare()

        num_sing = 0
        good_idx = 1

        do a = 1, 3
            lam = eigenvalues(a)
            D_scal = (3.0D0*lam - 2.0D0*I1)*lam + I2

            if (abs(D_scal) > EPS_TOL_SECOND) then
                H = H_tensors(a)                
                ! 2. Symmetrized fourth-order N_tensor (numerator) using Cayley-Hamilton-based expression
                N_tensor = M_tensor &
                           - (I1 - lam)*I4O4T &
                           + (I1 - lam)*I4O3T &
                           - 2.0D0 * (Sigma .tdotsym. I2O) &
                           + 2.0D0 * (Sigma .tdotsym. H) &
                           + 2.0D0 * (2.0D0*lam - I1) * (I2O .tdotsym. H) &
                           - 2.0D0 * (3.0D0*lam - I1) * (.tdotsym. H)
    
                res(a) = N_tensor / D_scal
                good_idx = a
            else 
                num_sing = num_sing + 1
                sing_idx(num_sing) = a
            end if
        end do

        ! Repair singular roots using orthogonal projections
        if (num_sing == 2) then
            ! Two repeated eigenvalues: Isotropic averaging of the degenerate subspace
            res(sing_idx(1)) = (-0.5D0) * res(good_idx)
            res(sing_idx(2)) = (-0.5D0) * res(good_idx)
        else if (num_sing == 3) then
            ! Three repeated eigenvalues: Isotropic distribution of the identity subspace
            res(1) = 0.0D0
            res(2) = 0.0D0
            res(3) = 0.0D0
        end if

    end function d2Eigenvalues_dTensor2

    pure subroutine spectral_decomposition(T, lam, V)
        !! Eigenvalues and orthonormal eigenvectors of T. The eigenvector of the eigenvalue
        !! farthest from the middle one is the normalized cross product of two rows of
        !! T - lam I (Kopp, 2008, Eq. 39); the other two come from one closed-form plane
        !! rotation (no iteration) that diagonalizes the 2x2 block of T in the plane normal to
        !! it, which stays accurate when they (nearly) coincide. That eigenvalue is at least
        !! half the eigenvalue range away from the other two, so the cross product needs no fallback.
        !! Called by the CPB06 and Yld2004-18p derivatives.
        use muscle_math_operations, only : eigenvals
        type(ten_3D2Osym), intent(in) :: T   !! Symmetric tensor
        real(real64), intent(out) :: lam(3)  !! Eigenvalues, in descending order
        real(real64), intent(out) :: V(3,3)  !! Eigenvectors: column a belongs to lam(a)

        real(real64) :: r1(3), r2(3), r3(3), c1(3), c2(3), c3(3), ca(3), u(3), w(3), tu(3), tw(3)
        real(real64) :: n1, n2, n3, na, p, q, r, theta, tg, cs, sn
        integer :: a, b

        lam = eigenvals(T)

        ! a: most separated eigenvalue; b, b+1: the remaining pair
        if (lam(1) - lam(2) >= lam(2) - lam(3)) then
            a = 1
            b = 2
        else
            a = 3
            b = 1
        end if

        ! Rows of T - lam(a) I and their pairwise cross products; the largest one is v_a
        r1 = [T%vals(1) - lam(a), T%vals(4), T%vals(6)]
        r2 = [T%vals(4), T%vals(2) - lam(a), T%vals(5)]
        r3 = [T%vals(6), T%vals(5), T%vals(3) - lam(a)]
        c1 = cross(r1, r2)
        c2 = cross(r2, r3)
        c3 = cross(r3, r1)
        n1 = c1(1)**2 + c1(2)**2 + c1(3)**2
        n2 = c2(1)**2 + c2(2)**2 + c2(3)**2
        n3 = c3(1)**2 + c3(2)**2 + c3(3)**2
        if (n1 >= n2 .and. n1 >= n3) then
            ca = c1
            na = n1
        else if (n2 >= n3) then
            ca = c2
            na = n2
        else
            ca = c3
            na = n3
        end if
        if (na <= 0.0D0) then
            ! T is a multiple of the identity: any orthonormal basis is an eigenbasis
            V = 0.0D0
            V(1,1) = 1.0D0
            V(2,2) = 1.0D0
            V(3,3) = 1.0D0
            return
        end if
        V(:,a) = ca/sqrt(na)

        ! Orthonormal basis (u, w) of the plane normal to V(:,a)
        if (abs(V(1,a)) > abs(V(2,a))) then
            u = [-V(3,a), 0.0D0, V(1,a)]/sqrt(V(1,a)**2 + V(3,a)**2)
        else
            u = [0.0D0, V(3,a), -V(2,a)]/sqrt(V(2,a)**2 + V(3,a)**2)
        end if
        w = cross(V(:,a), u)

        ! Plane rotation (cs, sn) that diagonalizes [[p, q], [q, r]], the restriction of T to (u, w)
        tu = [T%vals(1)*u(1) + T%vals(4)*u(2) + T%vals(6)*u(3), &
              T%vals(4)*u(1) + T%vals(2)*u(2) + T%vals(5)*u(3), &
              T%vals(6)*u(1) + T%vals(5)*u(2) + T%vals(3)*u(3)]
        tw = [T%vals(1)*w(1) + T%vals(4)*w(2) + T%vals(6)*w(3), &
              T%vals(4)*w(1) + T%vals(2)*w(2) + T%vals(5)*w(3), &
              T%vals(6)*w(1) + T%vals(5)*w(2) + T%vals(3)*w(3)]
        p = dot_product(u, tu)
        q = dot_product(u, tw)
        r = dot_product(w, tw)
        tg = 0.0D0
        if (abs(q) > 0.0D0) then
            theta = (r - p)/(2.0D0*q)
            tg = sign(1.0D0, theta)/(abs(theta) + sqrt(theta*theta + 1.0D0))
        end if
        cs = 1.0D0/sqrt(tg*tg + 1.0D0)
        sn = tg*cs
        ! Rotated eigenvalues: p - tg*q for cs*u - sn*w, and r + tg*q for sn*u + cs*w
        if (p - tg*q >= r + tg*q) then
            V(:,b) = cs*u - sn*w
            V(:,b+1) = sn*u + cs*w
        else
            V(:,b) = sn*u + cs*w
            V(:,b+1) = cs*u - sn*w
        end if
    end subroutine spectral_decomposition

    pure function spectral_gradient(V, dphi) result(res)
        !! Gradient of a symmetric function of the eigenvalues, phi(T) = phi(lam_1, lam_2, lam_3):
        !! dphi/dT = sum_a phi_a v_a (x) v_a, with phi_a = dphi/dlam_a supplied by the caller
        !! (de Souza Neto et al., 2008, Eq. A.7). Called by the CPB06 and Yld2004-18p gradients.
        real(real64), intent(in) :: V(3,3)   !! Eigenvectors from spectral_decomposition
        real(real64), intent(in) :: dphi(3)  !! phi_a = dphi/dlam_a
        type(ten_3D2Osym) :: res

        res%vals(1) = dphi(1)*V(1,1)*V(1,1) + dphi(2)*V(1,2)*V(1,2) + dphi(3)*V(1,3)*V(1,3)
        res%vals(2) = dphi(1)*V(2,1)*V(2,1) + dphi(2)*V(2,2)*V(2,2) + dphi(3)*V(2,3)*V(2,3)
        res%vals(3) = dphi(1)*V(3,1)*V(3,1) + dphi(2)*V(3,2)*V(3,2) + dphi(3)*V(3,3)*V(3,3)
        res%vals(4) = dphi(1)*V(1,1)*V(2,1) + dphi(2)*V(1,2)*V(2,2) + dphi(3)*V(1,3)*V(2,3)
        res%vals(5) = dphi(1)*V(2,1)*V(3,1) + dphi(2)*V(2,2)*V(3,2) + dphi(3)*V(2,3)*V(3,3)
        res%vals(6) = dphi(1)*V(1,1)*V(3,1) + dphi(2)*V(1,2)*V(3,2) + dphi(3)*V(1,3)*V(3,3)
    end function spectral_gradient

    pure function spectral_hessian(lam, V, dphi, d2phi) result(res)
        !! Hessian of phi(T) = phi(lam_1, lam_2, lam_3) (de Souza Neto et al., 2008, Box A.6;
        !! Miehe, 1998, Eq. 10):
        !! d2phi/dT2 = sum_ab phi_ab E_a (x) E_b + sum_(a<b) 2 theta_ab N_ab (x) N_ab,
        !! with E_a = v_a (x) v_a, N_ab = sym(v_a (x) v_b) and
        !! theta_ab = (phi_a - phi_b)/(lam_a - lam_b). For (nearly) repeated eigenvalues
        !! theta_ab takes its limit 1/2 (phi_aa + phi_bb) - phi_ab.
        !! Called by the CPB06 and Yld2004-18p Hessians.
        real(real64), intent(in) :: lam(3)      !! Eigenvalues from spectral_decomposition
        real(real64), intent(in) :: V(3,3)      !! Eigenvectors from spectral_decomposition
        real(real64), intent(in) :: dphi(3)     !! phi_a = dphi/dlam_a
        real(real64), intent(in) :: d2phi(3,3)  !! phi_ab = d2phi/dlam_a dlam_b
        type(ten_3D4O3sym) :: res

        type(ten_3D2Osym) :: E1, E2, E3, N
        real(real64) :: tol, theta

        E1 = sym_outer(V(:,1), V(:,1))
        E2 = sym_outer(V(:,2), V(:,2))
        E3 = sym_outer(V(:,3), V(:,3))
        tol = REPEATED_TOL*max(abs(lam(1)), abs(lam(3)))

        ! Terms of a = 1, then of the pairs (1,2) and (1,3), of a = 2, of the pair (2,3), of a = 3
        res = d2phi(1,1)*(.tdotsym. E1)
        if (abs(lam(1) - lam(2)) > tol) then
            theta = (dphi(1) - dphi(2))/(lam(1) - lam(2))
        else
            theta = 0.5D0*(d2phi(1,1) + d2phi(2,2)) - d2phi(1,2)
        end if
        N = sym_outer(V(:,1), V(:,2))
        res = res + (2.0D0*d2phi(1,2))*(E1 .tdotsym. E2) + (2.0D0*theta)*(.tdotsym. N)

        if (abs(lam(1) - lam(3)) > tol) then
            theta = (dphi(1) - dphi(3))/(lam(1) - lam(3))
        else
            theta = 0.5D0*(d2phi(1,1) + d2phi(3,3)) - d2phi(1,3)
        end if
        N = sym_outer(V(:,1), V(:,3))
        res = res + (2.0D0*d2phi(1,3))*(E1 .tdotsym. E3) + (2.0D0*theta)*(.tdotsym. N)

        res = res + d2phi(2,2)*(.tdotsym. E2)
        if (abs(lam(2) - lam(3)) > tol) then
            theta = (dphi(2) - dphi(3))/(lam(2) - lam(3))
        else
            theta = 0.5D0*(d2phi(2,2) + d2phi(3,3)) - d2phi(2,3)
        end if
        N = sym_outer(V(:,2), V(:,3))
        res = res + (2.0D0*d2phi(2,3))*(E2 .tdotsym. E3) + (2.0D0*theta)*(.tdotsym. N)

        res = res + d2phi(3,3)*(.tdotsym. E3)
    end function spectral_hessian

    pure function sym_outer(x, y) result(res)
        ! Symmetric part of the outer product of two vectors, sym(x (x) y); x = y gives x (x) x.
        ! Used by spectral_hessian
        real(real64), intent(in) :: x(3), y(3)
        type(ten_3D2Osym) :: res

        res%vals(1) = 0.5D0*(x(1)*y(1) + x(1)*y(1))
        res%vals(2) = 0.5D0*(x(2)*y(2) + x(2)*y(2))
        res%vals(3) = 0.5D0*(x(3)*y(3) + x(3)*y(3))
        res%vals(4) = 0.5D0*(x(1)*y(2) + y(1)*x(2))
        res%vals(5) = 0.5D0*(x(2)*y(3) + y(2)*x(3))
        res%vals(6) = 0.5D0*(x(1)*y(3) + y(1)*x(3))
    end function sym_outer

    pure function cross(x, y) result(res)
        ! Cross product x * y; used by spectral_decomposition
        real(real64), intent(in) :: x(3), y(3)
        real(real64) :: res(3)

        res = [x(2)*y(3) - x(3)*y(2), x(3)*y(1) - x(1)*y(3), x(1)*y(2) - x(2)*y(1)]
    end function cross

end module muscle_math_spectral_derivs