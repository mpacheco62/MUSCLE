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

    real(real64), parameter :: EPS_TOL = 1.0D-11
    real(real64), parameter :: EPS_TOL_SECOND = 1.0D-6
    !! Safe tolerance threshold to prevent division by zero at repeated roots (singularities)

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
        call M_tensor%init( &
            xxxx=2.0D0*Sigma%xx(), xxyy=0.0D0, xxzz=0.0D0, xxxy=Sigma%xy(), xxyz=0.0D0, xxxz=Sigma%xz(), &
            yyyy=2.0D0*Sigma%yy(), yyzz=0.0D0, yyxy=Sigma%xy(), yyyz=Sigma%yz(), yyxz=0.0D0, &
            zzzz=2.0D0*Sigma%zz(), zzxy=0.0D0, zzyz=Sigma%yz(), zzxz=Sigma%xz(), &
            xyxy=0.5D0*(Sigma%xx()+Sigma%yy()), xyyz=0.5D0*Sigma%xz(), xyxz=0.5D0*Sigma%yz(), &
            yzyz=0.5D0*(Sigma%yy()+Sigma%zz()), yzxz=0.5D0*Sigma%xy(), &
            xzxz=0.5D0*(Sigma%xx()+Sigma%zz()) &
        )

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

end module muscle_math_spectral_derivs