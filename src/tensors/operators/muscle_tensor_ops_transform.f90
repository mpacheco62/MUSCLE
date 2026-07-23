! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_ops_transform
    !! Module muscle_tensor_ops_transform
    !! ==================================
    !!
    !! This module defines the overloaded congruence transformation operator `.transform.`
    !! for the tensor engine.
    !!
    !! Mathematically, the congruence transformation of a symmetric tensor \(\mathbf{S}\)
    !! by a general tensor \(\mathbf{F}\) (such as the deformation gradient) is defined as:
    !! \[ \boldsymbol{\sigma} = \mathbf{F} \cdot \mathbf{S} \cdot \mathbf{F}^T \]
    !! \[ \sigma_{ij} = F_{iK} S_{KL} F_{jL} \]
    !!
    !! This operation is the fundamental algebraic building block for physical pull-back 
    !! and push-forward operations in finite strain continuum mechanics (e.g., transforming 
    !! 2nd Piola-Kirchhoff stress to Cauchy stress).
    !!
    !! Design Rationale
    !! ----------------
    !! This operator isolates the pure algebraic projection from any physical meaning 
    !! (like the Jacobian determinant \(J\)). By keeping this inside the tensor engine, 
    !! the loops are manually unrolled to operate directly on the 9-component (column-major) 
    !! and 6-component (Voigt) storage formats. This avoids the creation of temporary 
    !! 9-component intermediate tensors, ensuring maximum computational performance in 
    !! implicit/explicit FEA solvers.

    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensor_3d2o
    use muscle_tensor_3d2osym
    use muscle_tensor_3d4o2sym

    implicit none
    private

    ! =========================================================================
    ! PUBLIC INTERFACES
    ! =========================================================================
    
    public :: operator(.transform.)
    interface operator (.transform.)
        ! --- General Tensor (3D2O) transforming Symmetric Tensor (3D2Osym) ---
        module procedure transform_3D2O_3D2Osym
        module procedure transform_3D2O_3D4O2sym
    end interface

contains

    pure function transform_3D2O_3D2Osym(F, S) result(res)
        !! Computes the congruence transformation \(\mathbf{res} = \mathbf{F} \cdot \mathbf{S} \cdot \mathbf{F}^T\).
        !!
        !! This function fully unrolls the double dot product into scalar operations.
        !! It calculates the intermediate product \(\mathbf{T} = \mathbf{F} \cdot \mathbf{S}\), 
        !! and then immediately computes the symmetric part of \(\mathbf{T} \cdot \mathbf{F}^T\).
        !!
        !! Storage reference:
        !! - `F` (ten_3D2O): (11, 21, 31, 12, 22, 32, 13, 23, 33) -> Indices (1 to 9)
        !! - `S` (ten_3D2Osym): (11, 22, 33, 12, 23, 13) -> Indices (1 to 6)
        !! - `res` (ten_3D2Osym): (11, 22, 33, 12, 23, 13) -> Indices (1 to 6)
        implicit none
        type(ten_3D2O), intent(in) :: F
            !! The transformation tensor (e.g. Deformation Gradient), general 2nd order.
        type(ten_3D2Osym), intent(in) :: S
            !! The tensor to be transformed (e.g. PK2 Stress), symmetric 2nd order.
        type(ten_3D2Osym) :: res
            !! The resulting transformed tensor, symmetric 2nd order.

        real(real64) :: F11, F21, F31, F12, F22, F32, F13, F23, F33
        real(real64) :: S11, S22, S33, S12, S23, S13
        real(real64) :: T11, T12, T13, T21, T22, T23, T31, T32, T33

        ! 1. Map arrays to local scalars for optimal register allocation
        F11 = F%vals(1); F21 = F%vals(2); F31 = F%vals(3)
        F12 = F%vals(4); F22 = F%vals(5); F32 = F%vals(6)
        F13 = F%vals(7); F23 = F%vals(8); F33 = F%vals(9)

        S11 = S%vals(1); S22 = S%vals(2); S33 = S%vals(3)
        S12 = S%vals(4); S23 = S%vals(5); S13 = S%vals(6)

        ! 2. Compute Intermediate Tensor T = F * S
        ! Since S is symmetric, S_21 = S_12 (S12), S_32 = S_23 (S23), S_31 = S_13 (S13)
        T11 = F11*S11 + F12*S12 + F13*S13
        T12 = F11*S12 + F12*S22 + F13*S23
        T13 = F11*S13 + F12*S23 + F13*S33

        T21 = F21*S11 + F22*S12 + F23*S13
        T22 = F21*S12 + F22*S22 + F23*S23
        T23 = F21*S13 + F22*S23 + F23*S33

        T31 = F31*S11 + F32*S12 + F33*S13
        T32 = F31*S12 + F32*S22 + F33*S23
        T33 = F31*S13 + F32*S23 + F33*S33

        ! 3. Compute res = T * F^T (Only the 6 symmetric components are evaluated)
        ! res_ij = T_iK * F^T_Kj = T_iK * F_jK
        
        ! Normal components (xx, yy, zz)
        res%vals(1) = T11*F11 + T12*F12 + T13*F13
        res%vals(2) = T21*F21 + T22*F22 + T23*F23
        res%vals(3) = T31*F31 + T32*F32 + T33*F33

        ! Shear components (xy, yz, xz)
        res%vals(4) = T11*F21 + T12*F22 + T13*F23
        res%vals(5) = T21*F31 + T22*F32 + T23*F33
        res%vals(6) = T11*F31 + T12*F32 + T13*F33

    end function transform_3D2O_3D2Osym


    pure function transform_3D2O_3D4O2sym(A, C_in) result(res)
        !! Computes 4th order tensor transformation: res_ijkl = A_iI * A_jJ * A_kK * A_lL * C_IJKL
        type(ten_3D2O), intent(in)     :: A
        type(ten_3D4O2sym), intent(in) :: C_in
        type(ten_3D4O2sym)             :: res

        real(real64) :: A3(3,3), C4_in(3,3,3,3), C4_out(3,3,3,3)
        integer, parameter :: v_i(6) = (/1, 2, 3, 1, 2, 1/)
        integer, parameter :: v_j(6) = (/1, 2, 3, 2, 3, 3/)
        integer :: I, J, i_idx, j_idx, k_idx, l_idx
        integer :: I_mat, J_mat, K_mat, L_mat

        ! 1. Unpack A (ten_3D2O column-major) to 3x3 matrix
        A3(1,1)=A%vals(1); A3(2,1)=A%vals(2); A3(3,1)=A%vals(3)
        A3(1,2)=A%vals(4); A3(2,2)=A%vals(5); A3(3,2)=A%vals(6)
        A3(1,3)=A%vals(7); A3(2,3)=A%vals(8); A3(3,3)=A%vals(9)

        ! 2. Unpack C_in (6x6 Voigt matrix) to 3x3x3x3 tensor
        C4_in = 0.0D0
        do I = 1, 6
            do J = 1, 6
                i_idx = v_i(I); j_idx = v_j(I)
                k_idx = v_i(J); l_idx = v_j(J)

                ! Apply minor symmetries
                C4_in(i_idx, j_idx, k_idx, l_idx) = C_in%vals(I, J)
                C4_in(j_idx, i_idx, k_idx, l_idx) = C_in%vals(I, J)
                C4_in(i_idx, j_idx, l_idx, k_idx) = C_in%vals(I, J)
                C4_in(j_idx, i_idx, l_idx, k_idx) = C_in%vals(I, J)
            end do
        end do

        ! 3. 4-Index Contraction: C_out_ijkl = A_iI * A_jJ * A_kK * A_lL * C_IJKL
        C4_out = 0.0D0
        do i_idx = 1, 3
        do j_idx = 1, 3
        do k_idx = 1, 3
        do l_idx = 1, 3
            do I_mat = 1, 3
            do J_mat = 1, 3
            do K_mat = 1, 3
            do L_mat = 1, 3
                C4_out(i_idx, j_idx, k_idx, l_idx) = C4_out(i_idx, j_idx, k_idx, l_idx) + &
                    A3(i_idx, I_mat) * A3(j_idx, J_mat) * A3(k_idx, K_mat) * A3(l_idx, L_mat) * &
                    C4_in(I_mat, J_mat, K_mat, L_mat)
            end do
            end do
            end do
            end do
        end do
        end do
        end do
        end do

        ! 4. Pack 3x3x3x3 tensor back to 6x6 Voigt res%vals
        do I = 1, 6
            do J = 1, 6
                i_idx = v_i(I); j_idx = v_j(I)
                k_idx = v_i(J); l_idx = v_j(J)
                res%vals(I, J) = C4_out(i_idx, j_idx, k_idx, l_idx)
            end do
        end do

    end function transform_3D2O_3D4O2sym

end module muscle_tensor_ops_transform