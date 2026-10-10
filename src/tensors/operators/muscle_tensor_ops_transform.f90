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
    use muscle_tensor_2d2o
    use muscle_tensor_2d2osym
    use muscle_tensor_2d4o2sym
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
        ! --- 2D: General Tensor (2D2O) transforming 2D2Osym / 2D4O2sym ---
        module procedure transform_2D2O_2D2Osym
        module procedure transform_2D2O_2D4O2sym
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

    pure function transform_2D2O_2D2Osym(F, S) result(res)
        !! Computes the 2D congruence transformation \(\mathbf{res} = \mathbf{F} \cdot \mathbf{S} \cdot \mathbf{F}^T\).
        !!
        !! Since \(F_{13} = F_{23} = F_{31} = F_{32} = 0\) and \(S_{13} = S_{23} = 0\), the in-plane
        !! 2x2 blocks transform as matrices and the out-of-plane component as a scalar:
        !! \( res_{33} = F_{33}^2\,S_{33} \).
        !!
        !! Storage reference:
        !! - `F` (ten_2D2O): (11, 21, 12, 22, 33) -> Indices (1 to 5)
        !! - `S`, `res` (ten_2D2Osym): (11, 22, 33, 12) -> Indices (1 to 4)
        implicit none
        type(ten_2D2O), intent(in) :: F
            !! The transformation tensor (e.g. plane/axisymmetric deformation gradient).
        type(ten_2D2Osym), intent(in) :: S
            !! The tensor to be transformed (e.g. PK2 stress), symmetric 2nd order.
        type(ten_2D2Osym) :: res
            !! The resulting transformed tensor, symmetric 2nd order.

        real(real64) :: F11, F21, F12, F22, F33
        real(real64) :: S11, S22, S33, S12
        real(real64) :: T11, T12, T21, T22

        F11 = F%vals(1); F21 = F%vals(2); F12 = F%vals(3); F22 = F%vals(4); F33 = F%vals(5)
        S11 = S%vals(1); S22 = S%vals(2); S33 = S%vals(3); S12 = S%vals(4)

        ! Intermediate in-plane tensor T = F * S
        T11 = F11*S11 + F12*S12
        T12 = F11*S12 + F12*S22
        T21 = F21*S11 + F22*S12
        T22 = F21*S12 + F22*S22

        ! res = T * F^T  (res_ij = T_iK * F_jK)
        res%vals(1) = T11*F11 + T12*F12 ! xx
        res%vals(2) = T21*F21 + T22*F22 ! yy
        res%vals(3) = F33*F33*S33       ! zz
        res%vals(4) = T11*F21 + T12*F22 ! xy
    end function transform_2D2O_2D2Osym

    pure function transform_2D2O_2D4O2sym(A, C_in) result(res)
        !! Computes the 2D fourth-order tensor transformation
        !! \( res_{ijkl} = A_{iI} A_{jJ} A_{kK} A_{lL} C_{IJKL} \).
        !!
        !! Thanks to the minor symmetries, the transformation is evaluated directly on
        !! the 4x4 Voigt matrices (11, 22, 33, 12) as
        !! \[ \mathbf{R} = \mathbf{T}\,\mathbf{C}\,\mathbf{T}^T, \qquad
        !!    T_{(ij),(ab)} = A_{ia}A_{jb} + [a \neq b]\,A_{ib}A_{ja}, \]
        !! which, with \(A_{13} = A_{23} = A_{31} = A_{32} = 0\), reads
        !! \[ \mathbf{T} = \begin{bmatrix}
        !!    A_{11}^2 & A_{12}^2 & 0 & 2A_{11}A_{12} \\
        !!    A_{21}^2 & A_{22}^2 & 0 & 2A_{21}A_{22} \\
        !!    0 & 0 & A_{33}^2 & 0 \\
        !!    A_{11}A_{21} & A_{12}A_{22} & 0 & A_{11}A_{22} + A_{12}A_{21}
        !!    \end{bmatrix}. \]
        !! The two products exploit the zero pattern of \(\mathbf{T}\) (about 80
        !! multiplications instead of the \(4 \cdot 3^8\) of the full index contraction;
        !! ~150x faster at -O3).
        implicit none
        type(ten_2D2O), intent(in)     :: A
            !! The transformation tensor, general 2D 2nd order.
        type(ten_2D4O2sym), intent(in) :: C_in
            !! The fourth-order tensor to be transformed (minor symmetries).
        type(ten_2D4O2sym)             :: res
            !! The transformed fourth-order tensor.

        real(real64) :: W(4,4)
        real(real64) :: A11, A21, A12, A22, A33
        real(real64) :: T11, T12, T14, T21, T22, T24, T33, T41, T42, T44

        ! ten_2D2O storage: (11, 21, 12, 22, 33)
        A11 = A%vals(1); A21 = A%vals(2); A12 = A%vals(3); A22 = A%vals(4); A33 = A%vals(5)

        ! Non-zero entries of T
        T11 = A11*A11; T12 = A12*A12; T14 = 2.0D0*A11*A12
        T21 = A21*A21; T22 = A22*A22; T24 = 2.0D0*A21*A22
        T33 = A33*A33
        T41 = A11*A21; T42 = A12*A22; T44 = A11*A22 + A12*A21

        ! W = T * C
        W(1,:) = T11*C_in%vals(1,:) + T12*C_in%vals(2,:) + T14*C_in%vals(4,:)
        W(2,:) = T21*C_in%vals(1,:) + T22*C_in%vals(2,:) + T24*C_in%vals(4,:)
        W(3,:) = T33*C_in%vals(3,:)
        W(4,:) = T41*C_in%vals(1,:) + T42*C_in%vals(2,:) + T44*C_in%vals(4,:)

        ! R = W * T^T
        res%vals(:,1) = W(:,1)*T11 + W(:,2)*T12 + W(:,4)*T14
        res%vals(:,2) = W(:,1)*T21 + W(:,2)*T22 + W(:,4)*T24
        res%vals(:,3) = W(:,3)*T33
        res%vals(:,4) = W(:,1)*T41 + W(:,2)*T42 + W(:,4)*T44
    end function transform_2D2O_2D4O2sym

end module muscle_tensor_ops_transform
