! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco@usach.cl>

module muscle_tensor_ops_contraction_single
    !! This module defines the overloaded single contraction (or matrix multiplication) operator `*`
    !! for mixed second-order and identity tensors.
    !!
    !! Mathematically, the single contraction of two second-order tensors \(\mathbf{A}\) and \(\mathbf{B}\) 
    !! represents their standard dot product (matrix multiplication):
    !! \[ \mathbf{C} = \mathbf{A} \cdot \mathbf{B} \quad \implies \quad C_{ij} = A_{ik} B_{kj} \]
    !!
    !! Since the product of two symmetric second-order tensors is generally non-symmetric (unless they commute),
    !! these mixed operations safely return a general second-order tensor (`ten_3D2O`) with 9 components.

    use, intrinsic :: iso_fortran_env, only : real64
    
    ! Import necessary types
    use muscle_tensor_iden_2o
    use muscle_tensor_iden_2os
    use muscle_tensor_2d2osym
    use muscle_tensor_3d2o
    use muscle_tensor_3d2osym

    implicit none
    private

    ! =========================================================================
    ! PUBLIC INTERFACES
    ! =========================================================================
    
    public :: operator(*)
    interface operator (*)
        ! --- Tensor * Tensor (Single Contraction / Matrix Product) ---
        module procedure dot_3D2Osym_3D2Osym
        module procedure dot_3D2O_3D2Osym
        module procedure dot_3D2Osym_3D2O

        ! --- Identity * Tensor (Single Contraction) ---
        module procedure mul_2D2Osym_I2O
        module procedure mul_I2O_2D2Osym
        module procedure mul_3D2O_I2O
        module procedure mul_I2O_3D2O
        module procedure mul_3D2Osym_I2O
        module procedure mul_I2O_3D2Osym
        module procedure mul_I2OS_3D2O
        module procedure mul_3D2O_I2OS
    end interface

contains

    ! =========================================================================
    ! 1. TENSOR * TENSOR (SINGLE CONTRACTION)
    ! =========================================================================

    pure function dot_3D2Osym_3D2Osym(a, b) result(res)
        !! Computes the single tensor contraction (matrix product) between two 
        !! symmetric second-order tensors.
        !!
        !! Mathematically: \( res_{ij} = a_{ik} b_{kj} \)
        !!
        !! Note: The product of two symmetric tensors is only symmetric if they 
        !! commute. Thus, the result is safely returned as a general `ten_3D2O` tensor.
        !!
        !! Performance: Manual loop unrolling mapping the 6-component Voigt 
        !! arrays to a 9-component column-major array.
        implicit none
        type(ten_3D2Osym), intent(in) :: a
            !! First symmetric second-order tensor \(\mathbf{A}\) (6 components)
        type(ten_3D2Osym), intent(in) :: b
            !! Second symmetric second-order tensor \(\mathbf{B}\) (6 components)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\) (9 components)
        
        ! Column 1 (xx, yx, zx)
        res%vals(1) = a%vals(1)*b%vals(1) + a%vals(4)*b%vals(4) + a%vals(6)*b%vals(6)
        res%vals(2) = a%vals(4)*b%vals(1) + a%vals(2)*b%vals(4) + a%vals(5)*b%vals(6)
        res%vals(3) = a%vals(6)*b%vals(1) + a%vals(5)*b%vals(4) + a%vals(3)*b%vals(6)
        
        ! Column 2 (xy, yy, zy)
        res%vals(4) = a%vals(1)*b%vals(4) + a%vals(4)*b%vals(2) + a%vals(6)*b%vals(5)
        res%vals(5) = a%vals(4)*b%vals(4) + a%vals(2)*b%vals(2) + a%vals(5)*b%vals(5)
        res%vals(6) = a%vals(6)*b%vals(4) + a%vals(5)*b%vals(2) + a%vals(3)*b%vals(5)
        
        ! Column 3 (xz, yz, zz)
        res%vals(7) = a%vals(1)*b%vals(6) + a%vals(4)*b%vals(5) + a%vals(6)*b%vals(3)
        res%vals(8) = a%vals(4)*b%vals(6) + a%vals(2)*b%vals(5) + a%vals(5)*b%vals(3)
        res%vals(9) = a%vals(6)*b%vals(6) + a%vals(5)*b%vals(5) + a%vals(3)*b%vals(3)
    end function dot_3D2Osym_3D2Osym

    pure function dot_3D2O_3D2Osym(a, b) result(res)
        !! Computes the single tensor contraction between a general and a 
        !! symmetric second-order tensor.
        !!
        !! Mathematically: \( res_{ij} = a_{ik} b_{kj} \)
        !!
        !! Performance: Manual loop unrolling mapping the 9-component array 
        !! and the 6-component array directly.
        implicit none
        type(ten_3D2O), intent(in) :: a     
            !! General second-order tensor \(\mathbf{A}\) (9 components)
        type(ten_3D2Osym), intent(in) :: b  
            !! Symmetric second-order tensor \(\mathbf{B}\) (6 components)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\) (9 components)
        
        ! Column 1
        res%vals(1) = a%vals(1)*b%vals(1) + a%vals(4)*b%vals(4) + a%vals(7)*b%vals(6)
        res%vals(2) = a%vals(2)*b%vals(1) + a%vals(5)*b%vals(4) + a%vals(8)*b%vals(6)
        res%vals(3) = a%vals(3)*b%vals(1) + a%vals(6)*b%vals(4) + a%vals(9)*b%vals(6)
        
        ! Column 2
        res%vals(4) = a%vals(1)*b%vals(4) + a%vals(4)*b%vals(2) + a%vals(7)*b%vals(5)
        res%vals(5) = a%vals(2)*b%vals(4) + a%vals(5)*b%vals(2) + a%vals(8)*b%vals(5)
        res%vals(6) = a%vals(3)*b%vals(4) + a%vals(6)*b%vals(2) + a%vals(9)*b%vals(5)
        
        ! Column 3
        res%vals(7) = a%vals(1)*b%vals(6) + a%vals(4)*b%vals(5) + a%vals(7)*b%vals(3)
        res%vals(8) = a%vals(2)*b%vals(6) + a%vals(5)*b%vals(5) + a%vals(8)*b%vals(3)
        res%vals(9) = a%vals(3)*b%vals(6) + a%vals(6)*b%vals(5) + a%vals(9)*b%vals(3)
    end function dot_3D2O_3D2Osym

    pure function dot_3D2Osym_3D2O(a, b) result(res)
        !! Computes the single tensor contraction between a symmetric and a 
        !! general second-order tensor.
        !!
        !! Mathematically: \( res_{ij} = a_{ik} b_{kj} \)
        !!
        !! Performance: Manual loop unrolling.
        implicit none
        type(ten_3D2Osym), intent(in) :: a  
            !! Symmetric second-order tensor \(\mathbf{A}\) (6 components)
        type(ten_3D2O), intent(in) :: b     
            !! General second-order tensor \(\mathbf{B}\) (9 components)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\) (9 components)
        
        ! Column 1
        res%vals(1) = a%vals(1)*b%vals(1) + a%vals(4)*b%vals(2) + a%vals(6)*b%vals(3)
        res%vals(2) = a%vals(4)*b%vals(1) + a%vals(2)*b%vals(2) + a%vals(5)*b%vals(3)
        res%vals(3) = a%vals(6)*b%vals(1) + a%vals(5)*b%vals(2) + a%vals(3)*b%vals(3)
        
        ! Column 2
        res%vals(4) = a%vals(1)*b%vals(4) + a%vals(4)*b%vals(5) + a%vals(6)*b%vals(6)
        res%vals(5) = a%vals(4)*b%vals(4) + a%vals(2)*b%vals(5) + a%vals(5)*b%vals(6)
        res%vals(6) = a%vals(6)*b%vals(4) + a%vals(5)*b%vals(5) + a%vals(3)*b%vals(6)
        
        ! Column 3
        res%vals(7) = a%vals(1)*b%vals(7) + a%vals(4)*b%vals(8) + a%vals(6)*b%vals(9)
        res%vals(8) = a%vals(4)*b%vals(7) + a%vals(2)*b%vals(8) + a%vals(5)*b%vals(9)
        res%vals(9) = a%vals(6)*b%vals(7) + a%vals(5)*b%vals(8) + a%vals(3)*b%vals(9)
    end function dot_3D2Osym_3D2O

    ! =========================================================================
    ! 2. IDENTITY * TENSOR (SINGLE CONTRACTION)
    ! =========================================================================

    pure function mul_I2O_2D2Osym(I2, a) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{I} \cdot \mathbf{A} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
    end function mul_I2O_2D2Osym

    pure function mul_2D2Osym_I2O(a, I2) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{A} \cdot \mathbf{I} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
    end function mul_2D2Osym_I2O

    pure function mul_I2O_3D2O(I2, a) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{I} \cdot \mathbf{A} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
    end function mul_I2O_3D2O

    pure function mul_3D2O_I2O(a, I2) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{A} \cdot \mathbf{I} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
    end function mul_3D2O_I2O

    pure function mul_I2O_3D2Osym(I2, a) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{I} \cdot \mathbf{A} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
    end function mul_I2O_3D2Osym

    pure function mul_3D2Osym_I2O(a, I2) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{A} \cdot \mathbf{I} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
    end function mul_3D2Osym_I2O

    pure function mul_I2OS_3D2O(I2, a) result(res)
        !! Computes the scaled contraction \(\mathbf{res} = c\mathbf{I} \cdot \mathbf{A} = c\mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = I2%val * a%vals
    end function mul_I2OS_3D2O

    pure function mul_3D2O_I2OS(a, I2) result(res)
        !! Computes the scaled contraction \(\mathbf{res} = \mathbf{A} \cdot c\mathbf{I} = c\mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = I2%val * a%vals
    end function mul_3D2O_I2OS
end module muscle_tensor_ops_contraction_single