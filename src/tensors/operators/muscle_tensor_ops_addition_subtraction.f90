! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_ops_addition_subtraction
    !! This module consolidates all mixed-type addition and subtraction operators
    !! for second and fourth-order tensors, as well as standard and scaled identity tensors.
    !!
    !! It overloads the intrinsic operators:
    !! - `+` : Mixed-type addition.
    !! - `-` : Mixed-type subtraction and unary negation.
    !!
    !! Grouping these operations here reduces boilerplate interface code and centralizes
    !! the algebraic rules of the tensor engine.
    
    use, intrinsic :: iso_fortran_env, only : real64
    
    ! Import all necessary tensor and identity types
    use muscle_tensor_iden_2o
    use muscle_tensor_iden_2os
    use muscle_tensor_iden_4o3t
    use muscle_tensor_iden_4o3ts
    use muscle_tensor_iden_4o4t
    use muscle_tensor_iden_4o4ts
    use muscle_tensor_2d2osym
    use muscle_tensor_2d4o3sym
    use muscle_tensor_3d2o
    use muscle_tensor_3d2osym
    use muscle_tensor_3d4o2sym
    use muscle_tensor_3d4o3sym

    implicit none
    private

    ! =========================================================================
    ! PUBLIC INTERFACES
    ! =========================================================================
    
    public :: operator(+)
    interface operator (+)
        ! --- 2nd Order Identity + 2nd Order Tensor ---    
        module procedure sum_I2O_3D2Osym
        module procedure sum_3D2Osym_I2O
        module procedure sum_I2O_2D2Osym
        module procedure sum_2D2Osym_I2O
        module procedure sum_I2O_3D2O
        module procedure sum_3D2O_I2O

        ! --- Scaled 2nd Order Identity + 2nd Order Tensor ---
        module procedure sum_I2OS_3D2Osym
        module procedure sum_3D2Osym_I2OS
        module procedure sum_I2OS_2D2Osym
        module procedure sum_2D2Osym_I2OS
        module procedure sum_I2OS_3D2O
        module procedure sum_3D2O_I2OS

        ! --- 4th Order Identity + 4th Order Tensor ---
        module procedure sum_I4O3TS_3D4O2sym
        module procedure sum_3D4O2sym_I4O3TS
        module procedure sum_I4O3TS_3D4O3sym
        module procedure sum_3D4O3sym_I4O3TS
        module procedure sum_I4O4TS_3D4O2sym
        module procedure sum_3D4O2sym_I4O4TS
        module procedure sum_I4O4TS_3D4O3sym
        module procedure sum_3D4O3sym_I4O4TS
        module procedure sum_I4O3T_3D4O2sym
        module procedure sum_3D4O2sym_I4O3T
        module procedure sum_I4O3T_3D4O3sym
        module procedure sum_3D4O3sym_I4O3T
        module procedure sum_I4O4T_3D4O2sym
        module procedure sum_3D4O2sym_I4O4T
        module procedure sum_I4O4T_3D4O3sym
        module procedure sum_3D4O3sym_I4O4T

        ! --- Identity + Identity ---
        module procedure sum_I4O3TS_I4O4T
        module procedure sum_I4O4T_I4O3TS
        module procedure sum_I4O3TS_I4O4TS
        module procedure sum_I4O4TS_I4O3TS
        module procedure sum_I4O3T_I4O4T
        module procedure sum_I4O4T_I4O3T
        module procedure sum_I4O3T_I4O4TS
        module procedure sum_I4O4TS_I4O3T

        ! --- 2nd Order Tensor + 2nd Order Tensor ---
        module procedure sum_3D2O_3D2Osym
        module procedure sum_3D2Osym_3D2O
    end interface

    public :: operator(-)
    interface operator (-)
        ! --- 2nd Order Identity - 2nd Order Tensor ---
        module procedure sub_I2O_3D2Osym
        module procedure sub_3D2Osym_I2O
        module procedure sub_I2O_2D2Osym
        module procedure sub_2D2Osym_I2O
        module procedure sub_I2O_3D2O
        module procedure sub_3D2O_I2O
        
        ! --- Scaled 2nd Order Identity - 2nd Order Tensor ---
        module procedure sub_I2OS_3D2Osym
        module procedure sub_3D2Osym_I2OS
        module procedure sub_I2OS_2D2Osym
        module procedure sub_2D2Osym_I2OS
        module procedure sub_I2OS_3D2O
        module procedure sub_3D2O_I2OS

        ! ! --- 4th Order Identity - 4th Order Tensor ---
        module procedure sub_I4O3TS_3D4O2sym
        module procedure sub_3D4O2sym_I4O3TS
        module procedure sub_I4O3TS_3D4O3sym
        module procedure sub_3D4O3sym_I4O3TS
        module procedure sub_I4O4TS_3D4O2sym
        module procedure sub_3D4O2sym_I4O4TS
        module procedure sub_I4O4TS_3D4O3sym
        module procedure sub_3D4O3sym_I4O4TS
        module procedure sub_I4O3T_3D4O2sym
        module procedure sub_3D4O2sym_I4O3T
        module procedure sub_I4O3T_3D4O3sym
        module procedure sub_3D4O3sym_I4O3T
        module procedure sub_I4O4T_3D4O2sym
        module procedure sub_3D4O2sym_I4O4T
        module procedure sub_I4O4T_3D4O3sym
        module procedure sub_3D4O3sym_I4O4T

        ! --- Identity - Identity ---
        module procedure sub_I4O3TS_I4O4T
        module procedure sub_I4O4T_I4O3TS
        module procedure sub_I4O3TS_I4O4TS
        module procedure sub_I4O4TS_I4O3TS
        module procedure sub_I4O3T_I4O4T
        module procedure sub_I4O4T_I4O3T
        module procedure sub_I4O3T_I4O4TS
        module procedure sub_I4O4TS_I4O3T
    end interface

contains

    ! =========================================================================
    ! IMPLEMENTATIONS: ADDITION (+)
    ! =========================================================================

    ! --- 2nd Order Identity + 2nd Order Tensor ---    
    !*************************************************************************

    pure function sum_I2O_3D2Osym(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
    end function sum_I2O_3D2Osym

    pure function sum_3D2Osym_I2O(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
    end function sum_3D2Osym_I2O

    pure function sum_I2O_2D2Osym(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
    end function sum_I2O_2D2Osym

    pure function sum_2D2Osym_I2O(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
    end function sum_2D2Osym_I2O

    pure function sum_I2O_3D2O(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) + 1.0D0
        res%vals(5) = res%vals(5) + 1.0D0
        res%vals(9) = res%vals(9) + 1.0D0
    end function sum_I2O_3D2O

    pure function sum_3D2O_I2O(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) + 1.0D0
        res%vals(5) = res%vals(5) + 1.0D0
        res%vals(9) = res%vals(9) + 1.0D0
    end function sum_3D2O_I2O

    ! --- Scaled 2nd Order Identity - 2nd Order Tensor ---
    !*************************************************************************
    pure function sum_I2OS_3D2Osym(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = c\mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_I2OS_3D2Osym

    pure function sum_3D2Osym_I2OS(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_3D2Osym_I2OS

    pure function sum_I2OS_2D2Osym(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = c\mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_I2OS_2D2Osym

    pure function sum_2D2Osym_I2OS(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_2D2Osym_I2OS

    pure function sum_I2OS_3D2O(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = c\mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) + I2%val
        res%vals(5) = res%vals(5) + I2%val
        res%vals(9) = res%vals(9) + I2%val
    end function sum_I2OS_3D2O

    pure function sum_3D2O_I2OS(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) + I2%val
        res%vals(5) = res%vals(5) + I2%val
        res%vals(9) = res%vals(9) + I2%val
    end function sum_3D2O_I2OS

    ! --- 4th Order Identity + 4th Order Tensor ---
    !*************************************************************************

    pure function sum_I4O3TS_3D4O2sym(I4S, a) result(res)
        !! Computes the sum \(\mathbb{res} = c\mathbb{I}_{3T} + \mathbb{A}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I4S
            !! Scaled 4th-order identity tensor (\(c \cdot \delta_{ij}\delta_{kl}\)).
        type(ten_3D4O2sym), intent(in) :: a
            !! 4th-order tensor with minor symmetries (6x6 matrix).
        type(ten_3D4O2sym) :: res
            !! Resulting 4th-order tensor.
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + I4S%val
    end function sum_I4O3TS_3D4O2sym

    pure function sum_3D4O2sym_I4O3TS(a, I4S) result(res)
        !! Computes the sum \(\mathbb{res} = \mathbb{A} + c\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I4S
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + I4S%val
    end function sum_3D4O2sym_I4O3TS

    pure function sum_I4O3TS_3D4O3sym(I4S, a) result(res)
        !! Computes the sum \(\mathbb{res} = c\mathbb{I}_{3T} + \mathbb{A}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I4S
            !! Scaled 4th-order identity tensor (\(c \cdot \delta_{ij}\delta_{kl}\)).
        type(ten_3D4O3sym), intent(in) :: a
            !! Fully symmetric 4th-order tensor (21 components).
        type(ten_3D4O3sym) :: res
            !! Resulting 4th-order tensor.
        
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I4S%val
        res%vals(7:8) = res%vals(7:8) + I4S%val
        res%vals(12)  = res%vals(12)  + I4S%val
    end function sum_I4O3TS_3D4O3sym

    pure function sum_3D4O3sym_I4O3TS(a, I4S) result(res)
        !! Computes the sum \(\mathbb{res} = \mathbb{A} + c\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I4S
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I4S%val
        res%vals(7:8) = res%vals(7:8) + I4S%val
        res%vals(12)  = res%vals(12)  + I4S%val
    end function sum_3D4O3sym_I4O3TS

    pure function sum_I4O4TS_3D4O2sym(I4S, a) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}^S + \mathbb{A}\).
        implicit none
        type(iden_4O4TS),  intent(in) :: I4S
            !! Scaled fourth-order symmetric identity tensor.
        type(ten_3D4O2sym), intent(in) :: a
            !! Fourth-order tensor with minor symmetries (6x6 matrix).
        type(ten_3D4O2sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = a%vals
        ! Add components to the Voigt diagonal
        res%vals(1,1) = res%vals(1,1) + c
        res%vals(2,2) = res%vals(2,2) + c
        res%vals(3,3) = res%vals(3,3) + c
        res%vals(4,4) = res%vals(4,4) + half_c
        res%vals(5,5) = res%vals(5,5) + half_c
        res%vals(6,6) = res%vals(6,6) + half_c
    end function sum_I4O4TS_3D4O2sym

    pure function sum_3D4O2sym_I4O4TS(a, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{A} + c\mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
        type(iden_4O4TS),  intent(in) :: I4S
        type(ten_3D4O2sym)             :: res
        
        ! Commutative property: delegate to the primary implementation
        res = sum_I4O4TS_3D4O2sym(I4S, a)
    end function sum_3D4O2sym_I4O4TS

    pure function sum_I4O4TS_3D4O3sym(I4S, a) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}^S + \mathbf{A}\).
        implicit none
        type(iden_4O4TS),   intent(in) :: I4S
            !! Scaled fourth-order symmetric identity tensor.
        type(ten_3D4O3sym), intent(in) :: a
            !! Fully symmetric fourth-order tensor (21 components).
        type(ten_3D4O3sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = a%vals
        ! Add scaled identity components to the Voigt diagonal
        res%vals(1:3) = res%vals(1:3) + c
        res%vals(4:6) = res%vals(4:6) + half_c
    end function sum_I4O4TS_3D4O3sym

    pure function sum_3D4O3sym_I4O4TS(a, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} + c\mathbb{I}^S\).
        implicit none
        type(ten_3D4O3sym), intent(in) :: a
        type(iden_4O4TS),   intent(in) :: I4S
        type(ten_3D4O3sym)             :: res
        
        ! Commutative property: delegate to the primary implementation
        res = sum_I4O4TS_3D4O3sym(I4S, a)
    end function sum_3D4O3sym_I4O4TS

    pure function sum_I4O3T_3D4O2sym(I4, a) result(res)
        !! Computes the sum \(\mathbb{res} = \mathbb{I} + \mathbb{A}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
            !! Standard 4th-order identity tensor (\(\delta_{ij}\delta_{kl}\)).
        type(ten_3D4O2sym), intent(in) :: a
            !! 4th-order tensor with minor symmetries.
        type(ten_3D4O2sym) :: res
            !! Resulting 4th-order tensor.
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1.0D0
    end function sum_I4O3T_3D4O2sym

    pure function sum_3D4O2sym_I4O3T(a, I4) result(res)
        !! Computes the sum \(\mathbb{res} = \mathbb{A} + \mathbb{I}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1.0D0
    end function sum_3D4O2sym_I4O3T

    pure function sum_I4O3T_3D4O3sym(I4, a) result(res)
        !! Computes the sum \(\mathbb{res} = \mathbb{I} + \mathbb{A}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
            !! Standard 4th-order identity tensor (\(\delta_{ij}\delta_{kl}\)).
        type(ten_3D4O3sym), intent(in) :: a
            !! Fully symmetric 4th-order tensor (21 components).
        type(ten_3D4O3sym) :: res
            !! Resulting 4th-order tensor.
        
        res%vals = a%vals
        ! Apply identity to normal-normal interaction indices
        res%vals(1:3) = res%vals(1:3) + 1.0D0
        res%vals(7:8) = res%vals(7:8) + 1.0D0
        res%vals(12)  = res%vals(12)  + 1.0D0
    end function sum_I4O3T_3D4O3sym

    pure function sum_3D4O3sym_I4O3T(a, I4) result(res)
        !! Computes the sum \(\mathbb{res} = \mathbb{A} + \mathbb{I}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
        res%vals(7:8) = res%vals(7:8) + 1.0D0
        res%vals(12)  = res%vals(12)  + 1.0D0
    end function sum_3D4O3sym_I4O3T

    pure function sum_I4O4T_3D4O2sym(I4, a) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}^S + \mathbf{A}\).
        implicit none
        type(iden_4O4T),   intent(in) :: I4
            !! Standard fourth-order symmetric identity tensor.
        type(ten_3D4O2sym), intent(in) :: a
            !! Fourth-order tensor with minor symmetries.
        type(ten_3D4O2sym)             :: res

        res%vals = a%vals
        ! Update Voigt diagonal components
        res%vals(1,1) = res%vals(1,1) + 1.0D0
        res%vals(2,2) = res%vals(2,2) + 1.0D0
        res%vals(3,3) = res%vals(3,3) + 1.0D0
        res%vals(4,4) = res%vals(4,4) + 0.5D0
        res%vals(5,5) = res%vals(5,5) + 0.5D0
        res%vals(6,6) = res%vals(6,6) + 0.5D0
    end function sum_I4O4T_3D4O2sym

    pure function sum_3D4O2sym_I4O4T(a, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} + \mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
        type(iden_4O4T),   intent(in) :: I4
        type(ten_3D4O2sym)             :: res
        
        ! Commutative property: delegate to the other implementation
        res = sum_I4O4T_3D4O2sym(I4, a)
    end function sum_3D4O2sym_I4O4T

    pure function sum_I4O4T_3D4O3sym(I4, a) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}^S + \mathbf{A}\).
        implicit none
        type(iden_4O4T),    intent(in) :: I4
            !! Standard fourth-order symmetric identity tensor.
        type(ten_3D4O3sym), intent(in) :: a
            !! Fully symmetric fourth-order tensor (21 components).
        type(ten_3D4O3sym)             :: res

        res%vals = a%vals
        ! Add 1.0 to normal diagonals (11, 22, 33)
        res%vals(1:3) = res%vals(1:3) + 1.0D0
        ! Add 0.5 to shear diagonals (44, 55, 66 in Voigt-notation indices)
        res%vals(4:6) = res%vals(4:6) + 0.5D0
    end function sum_I4O4T_3D4O3sym

    pure function sum_3D4O3sym_I4O4T(a, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} + \mathbb{I}^S\).
        implicit none
        type(ten_3D4O3sym), intent(in) :: a
        type(iden_4O4T),    intent(in) :: I4
        type(ten_3D4O3sym)             :: res
        
        ! Delegate to ensure identical behavior and enable compiler inlining
        res = sum_I4O4T_3D4O3sym(I4, a)
    end function sum_3D4O3sym_I4O4T

    ! --- Identity + Identity ---
    !*************************************************************************

    pure function sum_I4O3TS_I4O4T(I3S, I4) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}_{3T} + \mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I3S
            !! Scaled Type-3 identity (\(c \cdot \delta_{ij}\delta_{kl}\)).
        type(iden_4O4T), intent(in) :: I4
            !! Standard Type-4 identity (Symmetric \(\mathbf{I}^S\)).
        type(ten_3D4O3sym) :: res
            !! Resulting 21-component symmetric tensor.
        
        res%vals = (/1.0D0 + I3S%val, 1.0D0 + I3S%val, 1.0D0 + I3S%val, & ! Diagonals
                     0.5D0, 0.5D0, 0.5D0,                               & ! Shear diagonals
                     I3S%val, I3S%val, 0.0D0, 0.0D0, 0.0D0, I3S%val,    & ! Coupling
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          & ! Others
                     0.0D0, 0.0D0, 0.0D0/)
    end function sum_I4O3TS_I4O4T

    pure function sum_I4O4T_I4O3TS(I4, I3S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{4T} + c\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4T), intent(in) :: I4
        type(iden_4O3TS), intent(in) :: I3S
        type(ten_3D4O3sym) :: res
        res = sum_I4O3TS_I4O4T(I3S, I4)
    end function sum_I4O4T_I4O3TS

    pure function sum_I4O3TS_I4O4TS(I3S, I4S) result(res)
        !! Computes \(\mathbb{res} = c_1\mathbb{I}_{3T} + c_2\mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I3S
            !! Scaled Type-3 identity (\(c_1 \cdot \delta_{ij}\delta_{kl}\)).
        type(iden_4O4TS), intent(in) :: I4S
            !! Scaled Type-4 identity (\(c_2 \cdot \mathbf{I}^S\)).
        type(ten_3D4O3sym) :: res
            !! Resulting 21-component symmetric tensor.
        
        res%vals = (/I4S%val + I3S%val, I4S%val + I3S%val, I4S%val + I3S%val, & ! Diagonals
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val,       & ! Shear diagonals
                     I3S%val, I3S%val, 0.0D0, 0.0D0, 0.0D0, I3S%val,          & ! Normal coupling
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,                & ! Others
                     0.0D0, 0.0D0, 0.0D0/)
    end function sum_I4O3TS_I4O4TS

    pure function sum_I4O4TS_I4O3TS(I4S, I3S) result(res)
        !! Computes \(\mathbb{res} = c_2\mathbb{I}_{4T} + c_1\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
        type(iden_4O3TS), intent(in) :: I3S
        type(ten_3D4O3sym) :: res
        res = sum_I4O3TS_I4O4TS(I3S, I4S)
    end function sum_I4O4TS_I4O3TS

    pure function sum_I4O3T_I4O4T(I3, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{3T} + \mathbb{I}_{4T}\).
        !! Result indices in Voigt: Diagonals = 2.0, Shears = 0.5, Coupling = 1.0.
        implicit none
        type(iden_4O3T), intent(in) :: I3
            !! Standard Type-3 identity (\(\delta_{ij}\delta_{kl}\)).
        type(iden_4O4T), intent(in) :: I4
            !! Standard Type-4 identity (Symmetric \(\mathbf{I}^S\)).
        type(ten_3D4O3sym) :: res
            !! Resulting 21-component symmetric tensor.
        
        res%vals = (/2.0D0, 2.0D0, 2.0D0, 0.5D0, 0.5D0, 0.5D0, &
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                     0.0D0, 0.0D0, 0.0D0/)
        ! Mapping: (1-3: diag, 4-6: shear, 7-8: coupling, 12: coupling, others: 0)
    end function sum_I4O3T_I4O4T

    pure function sum_I4O4T_I4O3T(I4, I3) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{4T} + \mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4T), intent(in) :: I4
        type(iden_4O3T), intent(in) :: I3
        type(ten_3D4O3sym) :: res
        res%vals = (/2.0D0, 2.0D0, 2.0D0, 0.5D0, 0.5D0, 0.5D0, &
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sum_I4O4T_I4O3T

        pure function sum_I4O3T_I4O4TS(I3, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{3T} + c\mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3T), intent(in) :: I3
            !! Standard Type-3 identity (\(\delta_{ij}\delta_{kl}\)).
        type(iden_4O4TS), intent(in) :: I4S
            !! Scaled Type-4 identity (\(c\mathbf{I}^S\)).
        type(ten_3D4O3sym) :: res
            !! Resulting 21-component symmetric tensor.
        
        res%vals = (/I4S%val + 1.0D0, I4S%val + 1.0D0, I4S%val + 1.0D0, & ! Diagonals
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val, & ! Shears
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0,          & ! Normal coupling
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          & ! Mixed
                     0.0D0, 0.0D0, 0.0D0/)                                ! Mixed
    end function sum_I4O3T_I4O4TS

    pure function sum_I4O4TS_I4O3T(I4S, I3) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}_{4T} + \mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
        type(iden_4O3T), intent(in) :: I3
        type(ten_3D4O3sym) :: res
        res%vals = (/I4S%val + 1.0D0, I4S%val + 1.0D0, I4S%val + 1.0D0, &
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val, &
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0,          &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sum_I4O4TS_I4O3T

    ! --- 2nd Order Identity + 2nd Order Tensor ---    
    !*************************************************************************

    pure function sum_3D2O_3D2Osym(a, b) result(res)
        !! Computes the addition of a general and a symmetric 2nd-order tensor.
        !!
        !! Mathematically: \( res_{ij} = A_{ij} + B_{ij} \)
        !! The symmetric tensor components are expanded to match the column-major 
        !! layout of the general tensor.
        implicit none
        type(ten_3D2O), intent(in) :: a
            !! General second-order tensor \(\mathbf{A}\)
        type(ten_3D2Osym), intent(in) :: b
            !! Symmetric second-order tensor \(\mathbf{B}\)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\)
        
        res%vals(1) = a%vals(1) + b%vals(1) ! xx
        res%vals(2) = a%vals(2) + b%vals(4) ! yx (sym: xy)
        res%vals(3) = a%vals(3) + b%vals(6) ! zx (sym: xz)
        res%vals(4) = a%vals(4) + b%vals(4) ! xy
        res%vals(5) = a%vals(5) + b%vals(2) ! yy
        res%vals(6) = a%vals(6) + b%vals(5) ! zy (sym: yz)
        res%vals(7) = a%vals(7) + b%vals(6) ! xz
        res%vals(8) = a%vals(8) + b%vals(5) ! yz
        res%vals(9) = a%vals(9) + b%vals(3) ! zz
    end function sum_3D2O_3D2Osym

    pure function sum_3D2Osym_3D2O(a, b) result(res)
        !! Computes the addition of a symmetric and a general 2nd-order tensor.
        !!
        !! Mathematically: \( res_{ij} = A_{ij} + B_{ij} \)
        !! The symmetric tensor components are expanded to match the column-major 
        !! layout of the general tensor.
        implicit none
        type(ten_3D2Osym), intent(in) :: a
            !! Symmetric second-order tensor \(\mathbf{A}\)
        type(ten_3D2O), intent(in) :: b
            !! General second-order tensor \(\mathbf{B}\)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\)
        
        res%vals(1) = a%vals(1) + b%vals(1) ! xx
        res%vals(2) = a%vals(4) + b%vals(2) ! yx (sym: xy)
        res%vals(3) = a%vals(6) + b%vals(3) ! zx (sym: xz)
        res%vals(4) = a%vals(4) + b%vals(4) ! xy
        res%vals(5) = a%vals(2) + b%vals(5) ! yy
        res%vals(6) = a%vals(5) + b%vals(6) ! zy (sym: yz)
        res%vals(7) = a%vals(6) + b%vals(7) ! xz
        res%vals(8) = a%vals(5) + b%vals(8) ! yz
        res%vals(9) = a%vals(3) + b%vals(9) ! zz
    end function sum_3D2Osym_3D2O


    ! =========================================================================
    ! IMPLEMENTATIONS: SUBTRACTION (-)
    ! =========================================================================
    
    ! --- 2nd Order Identity - 2nd Order Tensor ---    
    !*************************************************************************

    pure function sub_I2O_3D2Osym(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = -a%vals
        res%vals(1:3) = 1.0D0 + res%vals(1:3)
    end function sub_I2O_3D2Osym

    pure function sub_3D2Osym_I2O(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - 1.0D0
    end function sub_3D2Osym_I2O

    pure function sub_I2O_2D2Osym(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = -a%vals
        res%vals(1:3) = 1.0D0 + res%vals(1:3)
    end function sub_I2O_2D2Osym

    pure function sub_2D2Osym_I2O(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - 1.0D0
    end function sub_2D2Osym_I2O

    pure function sub_I2O_3D2O(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = -a%vals
        res%vals(1) = 1.0D0 + res%vals(1)
        res%vals(5) = 1.0D0 + res%vals(5)
        res%vals(9) = 1.0D0 + res%vals(9)
    end function sub_I2O_3D2O

    pure function sub_3D2O_I2O(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) - 1.0D0
        res%vals(5) = res%vals(5) - 1.0D0
        res%vals(9) = res%vals(9) - 1.0D0
    end function sub_3D2O_I2O

    ! --- Scaled 2nd Order Identity - 2nd Order Tensor ---
    !*************************************************************************

    pure function sub_I2OS_3D2Osym(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = c\mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = -a%vals
        res%vals(1:3) = I2%val + res%vals(1:3)
    end function sub_I2OS_3D2Osym

    pure function sub_3D2Osym_I2OS(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - I2%val
    end function sub_3D2Osym_I2OS

    pure function sub_I2OS_2D2Osym(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = c\mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = -a%vals
        res%vals(1:3) = I2%val + res%vals(1:3)
    end function sub_I2OS_2D2Osym

    pure function sub_2D2Osym_I2OS(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - I2%val
    end function sub_2D2Osym_I2OS

    pure function sub_I2OS_3D2O(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = c\mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = -a%vals
        res%vals(1) = I2%val + res%vals(1)
        res%vals(5) = I2%val + res%vals(5)
        res%vals(9) = I2%val + res%vals(9)
    end function sub_I2OS_3D2O

    pure function sub_3D2O_I2OS(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) - I2%val
        res%vals(5) = res%vals(5) - I2%val
        res%vals(9) = res%vals(9) - I2%val
    end function sub_3D2O_I2OS

    ! --- 4th Order Identity - 4th Order Tensor ---
    !*************************************************************************

    pure function sub_I4O3TS_3D4O2sym(I4S, a) result(res)
        !! Computes the subtraction \(\mathbb{res} = c\mathbb{I}_{3T} - \mathbb{A}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I4S
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = -a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + I4S%val
    end function sub_I4O3TS_3D4O2sym

    pure function sub_3D4O2sym_I4O3TS(a, I4S) result(res)
        !! Computes the subtraction \(\mathbb{res} = \mathbb{A} - c\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I4S
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) - I4S%val
    end function sub_3D4O2sym_I4O3TS

    pure function sub_I4O3TS_3D4O3sym(I4S, a) result(res)
        !! Computes the subtraction \(\mathbb{res} = c\mathbb{I}_{3T} - \mathbb{A}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I4S
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        
        res%vals = -a%vals
        res%vals(1:3) = res%vals(1:3) + I4S%val
        res%vals(7:8) = res%vals(7:8) + I4S%val
        res%vals(12)  = res%vals(12)  + I4S%val
    end function sub_I4O3TS_3D4O3sym

    pure function sub_3D4O3sym_I4O3TS(a, I4S) result(res)
        !! Computes the subtraction \(\mathbb{res} = \mathbb{A} - c\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I4S
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - I4S%val
        res%vals(7:8) = res%vals(7:8) - I4S%val
        res%vals(12)  = res%vals(12)  - I4S%val
    end function sub_3D4O3sym_I4O3TS

        pure function sub_I4O4TS_3D4O2sym(I4S, a) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}^S - \mathbb{A}\).
        implicit none
        type(iden_4O4TS),  intent(in) :: I4S
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = -a%vals
        res%vals(1,1) = res%vals(1,1) + c
        res%vals(2,2) = res%vals(2,2) + c
        res%vals(3,3) = res%vals(3,3) + c
        res%vals(4,4) = res%vals(4,4) + half_c
        res%vals(5,5) = res%vals(5,5) + half_c
        res%vals(6,6) = res%vals(6,6) + half_c
    end function sub_I4O4TS_3D4O2sym

    pure function sub_3D4O2sym_I4O4TS(a, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{A} - c\mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
        type(iden_4O4TS),  intent(in) :: I4S
        type(ten_3D4O2sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = a%vals
        res%vals(1,1) = res%vals(1,1) - c
        res%vals(2,2) = res%vals(2,2) - c
        res%vals(3,3) = res%vals(3,3) - c
        res%vals(4,4) = res%vals(4,4) - half_c
        res%vals(5,5) = res%vals(5,5) - half_c
        res%vals(6,6) = res%vals(6,6) - half_c
    end function sub_3D4O2sym_I4O4TS

    pure function sub_I4O4TS_3D4O3sym(I4S, a) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}^S - \mathbf{A}\).
        implicit none
        type(iden_4O4TS),   intent(in) :: I4S
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = -a%vals
        ! Add scaled identity components to the negated tensor
        res%vals(1:3) = res%vals(1:3) + c
        res%vals(4:6) = res%vals(4:6) + half_c
    end function sub_I4O4TS_3D4O3sym

    pure function sub_3D4O3sym_I4O4TS(a, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} - c\mathbb{I}^S\).
        implicit none
        type(ten_3D4O3sym), intent(in) :: a
        type(iden_4O4TS),   intent(in) :: I4S
        type(ten_3D4O3sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = a%vals
        ! Subtract scaled identity components from the diagonal
        res%vals(1:3) = res%vals(1:3) - c
        res%vals(4:6) = res%vals(4:6) - half_c
    end function sub_3D4O3sym_I4O4TS

    pure function sub_I4O3T_3D4O2sym(I4, a) result(res)
        !! Computes the subtraction \(\mathbb{res} = \mathbb{I} - \mathbb{A}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = -a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1.0D0
    end function sub_I4O3T_3D4O2sym

    pure function sub_3D4O2sym_I4O3T(a, I4) result(res)
        !! Computes the subtraction \(\mathbb{res} = \mathbb{A} - \mathbb{I}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) - 1.0D0
    end function sub_3D4O2sym_I4O3T

    pure function sub_I4O3T_3D4O3sym(I4, a) result(res)
        !! Computes the subtraction \(\mathbb{res} = \mathbb{I} - \mathbb{A}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        
        res%vals = -a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
        res%vals(7:8) = res%vals(7:8) + 1.0D0
        res%vals(12)  = res%vals(12)  + 1.0D0
    end function sub_I4O3T_3D4O3sym

    pure function sub_3D4O3sym_I4O3T(a, I4) result(res)
        !! Computes the subtraction \(\mathbb{res} = \mathbb{A} - \mathbb{I}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - 1.0D0
        res%vals(7:8) = res%vals(7:8) - 1.0D0
        res%vals(12)  = res%vals(12)  - 1.0D0
    end function sub_3D4O3sym_I4O3T

    pure function sub_I4O4T_3D4O2sym(I4, a) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}^S - \mathbf{A}\).
        implicit none
        type(iden_4O4T),   intent(in) :: I4
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym)             :: res

        res%vals = -a%vals
        res%vals(1,1) = res%vals(1,1) + 1.0D0
        res%vals(2,2) = res%vals(2,2) + 1.0D0
        res%vals(3,3) = res%vals(3,3) + 1.0D0
        res%vals(4,4) = res%vals(4,4) + 0.5D0
        res%vals(5,5) = res%vals(5,5) + 0.5D0
        res%vals(6,6) = res%vals(6,6) + 0.5D0
    end function sub_I4O4T_3D4O2sym

    pure function sub_3D4O2sym_I4O4T(a, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} - \mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
        type(iden_4O4T),   intent(in) :: I4
        type(ten_3D4O2sym)             :: res

        res%vals = a%vals
        res%vals(1,1) = res%vals(1,1) - 1.0D0
        res%vals(2,2) = res%vals(2,2) - 1.0D0
        res%vals(3,3) = res%vals(3,3) - 1.0D0
        res%vals(4,4) = res%vals(4,4) - 0.5D0
        res%vals(5,5) = res%vals(5,5) - 0.5D0
        res%vals(6,6) = res%vals(6,6) - 0.5D0
    end function sub_3D4O2sym_I4O4T

    pure function sub_I4O4T_3D4O3sym(I4, a) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}^S - \mathbf{A}\).
        implicit none
        type(iden_4O4T),    intent(in) :: I4
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym)             :: res

        ! Compute -A
        res%vals = -a%vals
        ! Add Identity components: res = I + (-A)
        res%vals(1:3) = res%vals(1:3) + 1.0D0
        res%vals(4:6) = res%vals(4:6) + 0.5D0
    end function sub_I4O4T_3D4O3sym

    pure function sub_3D4O3sym_I4O4T(a, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} - \mathbb{I}^S\).
        implicit none
        type(ten_3D4O3sym), intent(in) :: a
        type(iden_4O4T),    intent(in) :: I4
        type(ten_3D4O3sym)             :: res

        res%vals = a%vals
        ! Subtract identity components from diagonal
        res%vals(1:3) = res%vals(1:3) - 1.0D0
        res%vals(4:6) = res%vals(4:6) - 0.5D0
    end function sub_3D4O3sym_I4O4T

    ! --- Identity - Identity ---
    !*************************************************************************

    pure function sub_I4O3TS_I4O4T(I3S, I4) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}_{3T} - \mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I3S
        type(iden_4O4T), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/I3S%val - 1.0D0, I3S%val - 1.0D0, I3S%val - 1.0D0, &
                     -0.5D0, -0.5D0, -0.5D0,                            &
                     I3S%val, I3S%val, 0.0D0, 0.0D0, 0.0D0, I3S%val,    &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O3TS_I4O4T

    pure function sub_I4O4T_I4O3TS(I4, I3S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{4T} - c\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4T), intent(in) :: I4
        type(iden_4O3TS), intent(in) :: I3S
        type(ten_3D4O3sym) :: res
        res%vals = (/1.0D0 - I3S%val, 1.0D0 - I3S%val, 1.0D0 - I3S%val, &
                     0.5D0, 0.5D0, 0.5D0,                               &
                     -I3S%val, -I3S%val, 0.0D0, 0.0D0, 0.0D0, -I3S%val, &
                      0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,         &
                      0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O4T_I4O3TS

    pure function sub_I4O3TS_I4O4TS(I3S, I4S) result(res)
        !! Computes \(\mathbb{res} = c_1\mathbb{I}_{3T} - c_2\mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I3S
        type(iden_4O4TS), intent(in) :: I4S
        type(ten_3D4O3sym) :: res
        res%vals = (/I3S%val - I4S%val, I3S%val - I4S%val, I3S%val - I4S%val, &
                     -0.5D0 * I4S%val, -0.5D0 * I4S%val, -0.5D0 * I4S%val,    &
                     I3S%val, I3S%val, 0.0D0, 0.0D0, 0.0D0, I3S%val,          &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,                &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O3TS_I4O4TS

    pure function sub_I4O4TS_I4O3TS(I4S, I3S) result(res)
        !! Computes \(\mathbb{res} = c_2\mathbb{I}_{4T} - c_1\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
        type(iden_4O3TS), intent(in) :: I3S
        type(ten_3D4O3sym) :: res
        res%vals = (/I4S%val - I3S%val, I4S%val - I3S%val, I4S%val - I3S%val, &
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val,       &
                     -I3S%val, -I3S%val, 0.0D0, 0.0D0, 0.0D0, -I3S%val,       &
                      0.0D0,  0.0D0, 0.0D0, 0.0D0, 0.0D0,  0.0D0,             &
                      0.0D0,  0.0D0, 0.0D0/)
    end function sub_I4O4TS_I4O3TS

    pure function sub_I4O3T_I4O4T(I3, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{3T} - \mathbb{I}_{4T}\).
        !! Result indices in Voigt: Diagonals = 0.0, Shears = -0.5, Coupling = 1.0.
        implicit none
        type(iden_4O3T), intent(in) :: I3
        type(iden_4O4T), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/0.0D0, 0.0D0, 0.0D0, -0.5D0, -0.5D0, -0.5D0, &
                     1.0D0, 1.0D0, 0.0D0,  0.0D0,  0.0D0,  1.0D0, &
                     0.0D0, 0.0D0, 0.0D0,  0.0D0,  0.0D0,  0.0D0, &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O3T_I4O4T

    pure function sub_I4O4T_I4O3T(I4, I3) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{4T} - \mathbb{I}_{3T}\).
        !! Result indices in Voigt: Diagonals = 0.0, Shears = 0.5, Coupling = -1.0.
        implicit none
        type(iden_4O4T), intent(in) :: I4
        type(iden_4O3T), intent(in) :: I3
        type(ten_3D4O3sym) :: res
        res%vals = (/ 0.0D0,  0.0D0, 0.0D0, 0.5D0, 0.5D0,  0.5D0, &
                     -1.0D0, -1.0D0, 0.0D0, 0.0D0, 0.0D0, -1.0D0, &
                      0.0D0,  0.0D0, 0.0D0, 0.0D0, 0.0D0,  0.0D0, &
                      0.0D0,  0.0D0, 0.0D0/)
    end function sub_I4O4T_I4O3T

    pure function sub_I4O3T_I4O4TS(I3, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{3T} - c\mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3T), intent(in) :: I3
        type(iden_4O4TS), intent(in) :: I4S
        type(ten_3D4O3sym) :: res
        res%vals = (/1.0D0 - I4S%val, 1.0D0 - I4S%val, 1.0D0 - I4S%val, &
                     -0.5D0 * I4S%val, -0.5D0 * I4S%val, -0.5D0 * I4S%val, &
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0,          &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O3T_I4O4TS

    pure function sub_I4O4TS_I4O3T(I4S, I3) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}_{4T} - \mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
        type(iden_4O3T), intent(in) :: I3
        type(ten_3D4O3sym) :: res
        res%vals = (/I4S%val - 1.0D0, I4S%val - 1.0D0, I4S%val - 1.0D0, &
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val, &
                     -1.0D0, -1.0D0, 0.0D0, 0.0D0, 0.0D0, -1.0D0,       &
                      0.0D0,  0.0D0, 0.0D0, 0.0D0, 0.0D0,  0.0D0,       &
                      0.0D0,  0.0D0, 0.0D0/)
    end function sub_I4O4TS_I4O3T

    ! --- 2nd Order Identity + 2nd Order Tensor ---    
    !*************************************************************************

    pure function sub_3D2O_3D2Osym(a, b) result(res)
        !! Computes the subtraction of a symmetric tensor from a general tensor.
        !!
        !! Mathematically: \( res_{ij} = A_{ij} - B_{ij} \)
        implicit none
        type(ten_3D2O), intent(in) :: a     
            !! General second-order tensor \(\mathbf{A}\)
        type(ten_3D2Osym), intent(in) :: b  
            !! Symmetric second-order tensor \(\mathbf{B}\)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\)
        
        res%vals(1) = a%vals(1) - b%vals(1) ! xx
        res%vals(2) = a%vals(2) - b%vals(4) ! yx - xy
        res%vals(3) = a%vals(3) - b%vals(6) ! zx - xz
        res%vals(4) = a%vals(4) - b%vals(4) ! xy
        res%vals(5) = a%vals(5) - b%vals(2) ! yy
        res%vals(6) = a%vals(6) - b%vals(5) ! zy - yz
        res%vals(7) = a%vals(7) - b%vals(6) ! xz
        res%vals(8) = a%vals(8) - b%vals(5) ! yz
        res%vals(9) = a%vals(9) - b%vals(3) ! zz
    end function sub_3D2O_3D2Osym

    pure function sub_3D2Osym_3D2O(a, b) result(res)
        !! Computes the subtraction of a general tensor from a symmetric tensor.
        !!
        !! Mathematically: \( res_{ij} = A_{ij} - B_{ij} \)
        implicit none
        type(ten_3D2Osym), intent(in) :: a  
            !! Symmetric second-order tensor \(\mathbf{A}\)
        type(ten_3D2O), intent(in) :: b     
            !! General second-order tensor \(\mathbf{B}\)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\)
        
        res%vals(1) = a%vals(1) - b%vals(1) ! xx
        res%vals(2) = a%vals(4) - b%vals(2) ! yx (sym: xy)
        res%vals(3) = a%vals(6) - b%vals(3) ! zx (sym: xz)
        res%vals(4) = a%vals(4) - b%vals(4) ! xy
        res%vals(5) = a%vals(2) - b%vals(5) ! yy
        res%vals(6) = a%vals(5) - b%vals(6) ! zy (sym: yz)
        res%vals(7) = a%vals(6) - b%vals(7) ! xz
        res%vals(8) = a%vals(5) - b%vals(8) ! yz
        res%vals(9) = a%vals(3) - b%vals(9) ! zz
    end function sub_3D2Osym_3D2O

end module muscle_tensor_ops_addition_subtraction