module mod_ops_contraction_double
    !! This module defines the overloaded double contraction operator `.ddot.` (tensor inner product)
    !! for mixed-type interactions in the 3D and 2D tensor engine.
    !!
    !! Mathematically, the double contraction of two tensors \(\mathbf{A}\) and \(\mathbf{B}\) 
    !! is represented as:
    !! \[ \alpha = \mathbf{A} : \mathbf{B} = A_{ij} B_{ij} \]
    !! or for a fourth-order tensor \(\mathbb{C}\) and a second-order tensor \(\mathbf{A}\):
    !! \[ \mathbf{res} = \mathbb{C} : \mathbf{A} \quad \implies \quad res_{ij} = C_{ijkl} A_{kl} \]
    !!
    !! This module centralizes these operations to eliminate interface redundancy.
    
    use, intrinsic :: iso_fortran_env, only : real64
    
    ! Import necessary types
    use mod_iden_2O
    use mod_iden_2OS
    use mod_iden_4O4T
    use mod_iden_4O4TS
    use mod_ten_2D2Osym
    use mod_ten_2D4O3sym
    use mod_ten_3D2O
    use mod_ten_3D2Osym
    use mod_ten_3D4O2sym
    use mod_ten_3D4O3sym

    implicit none
    private

    ! =========================================================================
    ! PUBLIC INTERFACES
    ! =========================================================================
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        ! --- 2nd Order Identity : 2nd Order Tensor ---
        module procedure ddot_I2O_2D2Osym
        module procedure ddot_2D2Osym_I2O
        module procedure ddot_I2O_3D2O
        module procedure ddot_3D2O_I2O
        module procedure ddot_I2O_3D2Osym
        module procedure ddot_3D2Osym_I2O

        ! --- Scaled 2nd Order Identity : 2nd Order Tensor ---
        module procedure ddot_I2OS_2D2Osym
        module procedure ddot_2D2Osym_I2OS
        module procedure ddot_I2OS_3D2O
        module procedure ddot_3D2O_I2OS
        module procedure ddot_I2OS_3D2Osym
        module procedure ddot_3D2Osym_I2OS

        ! --- 2nd Order Tensor : 2nd Order Tensor ---
        module procedure ddot_3D2O_3D2Osym
        module procedure ddot_3D2Osym_3D2O

        ! --- 4th Order Symmetric Identity : 2nd Order Tensor ---
        module procedure ddot_I4O4T_3D2O
        module procedure ddot_3D2O_I4O4T

        ! --- Scaled 4th Order Symmetric Identity : 2nd Order Tensor ---
        module procedure ddot_I4O4TS_3D2O
        module procedure ddot_3D2O_I4O4TS
        module procedure ddot_I4O4TS_3D2Osym
        module procedure ddot_3D2Osym_I4O4TS

        ! --- 4th Order Tensor : 2nd Order Tensor ---
        module procedure ddot_2D4O3sym_2D2Osym
        module procedure ddot_2D2Osym_2D4O3sym
        module procedure ddot_3D4O2sym_3D2Osym
        module procedure ddot_3D2Osym_3D4O2sym
        module procedure ddot_3D4O3sym_3D2Osym
        module procedure ddot_3D2Osym_3D4O3sym

        ! --- 4th Order Tensor : 4th Order Tensor ---
        module procedure ddot_3D4O3sym_3D4O3sym
    end interface

contains

    ! =========================================================================
    ! 1. IDENTITIES (2ND ORDER) : TENSORS (2ND ORDER)
    ! =========================================================================

    pure function ddot_I2O_2D2Osym(I2, b) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{I} : \mathbf{B}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: b
            !! The symmetric 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! The resulting scalar (trace of \(\mathbf{B}\)).
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_I2O_2D2Osym

    pure function ddot_2D2Osym_I2O(b, I2) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{B} : \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: b
            !! The symmetric 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! The resulting scalar (trace of \(\mathbf{B}\)).
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_2D2Osym_I2O

    pure function ddot_I2O_3D2O(I2, b) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{I} : \mathbf{B}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(5) + b%vals(9)
    end function ddot_I2O_3D2O

    pure function ddot_3D2O_I2O(b, I2) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{B} : \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(5) + b%vals(9)
    end function ddot_3D2O_I2O

    pure function ddot_I2O_3D2Osym(I2, b) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{I} : \mathbf{B}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_I2O_3D2Osym

    pure function ddot_3D2Osym_I2O(b, I2) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{B} : \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_3D2Osym_I2O

    ! =========================================================================
    ! 2. SCALED IDENTITIES (2ND ORDER) : TENSORS (2ND ORDER)
    ! =========================================================================

    pure function ddot_I2OS_2D2Osym(I2, b) result(res)
        !! Computes the scaled double contraction \(c\mathbf{I} : \mathbf{B} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
            !! Scaled identity tensor \(c\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: b
            !! Symmetric 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! Scalar result.
        res = I2%val * (b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_I2OS_2D2Osym

    pure function ddot_2D2Osym_I2OS(b, I2) result(res)
        !! Computes the scaled double contraction \(\mathbf{B} : c\mathbf{I} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: b
        real(real64) :: res
        res = I2%val * (b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_2D2Osym_I2OS

    pure function ddot_I2OS_3D2O(I2, b) result(res)
        !! Computes the scaled double contraction \(c\mathbf{I} : \mathbf{B} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
            !! Scaled identity tensor \(c\mathbf{I}\).
        type(ten_3D2O), intent(in) :: b
            !! General 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! Scalar result.
        res = I2%val * (b%vals(1) + b%vals(5) + b%vals(9))
    end function ddot_I2OS_3D2O

    pure function ddot_3D2O_I2OS(b, I2) result(res)
        !! Computes the scaled double contraction \(\mathbf{B} : c\mathbf{I} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: b
        real(real64) :: res
        res = I2%val * (b%vals(1) + b%vals(5) + b%vals(9))
    end function ddot_3D2O_I2OS

    pure function ddot_I2OS_3D2Osym(I2, b) result(res)
        !! Computes the scaled double contraction \(c\mathbf{I} : \mathbf{B} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
            !! Scaled identity tensor \(c\mathbf{I}\).
        type(ten_3D2Osym), intent(in) :: b
            !! Symmetric 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! Scalar result.
        res = I2%val * (b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_I2OS_3D2Osym

    pure function ddot_3D2Osym_I2OS(b, I2) result(res)
        !! Computes the scaled double contraction \(\mathbf{B} : c\mathbf{I} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = I2%val * (b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_3D2Osym_I2OS

    ! =========================================================================
    ! 3. TENSORS (2ND ORDER) : TENSORS (2ND ORDER)
    ! =========================================================================

    pure function ddot_3D2O_3D2Osym(a, b) result(res)
        !! Computes the double contraction between a general and a symmetric tensor.
        !!
        !! Mathematically: \( \alpha = A_{ij} B_{ij} \)
        !! Result is a scalar `real(real64)`.
        !!
        !! Since B is symmetric (\(B_{ij} = B_{ji}\)), we can group terms:
        !! \(\mathbf{A}:\mathbf{B} = A_{11} B_{11} + A_{22} B_{22} + A_{33} B_{33} 
        !!     + (A_{12} + A_{21}) B_{12} + (A_{23} + A_{32}) B_{23} + (A_{13} + A_{31}) B_{13}\)
        implicit none
        type(ten_3D2O), intent(in) :: a     
            !! General second-order tensor \(\mathbf{A}\)
        type(ten_3D2Osym), intent(in) :: b  
            !! Symmetric second-order tensor \(\mathbf{B}\)
        real(real64) :: res
            !! Scalar result \(\alpha\)
        
        res = a%vals(1)*b%vals(1) + a%vals(5)*b%vals(2) + a%vals(9)*b%vals(3) &
            + (a%vals(4) + a%vals(2))*b%vals(4) &  ! (xy + yx) * xy_sym
            + (a%vals(8) + a%vals(6))*b%vals(5) &  ! (yz + zy) * yz_sym
            + (a%vals(7) + a%vals(3))*b%vals(6)    ! (xz + zx) * xz_sym
    end function ddot_3D2O_3D2Osym

    pure function ddot_3D2Osym_3D2O(a, b) result(res)
        !! Computes the double contraction between a symmetric and a general tensor.
        !!
        !! Mathematically: \( \alpha = A_{ij} B_{ij} \)
        !! Due to the commutative property of the inner product, this is identical
        !! to the computation in `ddot_3D2O_3D2Osym`.
        implicit none
        type(ten_3D2Osym), intent(in) :: a  
            !! Symmetric second-order tensor \(\mathbf{A}\)
        type(ten_3D2O), intent(in) :: b     
            !! General second-order tensor \(\mathbf{B}\)
        real(real64) :: res
            !! Scalar result \(\alpha\)
        
        res = a%vals(1)*b%vals(1) + a%vals(2)*b%vals(5) + a%vals(3)*b%vals(9) &
            + a%vals(4)*(b%vals(4) + b%vals(2)) &  ! xy_sym * (xy + yx)
            + a%vals(5)*(b%vals(8) + b%vals(6)) &  ! yz_sym * (yz + zy)
            + a%vals(6)*(b%vals(7) + b%vals(3))    ! xz_sym * (xz + zx)
    end function ddot_3D2Osym_3D2O

    ! =========================================================================
    ! 4. 4TH ORDER SYMMETRIC IDENTITIES : 2ND ORDER TENSORS
    ! =========================================================================

    pure function ddot_I4O4T_3D2O(I4, a) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbb{I}^S : \mathbf{A}\).
        !!
        !! This results in the symmetric part of tensor \(\mathbf{A}\).
        !!
        !! Mapping logic from 9-component (column-major) to 6-component (Voigt):
        !! - res(11) = A(11)
        !! - res(22) = A(22)
        !! - res(33) = A(33)
        !! - res(12) = 0.5 * (A(12) + A(21))
        !! - res(23) = 0.5 * (A(23) + A(32))
        !! - res(13) = 0.5 * (A(13) + A(31))
        implicit none
        type(iden_4O4T), intent(in) :: I4
            !! Standard fourth-order symmetric identity tensor.
        type(ten_3D2O),  intent(in) :: a
            !! General second-order tensor (9 components, column-major).
        type(ten_3D2Osym)           :: res
            !! Resulting symmetric second-order tensor (6 components, Voigt).

        ! Diagonal components
        res%vals(1) = a%vals(1) ! xx
        res%vals(2) = a%vals(5) ! yy
        res%vals(3) = a%vals(9) ! zz
        
        ! Symmetric off-diagonal components: 0.5 * (A_ij + A_ji)
        ! Note: ten_3D2O storage is (11, 21, 31, 12, 22, 32, 13, 23, 33)
        res%vals(4) = 0.5D0 * (a%vals(2) + a%vals(4)) ! xy part: (21 + 12)/2
        res%vals(5) = 0.5D0 * (a%vals(6) + a%vals(8)) ! yz part: (32 + 23)/2
        res%vals(6) = 0.5D0 * (a%vals(3) + a%vals(7)) ! xz part: (31 + 13)/2
        
    end function ddot_I4O4T_3D2O

    pure function ddot_3D2O_I4O4T(a, I4) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbf{A} : \mathbb{I}^S\).
        !!
        !! Equivalent to \(\mathbb{I}^S : \mathbf{A}\) due to the symmetry of the operator.
        implicit none
        type(ten_3D2O),  intent(in) :: a
        type(iden_4O4T), intent(in) :: I4
        type(ten_3D2Osym)           :: res
        
        ! Delegate to the primary implementation
        res = ddot_I4O4T_3D2O(I4, a)
        
    end function ddot_3D2O_I4O4T

    ! =========================================================================
    ! 4. SCALED 4TH ORDER SYMMETRIC IDENTITIES : 2ND ORDER TENSORS
    ! =========================================================================

    pure function ddot_I4O4TS_3D2O(I4S, a) result(res)
        !! Computes the scaled double contraction product: \(\mathbf{res} = c\mathbb{I}^S : \mathbf{A}\).
        !!
        !! Mathematically, this yields the scaled symmetric part of tensor \(\mathbf{A}\).
        !!
        !! Mapping logic from 9-component (column-major) to 6-component (Voigt):
        !! - res(1) = c * A(11)
        !! - res(2) = c * A(22)
        !! - res(3) = c * A(33)
        !! - res(4) = 0.5 * c * (A(12) + A(21))
        !! - res(5) = 0.5 * c * (A(23) + A(32))
        !! - res(6) = 0.5 * c * (A(13) + A(31))
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
            !! Scaled fourth-order symmetric identity tensor (\(c\mathbb{I}^S\)).
        type(ten_3D2O),   intent(in) :: a
            !! General second-order tensor (9 components, column-major).
        type(ten_3D2Osym)           :: res
            !! Resulting symmetric second-order tensor (6 components, Voigt).
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        ! Diagonal components: c * A_ii
        res%vals(1) = a%vals(1) * c ! xx
        res%vals(2) = a%vals(5) * c ! yy
        res%vals(3) = a%vals(9) * c ! zz
        
        ! Symmetric off-diagonal components: 0.5 * c * (A_ij + A_ji)
        ! Note indices for ten_3D2O: (11:1, 21:2, 31:3, 12:4, 22:5, 32:6, 13:7, 23:8, 33:9)
        res%vals(4) = (a%vals(2) + a%vals(4)) * half_c ! xy part
        res%vals(5) = (a%vals(6) + a%vals(8)) * half_c ! yz part
        res%vals(6) = (a%vals(3) + a%vals(7)) * half_c ! xz part
        
    end function ddot_I4O4TS_3D2O

    pure function ddot_3D2O_I4O4TS(a, I4S) result(res)
        !! Computes the scaled double contraction product: \(\mathbf{res} = \mathbf{A} : c\mathbb{I}^S\).
        !!
        !! Equivalent to \(c\mathbb{I}^S : \mathbf{A}\) due to symmetry of the contraction.
        implicit none
        type(ten_3D2O),   intent(in) :: a
        type(iden_4O4TS), intent(in) :: I4S
        type(ten_3D2Osym)           :: res
        
        ! Delegate to the primary implementation to ensure consistency
        res = ddot_I4O4TS_3D2O(I4S, a)
        
    end function ddot_3D2O_I4O4TS

    pure function ddot_I4O4TS_3D2Osym(I4S, a) result(res)
        !! Computes the scaled double contraction product: \(\mathbf{res} = c\mathbb{I}^S : \mathbf{A}\).
        !!
        !! For a symmetric tensor \(\mathbf{A}\), this is equivalent to \(c \cdot \mathbf{A}\).
        implicit none
        type(iden_4O4TS),  intent(in) :: I4S
            !! Scaled fourth-order symmetric identity tensor (\(c\mathbb{I}^S\)).
        type(ten_3D2Osym), intent(in) :: a
            !! Symmetric second-order tensor (6 components, Voigt).
        type(ten_3D2Osym)             :: res
            !! Resulting symmetric second-order tensor.

        res%vals = I4S%val * a%vals
    end function ddot_I4O4TS_3D2Osym

    pure function ddot_3D2Osym_I4O4TS(a, I4S) result(res)
        !! Computes the scaled double contraction product: \(\mathbf{res} = \mathbf{A} : c\mathbb{I}^S\).
        !!
        !! Equivalent to \(c\mathbb{I}^S : \mathbf{A}\).
        implicit none
        type(ten_3D2Osym), intent(in) :: a
        type(iden_4O4TS),  intent(in) :: I4S
        type(ten_3D2Osym)             :: res
        
        ! Delegate to the primary implementation
        res = ddot_I4O4TS_3D2Osym(I4S, a)
    end function ddot_3D2Osym_I4O4TS

    ! =========================================================================
    ! 5. 4TH ORDER TENSORS : 2ND ORDER TENSORS
    ! =========================================================================

    pure function ddot_2D4O3sym_2D2Osym(a, b) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbb{A} : \mathbf{b}\).
        !!
        !! Mathematically: \( res_{ij} = A_{ijkl} b_{kl} \).
        !!
        !! Performance Optimization:
        !! -------------------------
        !! Performs manual loop unrolling of the contraction using Voigt notation 
        !! (4x4 matrix contracted with a 4x1 vector). Shear terms are implicitly 
        !! multiplied by 2 due to the nature of symmetric tensor contractions.
        !!
        !! Voigt Mapping Reference:
        !! ```
        !!  | ( 1:1111) ( 5:1122) ( 8:1133) (10:1112) |
        !!  | ( 5:2211) ( 2:2222) ( 6:2233) ( 9:2212) |
        !!  | ( 8:3311) ( 6:3322) ( 3:3333) ( 7:3312) |
        !!  | (10:1211) ( 9:1222) ( 7:1233) ( 4:1212) |
        !! ```
        implicit none
        class(ten_2D4O3sym), intent(in) :: a
            !! Fourth-order fully symmetric tensor \(\mathbb{A}\)
        class(ten_2D2Osym), intent(in) :: b
            !! Second-order symmetric tensor \(\mathbf{b}\)
        type(ten_2D2Osym) :: res
            !! Resulting second-order symmetric tensor \(\mathbf{res}\)
        
        real(real64) :: a1111b11, a1122b22, a1133b33, a1112b12
        real(real64) :: a2211b11, a2222b22, a2233b33, a2212b12
        real(real64) :: a3311b11, a3322b22, a3333b33, a3312b12
        real(real64) :: a1112b11, a2212b22, a3312b33, a1212b12

        a1111b11 = a%vals( 1)*b%vals(1)
        a1122b22 = a%vals( 5)*b%vals(2)
        a1133b33 = a%vals( 8)*b%vals(3)
        a1112b12 = a%vals(10)*b%vals(4)

        a2211b11 = a%vals( 5)*b%vals(1)
        a2222b22 = a%vals( 2)*b%vals(2)
        a2233b33 = a%vals( 6)*b%vals(3)
        a2212b12 = a%vals( 9)*b%vals(4)

        a3311b11 = a%vals( 8)*b%vals(1)
        a3322b22 = a%vals( 6)*b%vals(2)
        a3333b33 = a%vals( 3)*b%vals(3)
        a3312b12 = a%vals( 7)*b%vals(4)

        a1112b11 = a%vals(10)*b%vals(1)
        a2212b22 = a%vals( 9)*b%vals(2)
        a3312b33 = a%vals( 7)*b%vals(3)
        a1212b12 = a%vals( 4)*b%vals(4)
        
        res%vals(1) = a1111b11 + a1122b22 + a1133b33 + 2.0D0 * a1112b12
        res%vals(2) = a2211b11 + a2222b22 + a2233b33 + 2.0D0 * a2212b12
        res%vals(3) = a3311b11 + a3322b22 + a3333b33 + 2.0D0 * a3312b12
        res%vals(4) = a1112b11 + a2212b22 + a3312b33 + 2.0D0 * a1212b12
        
    end function ddot_2D4O3sym_2D2Osym

    pure function ddot_2D2Osym_2D4O3sym(b, a) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbf{b} : \mathbb{A}\).
        !!
        !! Mathematically: \( res_{ij} = b_{kl} A_{klij} \).
        !! Due to the major symmetry of \(\mathbb{A}\), this is equivalent to \(\mathbb{A} : \mathbf{b}\).
        !!
        !! Performance Optimization:
        !! -------------------------
        !! Performs manual loop unrolling of the contraction using Voigt notation 
        !! (1x4 vector contracted with a 4x4 matrix).
        !!
        !! Voigt Mapping Reference:
        !! ```
        !!  | ( 1:1111) ( 5:1122) ( 8:1133) (10:1112) |
        !!  | ( 5:2211) ( 2:2222) ( 6:2233) ( 9:2212) |
        !!  | ( 8:3311) ( 6:3322) ( 3:3333) ( 7:3312) |
        !!  | (10:1211) ( 9:1222) ( 7:1233) ( 4:1212) |
        !! ```
        implicit none
        class(ten_2D2Osym), intent(in) :: b
            !! Second-order symmetric tensor \(\mathbf{b}\)
        class(ten_2D4O3sym), intent(in) :: a
            !! Fourth-order fully symmetric tensor \(\mathbb{A}\)
        type(ten_2D2Osym) :: res
            !! Resulting second-order symmetric tensor \(\mathbf{res}\)
        
        real(real64) :: a1111b11, a2211b22, a3311b33, a1211b12
        real(real64) :: a1122b11, a2222b22, a3322b33, a1222b12
        real(real64) :: a1133b11, a2233b22, a3333b33, a1233b12
        real(real64) :: a1112b11, a2212b22, a3312b33, a1212b12

        a1111b11 = a%vals( 1)*b%vals(1)
        a2211b22 = a%vals( 5)*b%vals(2)
        a3311b33 = a%vals( 8)*b%vals(3)
        a1211b12 = a%vals(10)*b%vals(4)

        a1122b11 = a%vals( 5)*b%vals(1)
        a2222b22 = a%vals( 2)*b%vals(2)
        a3322b33 = a%vals( 6)*b%vals(3)
        a1222b12 = a%vals( 9)*b%vals(4)

        a1133b11 = a%vals( 8)*b%vals(1)
        a2233b22 = a%vals( 6)*b%vals(2)
        a3333b33 = a%vals( 3)*b%vals(3)
        a1233b12 = a%vals( 7)*b%vals(4)

        a1112b11 = a%vals(10)*b%vals(1)
        a2212b22 = a%vals( 9)*b%vals(2)
        a3312b33 = a%vals( 7)*b%vals(3)
        a1212b12 = a%vals( 4)*b%vals(4)
        
        res%vals(1) = a1111b11 + a2211b22 + a3311b33 + 2.0D0 * a1211b12
        res%vals(2) = a1122b11 + a2222b22 + a3322b33 + 2.0D0 * a1222b12
        res%vals(3) = a1133b11 + a2233b22 + a3333b33 + 2.0D0 * a1233b12
        res%vals(4) = a1112b11 + a2212b22 + a3312b33 + 2.0D0 * a1212b12
   
    end function ddot_2D2Osym_2D4O3sym

    pure function ddot_3D4O2sym_3D2Osym(a, b) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbb{A} : \mathbf{b}\).
        !!
        !! Mathematically: \(res_{ij} = A_{ijkl} b_{kl}\).
        !! In Voigt notation, this is a matrix-vector product: \([\text{res}]_I = [\mathbb{A}]_{IJ} [\mathbf{b}]_J\).
        !! Shear components of \(\mathbf{b}\) (indices 4, 5, 6) are weighted by a factor of 2
        !! in the contraction.
        !!
        !! This implementation is highly efficient as it performs a column-wise sum,
        !! which aligns with Fortran's column-major memory layout.
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
            !! The fourth-order minor-symmetric tensor \(\mathbb{A}\) (stored as a 6x6 matrix).
        type(ten_3D2Osym), intent(in) :: b
            !! The second-order symmetric tensor \(\mathbf{b}\) (stored as a 6-component vector).
        type(ten_3D2Osym) :: res
            !! The resulting second-order symmetric tensor \(\mathbf{res}\).
        
        res%vals(:) =   a%vals(:,1)*b%vals(1) + a%vals(:,2)*b%vals(2) + a%vals(:,3)*b%vals(3) +    &
                      2.0D0*a%vals(:,4)*b%vals(4) + 2.0D0*a%vals(:,5)*b%vals(5) + 2.0D0*a%vals(:,6)*b%vals(6)
        
    end function ddot_3D4O2sym_3D2Osym

    pure function ddot_3D2Osym_3D4O2sym(b, a) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbf{b} : \mathbb{A}\).
        !!
        !! Mathematically: \(res_{ij} = b_{kl} A_{klij}\).
        !! In Voigt notation, this is a vector-matrix product: \([\text{res}]_I = [\mathbf{b}]_J [\mathbb{A}]_{JI}\).
        !! This is implemented using the intrinsic `matmul` for optimal performance,
        !! which typically maps to a highly optimized BLAS DGEMV routine.
        !! Shear components of \(\mathbf{b}\) (indices 4, 5, 6) are weighted by 2.
        !!
        implicit none
        type(ten_3D2Osym), intent(in) :: b
            !! The second-order symmetric tensor \(\mathbf{b}\) (stored as a 6-component vector).
        type(ten_3D4O2sym), intent(in) :: a
            !! The fourth-order minor-symmetric tensor \(\mathbb{A}\) (stored as a 6x6 matrix).
        type(ten_3D2Osym) :: res
            !! The resulting second-order symmetric tensor \(\mathbf{res}\).
        
        real(real64), dimension(6) :: b_weighted
        
        ! Create the weighted Voigt vector for contraction.
        b_weighted(1:3) = b%vals(1:3)
        b_weighted(4:6) = 2.0D0 * b%vals(4:6)

        ! Perform vector-matrix multiplication using matmul.
        ! matmul(vector, matrix) performs v * A, which is the correct operation here.
        res%vals = matmul(b_weighted, a%vals)

    end function ddot_3D2Osym_3D4O2sym

    pure function ddot_3D4O3sym_3D2Osym(a, b) result(res)
        !! Computes the double contraction product: res = A : b.
        !! This is a manually unrolled matrix-vector multiplication in Voigt space
        !! to avoid temporary array allocation and maximize performance.
        implicit none
        class(ten_3D4O3sym), intent(in) :: a
            !! The fully symmetric 4th-order tensor A (21 components).
        class(ten_3D2Osym), intent(in) :: b
            !! The symmetric 2nd-order tensor b (6 components).
        type(ten_3D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor.
        
        real(real64) :: b1, b2, b3, b4, b5, b6

        ! Use weighted components for shear terms to account for the factor of 2
        b1 = b%vals(1); b2 = b%vals(2); b3 = b%vals(3)
        b4 = 2.0D0 * b%vals(4)
        b5 = 2.0D0 * b%vals(5)
        b6 = 2.0D0 * b%vals(6)
        
        ! Explicit matrix-vector product using the 21-component storage scheme
        ! res(1) = A(1,J) * b(J)
        res%vals(1) = a%vals(1)*b1 + a%vals(7)*b2 + a%vals(12)*b3 + a%vals(16)*b4 + a%vals(19)*b5 + a%vals(21)*b6
        ! res(2) = A(2,J) * b(J)
        res%vals(2) = a%vals(7)*b1 + a%vals(2)*b2 + a%vals(8)*b3  + a%vals(13)*b4 + a%vals(17)*b5 + a%vals(20)*b6
        ! res(3) = A(3,J) * b(J)
        res%vals(3) = a%vals(12)*b1+ a%vals(8)*b2 + a%vals(3)*b3  + a%vals(9)*b4  + a%vals(14)*b5 + a%vals(18)*b6
        ! res(4) = A(4,J) * b(J)
        res%vals(4) = a%vals(16)*b1+ a%vals(13)*b2+ a%vals(9)*b3  + a%vals(4)*b4  + a%vals(10)*b5 + a%vals(15)*b6
        ! res(5) = A(5,J) * b(J)
        res%vals(5) = a%vals(19)*b1+ a%vals(17)*b2+ a%vals(14)*b3 + a%vals(10)*b4 + a%vals(5)*b5  + a%vals(11)*b6
        ! res(6) = A(6,J) * b(J)
        res%vals(6) = a%vals(21)*b1+ a%vals(20)*b2+ a%vals(18)*b3 + a%vals(15)*b4 + a%vals(11)*b5 + a%vals(6)*b6

    end function ddot_3D4O3sym_3D2Osym

    pure function ddot_3D2Osym_3D4O3sym(b, a) result(res)
        !! Computes the double contraction product: res = b : A.
        !! This is a manually unrolled matrix-vector multiplication in Voigt space
        !! to avoid temporary array allocation and maximize performance.
        implicit none
        class(ten_3D4O3sym), intent(in) :: a
            !! The fully symmetric 4th-order tensor A (21 components).
        class(ten_3D2Osym), intent(in) :: b
            !! The symmetric 2nd-order tensor b (6 components).
        type(ten_3D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor.
        
        real(real64) :: b1, b2, b3, b4, b5, b6

        ! Use weighted components for shear terms to account for the factor of 2
        b1 = b%vals(1); b2 = b%vals(2); b3 = b%vals(3)
        b4 = 2.0D0 * b%vals(4)
        b5 = 2.0D0 * b%vals(5)
        b6 = 2.0D0 * b%vals(6)
        
        ! Explicit matrix-vector product using the 21-component storage scheme
        ! res(1) = A(1,J) * b(J)
        res%vals(1) = a%vals(1)*b1 + a%vals(7)*b2 + a%vals(12)*b3 + a%vals(16)*b4 + a%vals(19)*b5 + a%vals(21)*b6
        ! res(2) = A(2,J) * b(J)
        res%vals(2) = a%vals(7)*b1 + a%vals(2)*b2 + a%vals(8)*b3  + a%vals(13)*b4 + a%vals(17)*b5 + a%vals(20)*b6
        ! res(3) = A(3,J) * b(J)
        res%vals(3) = a%vals(12)*b1+ a%vals(8)*b2 + a%vals(3)*b3  + a%vals(9)*b4  + a%vals(14)*b5 + a%vals(18)*b6
        ! res(4) = A(4,J) * b(J)
        res%vals(4) = a%vals(16)*b1+ a%vals(13)*b2+ a%vals(9)*b3  + a%vals(4)*b4  + a%vals(10)*b5 + a%vals(15)*b6
        ! res(5) = A(5,J) * b(J)
        res%vals(5) = a%vals(19)*b1+ a%vals(17)*b2+ a%vals(14)*b3 + a%vals(10)*b4 + a%vals(5)*b5  + a%vals(11)*b6
        ! res(6) = A(6,J) * b(J)
        res%vals(6) = a%vals(21)*b1+ a%vals(20)*b2+ a%vals(18)*b3 + a%vals(15)*b4 + a%vals(11)*b5 + a%vals(6)*b6

    end function ddot_3D2Osym_3D4O3sym

    ! =========================================================================
    ! 6. 4TH ORDER TENSORS : 4TH ORDER TENSORS
    ! =========================================================================

    pure function ddot_3D4O3sym_3D4O3sym(aT, bT) result(res)
        !! Computes the double contraction of two fully symmetric 4th-order tensors.
        !! \([\mathbf{res}]_{ijmn} = [\mathbf{A}]_{ijkl} [\mathbf{B}]_{klmn}\).
        !!
        !! This implementation is fully unrolled for maximum performance, avoiding temporary arrays
        !! and leveraging compiler optimizations like register allocation and SIMD.
        !!
        implicit none
        class(ten_3D4O3sym), intent(in) :: aT
        class(ten_3D4O3sym), intent(in) :: bT
        type(ten_3D4O2sym) :: res
        
        real(real64) :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19
        real(real64) :: x20, x21, x22, x23, x24, x25, x26, x27, x28, x29

        x0 = aT%vals(7)*bT%vals(7)
        x1 = aT%vals(12)*bT%vals(12)
        x2 = aT%vals(16)*bT%vals(16)
        x3 = aT%vals(21)*bT%vals(21)
        x4 = aT%vals(19)*bT%vals(19)
        x5 = 2*aT%vals(16)
        x6 = 2*aT%vals(21)
        x7 = 2*aT%vals(19)
        x8 = 2*aT%vals(13)
        x9 = 2*aT%vals(20)
        x10 = 2*aT%vals(17)
        x11 = aT%vals(8)*bT%vals(8)
        x12 = aT%vals(13)*bT%vals(13)
        x13 = aT%vals(20)*bT%vals(20)
        x14 = aT%vals(17)*bT%vals(17)
        x15 = 2*aT%vals(9)
        x16 = 2*aT%vals(18)
        x17 = 2*aT%vals(14)
        x18 = aT%vals(9)*bT%vals(9)
        x19 = aT%vals(18)*bT%vals(18)
        x20 = aT%vals(14)*bT%vals(14)
        x21 = 2*aT%vals(4)
        x22 = 2*aT%vals(15)
        x23 = 2*aT%vals(10)
        x24 = bT%vals(15)*x22
        x25 = bT%vals(10)*x23
        x26 = 2*aT%vals(11)
        x27 = 2*aT%vals(5)
        x28 = bT%vals(11)*x26
        x29 = 2*aT%vals(6)
        res%vals(1, 1) = aT%vals(1)*bT%vals(1) + x0 + x1 + 2*x2 + 2*x3 + 2*x4
        res%vals(2, 1) = aT%vals(7)*bT%vals(1) + aT%vals(2)*bT%vals(7) + aT%vals(8)*bT%vals(12) + bT%vals(16)*x8 + &
                         bT%vals(21)*x9 + bT%vals(19)*x10
        res%vals(3, 1) = aT%vals(12)*bT%vals(1) + aT%vals(8)*bT%vals(7) + aT%vals(3)*bT%vals(12) + bT%vals(16)*x15 + &
                         bT%vals(21)*x16 + bT%vals(19)*x17
        res%vals(4, 1) = aT%vals(16)*bT%vals(1) + aT%vals(13)*bT%vals(7) + aT%vals(9)*bT%vals(12) + bT%vals(16)*x21 + &
                         bT%vals(21)*x22 + bT%vals(19)*x23
        res%vals(5, 1) = aT%vals(19)*bT%vals(1) + aT%vals(17)*bT%vals(7) + aT%vals(14)*bT%vals(12) + bT%vals(16)*x23 + &
                         bT%vals(21)*x26 + bT%vals(19)*x27
        res%vals(6, 1) = aT%vals(21)*bT%vals(1) + aT%vals(20)*bT%vals(7) + aT%vals(18)*bT%vals(12) + bT%vals(16)*x22 + &
                         bT%vals(21)*x29 + bT%vals(19)*x26
        res%vals(1, 2) = aT%vals(1)*bT%vals(7) + aT%vals(7)*bT%vals(2) + aT%vals(12)*bT%vals(8) + bT%vals(13)*x5 + &
                         bT%vals(20)*x6 + bT%vals(17)*x7
        res%vals(2, 2) = aT%vals(2)*bT%vals(2) + x0 + x11 + 2*x12 + 2*x13 + 2*x14
        res%vals(3, 2) = aT%vals(12)*bT%vals(7) + aT%vals(8)*bT%vals(2) + aT%vals(3)*bT%vals(8) + bT%vals(13)*x15 + &
                         bT%vals(20)*x16 + bT%vals(17)*x17
        res%vals(4, 2) = aT%vals(16)*bT%vals(7) + aT%vals(13)*bT%vals(2) + aT%vals(9)*bT%vals(8) + bT%vals(13)*x21 + &
                         bT%vals(20)*x22 + bT%vals(17)*x23
        res%vals(5, 2) = aT%vals(19)*bT%vals(7) + aT%vals(17)*bT%vals(2) + aT%vals(14)*bT%vals(8) + bT%vals(13)*x23 + &
                         bT%vals(20)*x26 + bT%vals(17)*x27
        res%vals(6, 2) = aT%vals(21)*bT%vals(7) + aT%vals(20)*bT%vals(2) + aT%vals(18)*bT%vals(8) + bT%vals(13)*x22 + &
                         bT%vals(20)*x29 + bT%vals(17)*x26
        res%vals(1, 3) = aT%vals(1)*bT%vals(12) + aT%vals(7)*bT%vals(8) + aT%vals(12)*bT%vals(3) + bT%vals(9)*x5 + &
                         bT%vals(18)*x6 + bT%vals(14)*x7
        res%vals(2, 3) = aT%vals(7)*bT%vals(12) + aT%vals(2)*bT%vals(8) + aT%vals(8)*bT%vals(3) + bT%vals(9)*x8 + &
                         bT%vals(18)*x9 + bT%vals(14)*x10
        res%vals(3, 3) = aT%vals(3)*bT%vals(3) + x1 + x11 + 2*x18 + 2*x19 + 2*x20
        res%vals(4, 3) = aT%vals(16)*bT%vals(12) + aT%vals(13)*bT%vals(8) + aT%vals(9)*bT%vals(3) + bT%vals(9)*x21 + &
                         bT%vals(18)*x22 + bT%vals(14)*x23
        res%vals(5, 3) = aT%vals(19)*bT%vals(12) + aT%vals(17)*bT%vals(8) + aT%vals(14)*bT%vals(3) + bT%vals(9)*x23 + &
                         bT%vals(18)*x26 + bT%vals(14)*x27
        res%vals(6, 3) = aT%vals(21)*bT%vals(12) + aT%vals(20)*bT%vals(8) + aT%vals(18)*bT%vals(3) + bT%vals(9)*x22 + &
                         bT%vals(18)*x29 + bT%vals(14)*x26
        res%vals(1, 4) = aT%vals(1)*bT%vals(16) + aT%vals(7)*bT%vals(13) + aT%vals(12)*bT%vals(9) + bT%vals(4)*x5 + &
                         bT%vals(15)*x6 + bT%vals(10)*x7
        res%vals(2, 4) = aT%vals(7)*bT%vals(16) + aT%vals(2)*bT%vals(13) + aT%vals(8)*bT%vals(9) + bT%vals(4)*x8 + &
                         bT%vals(15)*x9 + bT%vals(10)*x10
        res%vals(3, 4) = aT%vals(12)*bT%vals(16) + aT%vals(8)*bT%vals(13) + aT%vals(3)*bT%vals(9) + bT%vals(4)*x15 + &
                         bT%vals(15)*x16 + bT%vals(10)*x17
        res%vals(4, 4) = bT%vals(4)*x21 + x12 + x18 + x2 + x24 + x25
        res%vals(5, 4) = aT%vals(19)*bT%vals(16) + aT%vals(17)*bT%vals(13) + aT%vals(14)*bT%vals(9) + bT%vals(4)*x23 + &
                         bT%vals(15)*x26 + bT%vals(10)*x27
        res%vals(6, 4) = aT%vals(21)*bT%vals(16) + aT%vals(20)*bT%vals(13) + aT%vals(18)*bT%vals(9) + bT%vals(4)*x22 + &
                         bT%vals(15)*x29 + bT%vals(10)*x26
        res%vals(1, 5) = aT%vals(1)*bT%vals(19) + aT%vals(7)*bT%vals(17) + aT%vals(12)*bT%vals(14) + bT%vals(10)*x5 + &
                         bT%vals(11)*x6 + bT%vals(5)*x7
        res%vals(2, 5) = aT%vals(7)*bT%vals(19) + aT%vals(2)*bT%vals(17) + aT%vals(8)*bT%vals(14) + bT%vals(10)*x8 + &
                         bT%vals(11)*x9 + bT%vals(5)*x10
        res%vals(3, 5) = aT%vals(12)*bT%vals(19) + aT%vals(8)*bT%vals(17) + aT%vals(3)*bT%vals(14) + bT%vals(10)*x15 + &
                         bT%vals(11)*x16 + bT%vals(5)*x17
        res%vals(4, 5) = aT%vals(16)*bT%vals(19) + aT%vals(13)*bT%vals(17) + aT%vals(9)*bT%vals(14) + bT%vals(10)*x21 + &
                         bT%vals(11)*x22 + bT%vals(5)*x23
        res%vals(5, 5) = bT%vals(5)*x27 + x14 + x20 + x25 + x28 + x4
        res%vals(6, 5) = aT%vals(21)*bT%vals(19) + aT%vals(20)*bT%vals(17) + aT%vals(18)*bT%vals(14) + bT%vals(10)*x22 + &
                         bT%vals(11)*x29 + bT%vals(5)*x26
        res%vals(1, 6) = aT%vals(1)*bT%vals(21) + aT%vals(7)*bT%vals(20) + aT%vals(12)*bT%vals(18) + bT%vals(15)*x5 + &
                         bT%vals(6)*x6 + bT%vals(11)*x7
        res%vals(2, 6) = aT%vals(7)*bT%vals(21) + aT%vals(2)*bT%vals(20) + aT%vals(8)*bT%vals(18) + bT%vals(15)*x8 + &
                         bT%vals(6)*x9 + bT%vals(11)*x10
        res%vals(3, 6) = aT%vals(12)*bT%vals(21) + aT%vals(8)*bT%vals(20) + aT%vals(3)*bT%vals(18) + bT%vals(15)*x15 + &
                         bT%vals(6)*x16 + bT%vals(11)*x17
        res%vals(4, 6) = aT%vals(16)*bT%vals(21) + aT%vals(13)*bT%vals(20) + aT%vals(9)*bT%vals(18) + bT%vals(15)*x21 + &
                         bT%vals(6)*x22 + bT%vals(11)*x23
        res%vals(5, 6) = aT%vals(19)*bT%vals(21) + aT%vals(17)*bT%vals(20) + aT%vals(14)*bT%vals(18) + bT%vals(15)*x23 + &
                         bT%vals(6)*x26 + bT%vals(11)*x27
        res%vals(6, 6) = bT%vals(6)*x29 + x13 + x19 + x24 + x28 + x3

    end function ddot_3D4O3sym_3D4O3sym

end module mod_ops_contraction_double