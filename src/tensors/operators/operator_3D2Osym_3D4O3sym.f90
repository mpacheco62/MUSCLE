module mod_operator_3D2Osym_3D4O3sym
    !! Module mod_operator_3D2Osym_3D4O3sym
    !! ====================================
    !!
    !! Defines mixed algebraic operations between 3D symmetric second-order 
    !! tensors (`ten_3D2Osym`) and 3D fully symmetric fourth-order tensors 
    !! (`ten_3D4O3sym`).
    !!
    !! This module provides optimized and mathematically robust implementations for:
    !! - Double contraction (`.ddot.`).
    !! - Symmetrized dyadic product (`.tdotsym.`), overloaded for binary and unary-like use.
    !!
    !! Performance Note
    !! ----------------
    !! The double contraction `ddot` is manually unrolled to eliminate temporary
    !! array creation, maximizing performance in tight loops by reducing memory
    !! overhead and giving the compiler full optimization visibility.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `.ddot.` : Double tensor contraction.
    !!   - \(\mathbf{C} = \mathbb{A} : \mathbf{B}\) (`3D4O3sym .ddot. 3D2Osym` -> `3D2Osym`)
    !!
    !! - `.tdotsym.` : Symmetrized dyadic product (overloaded operator).
    !!   - **Binary Use**: \(\mathbb{C} = \frac{1}{2} (\mathbf{a} \otimes \mathbf{b} + \mathbf{b} \otimes \mathbf{a})\) (called as `a .tdotsym. b`)
    !!   - **Unary-like Use**: \(\mathbb{C} = \mathbf{a} \otimes \mathbf{a}\) (called as `.tdotsym. a`)
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_ten_3D2Osym
    use mod_ten_3D4O3sym
    implicit none
    private

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D4O3sym_3D2Osym
        module procedure ddot_3D2Osym_3D4O3sym
    end interface

    public :: operator(.tdotsym.)
    interface operator (.tdotsym.)
        module procedure tdotsym_3D2Osym_3D2Osym
        module procedure tdotsym_3D2Osym
    end interface

contains

    ! =========================================================================
    ! DOUBLE CONTRACTION
    ! =========================================================================

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
    ! SYMMETRIZED DYADIC PRODUCT (.tdotsym.)
    ! =========================================================================

    pure function tdotsym_3D2Osym_3D2Osym(a, b) result(res)
        !! Computes the symmetrized dyadic product: \(\mathbb{C} = \frac{1}{2} (\mathbf{a} \otimes \mathbf{b} + \mathbf{b} \otimes \mathbf{a})\).
        !! The result is always a fully symmetric 4th-order tensor (`ten_3D4O3sym`).
        !! When called as `a .tdotsym. a`, a modern compiler will optimize this to `a ⊗ a`.
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        class(ten_3D2Osym), intent(in) :: b
        type(ten_3D4O3sym) :: res
        
        res%vals(1)  = 0.5D0 * (a%vals(1)*b%vals(1) + b%vals(1)*a%vals(1))
        res%vals(2)  = 0.5D0 * (a%vals(2)*b%vals(2) + b%vals(2)*a%vals(2))
        res%vals(3)  = 0.5D0 * (a%vals(3)*b%vals(3) + b%vals(3)*a%vals(3))
        res%vals(4)  = 0.5D0 * (a%vals(4)*b%vals(4) + b%vals(4)*a%vals(4))
        res%vals(5)  = 0.5D0 * (a%vals(5)*b%vals(5) + b%vals(5)*a%vals(5))
        res%vals(6)  = 0.5D0 * (a%vals(6)*b%vals(6) + b%vals(6)*a%vals(6))
        res%vals(7)  = 0.5D0 * (a%vals(1)*b%vals(2) + b%vals(1)*a%vals(2))
        res%vals(8)  = 0.5D0 * (a%vals(2)*b%vals(3) + b%vals(2)*a%vals(3))
        res%vals(9)  = 0.5D0 * (a%vals(3)*b%vals(4) + b%vals(3)*a%vals(4))
        res%vals(10) = 0.5D0 * (a%vals(4)*b%vals(5) + b%vals(4)*a%vals(5))
        res%vals(11) = 0.5D0 * (a%vals(5)*b%vals(6) + b%vals(5)*a%vals(6))
        res%vals(12) = 0.5D0 * (a%vals(1)*b%vals(3) + b%vals(1)*a%vals(3))
        res%vals(13) = 0.5D0 * (a%vals(2)*b%vals(4) + b%vals(2)*a%vals(4))
        res%vals(14) = 0.5D0 * (a%vals(3)*b%vals(5) + b%vals(3)*a%vals(5))
        res%vals(15) = 0.5D0 * (a%vals(4)*b%vals(6) + b%vals(4)*a%vals(6))
        res%vals(16) = 0.5D0 * (a%vals(1)*b%vals(4) + b%vals(1)*a%vals(4))
        res%vals(17) = 0.5D0 * (a%vals(2)*b%vals(5) + b%vals(2)*a%vals(5))
        res%vals(18) = 0.5D0 * (a%vals(3)*b%vals(6) + b%vals(3)*a%vals(6))
        res%vals(19) = 0.5D0 * (a%vals(1)*b%vals(5) + b%vals(1)*a%vals(5))
        res%vals(20) = 0.5D0 * (a%vals(2)*b%vals(6) + b%vals(2)*a%vals(6))
        res%vals(21) = 0.5D0 * (a%vals(1)*b%vals(6) + b%vals(1)*a%vals(6))

    end function tdotsym_3D2Osym_3D2Osym

    pure function tdotsym_3D2Osym(a) result(res)
        !! Computes the direct dyadic product of a tensor with itself: \(\mathbb{C} = \mathbf{a} \otimes \mathbf{a}\).
        !! This function provides a syntactically "unary" way to perform the dyadic product.
        !! Mathematically: \( C_{ijkl} = a_{ij} a_{kl} \).
        implicit none
        class(ten_3D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor to be multiplied by itself.
        type(ten_3D4O3sym) :: res
            !! The resulting fully symmetric 4th-order tensor.
        
        ! Compute C_IJ = a_I * a_J for the 21 unique components
        res%vals(1)  = a%vals(1)*a%vals(1)
        res%vals(2)  = a%vals(2)*a%vals(2)
        res%vals(3)  = a%vals(3)*a%vals(3)
        res%vals(4)  = a%vals(4)*a%vals(4)
        res%vals(5)  = a%vals(5)*a%vals(5)
        res%vals(6)  = a%vals(6)*a%vals(6)
        res%vals(7)  = a%vals(1)*a%vals(2)
        res%vals(8)  = a%vals(2)*a%vals(3)
        res%vals(9)  = a%vals(3)*a%vals(4)
        res%vals(10) = a%vals(4)*a%vals(5)
        res%vals(11) = a%vals(5)*a%vals(6)
        res%vals(12) = a%vals(1)*a%vals(3)
        res%vals(13) = a%vals(2)*a%vals(4)
        res%vals(14) = a%vals(3)*a%vals(5)
        res%vals(15) = a%vals(4)*a%vals(6)
        res%vals(16) = a%vals(1)*a%vals(4)
        res%vals(17) = a%vals(2)*a%vals(5)
        res%vals(18) = a%vals(3)*a%vals(6)
        res%vals(19) = a%vals(1)*a%vals(5)
        res%vals(20) = a%vals(2)*a%vals(6)
        res%vals(21) = a%vals(1)*a%vals(6)

    end function tdotsym_3D2Osym

end module mod_operator_3D2Osym_3D4O3sym