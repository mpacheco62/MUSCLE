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

    public :: operator(.tdotsym.)
    interface operator (.tdotsym.)
        module procedure tdotsym_3D2Osym_3D2Osym
        module procedure tdotsym_3D2Osym
    end interface

contains

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