module mod_operator_3D2Osym_3D4O2sym
    !! Module mod_operator_3D2Osym_3D4O2sym
    !! ====================================
    !!
    !! Defines mixed algebraic operations involving 3D symmetric second-order 
    !! tensors (`ten_3D2Osym`) and 3D fourth-order tensors with only minor 
    !! symmetries (`ten_3D4O2sym`).
    !!
    !! Since `ten_3D4O2sym` lacks major symmetry (its 6x6 Voigt matrix is not
    !! necessarily symmetric), the order of operands in contractions matters.
    !! This module provides correct and optimized implementations for these operations.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `.ddot.` : Double tensor contraction.
    !!   - \(\mathbf{C} = \mathbb{A} : \mathbf{B}\) (`3D4O2sym .ddot. 3D2Osym` -> `3D2Osym`)
    !!   - \(\mathbf{C} = \mathbf{B} : \mathbb{A}\) (`3D2Osym .ddot. 3D4O2sym` -> `3D2Osym`)
    !!
    !! - `.tdot.` : Dyadic tensor product (Outer product).
    !!   - \(\mathbb{C} = \mathbf{a} \otimes \mathbf{b}\) (`3D2Osym .tdot. 3D2Osym` -> `3D4O2sym`)
    !!     *The result correctly has only minor symmetries.*
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_ten_3D2Osym
    use mod_ten_3D4O2sym
    implicit none
    private

    public :: operator(.tdot.)
    interface operator (.tdot.)
        module procedure tdot_3D2Osym_3D2Osym
    end interface

contains
    ! =========================================================================
    ! DYADIC PRODUCT (.tdot.)
    ! =========================================================================
    
    pure function tdot_3D2Osym_3D2Osym(a, b) result(res)
        !! Computes the dyadic (outer) tensor product: \(\mathbb{C} = \mathbf{a} \otimes \mathbf{b}\).
        !!
        !! Mathematically: \( C_{ijkl} = a_{ij} b_{kl} \).
        !!
        !! The result is a fourth-order tensor that only has minor symmetries, even if
        !! \(\mathbf{a}\) and \(\mathbf{b}\) are symmetric. Major symmetry (\(C_{ijkl} = C_{klij}\))
        !! is only guaranteed if \(\mathbf{a}\) and \(\mathbf{b}\) are proportional.
        !! The result is therefore correctly returned as a `ten_3D4O2sym` type.
        !!
        !! In Voigt notation, this operation constructs a 6x6 matrix by taking the
        !! outer product of the 6-component Voigt vectors of \(\mathbf{a}\) and \(\mathbf{b}\).
        !!
        implicit none
        type(ten_3D2Osym), intent(in) :: a  
            !! The first second-order symmetric tensor \(\mathbf{a}\).
        type(ten_3D2Osym), intent(in) :: b  
            !! The second second-order symmetric tensor \(\mathbf{b}\).
        type(ten_3D4O2sym) :: res         
            !! The resulting fourth-order tensor with minor symmetries \(\mathbb{C}\).

        ! Build the 6x6 Voigt matrix C_IJ = a_I * b_J by constructing each row.
        res%vals(1,:) = a%vals(1) * b%vals(:)
        res%vals(2,:) = a%vals(2) * b%vals(:)
        res%vals(3,:) = a%vals(3) * b%vals(:)
        res%vals(4,:) = a%vals(4) * b%vals(:)
        res%vals(5,:) = a%vals(5) * b%vals(:)
        res%vals(6,:) = a%vals(6) * b%vals(:)

    end function tdot_3D2Osym_3D2Osym

end module mod_operator_3D2Osym_3D4O2sym