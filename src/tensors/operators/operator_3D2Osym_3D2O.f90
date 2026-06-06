module mod_operator_3D2Osym_3D2O
    !! Module mod_operator_3D2Osym_3D2O
    !! ================================
    !!
    !! Defines mixed algebraic operations involving 3D symmetric second-order 
    !! tensors (`ten_3D2Osym`) and general (non-symmetric) 3D second-order 
    !! tensors (`ten_3D2O`).
    !!
    !! Since operations between a symmetric and a non-symmetric tensor generally 
    !! break symmetry, the resulting tensors from these operations are mostly 
    !! cast into the general `ten_3D2O` format (9 components, column-major).
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition.
    !!   - `3D2O + 3D2Osym` -> `3D2O`
    !!   - `3D2Osym + 3D2O` -> `3D2O`
    !!
    !! - `-` : Subtraction.
    !!   - `3D2O - 3D2Osym` -> `3D2O`
    !!   - `3D2Osym - 3D2O` -> `3D2O`
    !!
    !! - `*` : Single tensor contraction (Matrix dot product).
    !!   - \(\mathbf{C} = \mathbf{A} \cdot \mathbf{B}\)
    !!   - `3D2O * 3D2Osym` -> `3D2O`
    !!   - `3D2Osym * 3D2O` -> `3D2O`
    !!   - `3D2Osym * 3D2Osym` -> `3D2O` (Since the product of two symmetric tensors is generally non-symmetric).
    !!
    !! - `.ddot.` : Double tensor contraction (Frobenius inner product).
    !!   - \(\alpha = \mathbf{A} : \mathbf{B}\) -> `real(real64)`
    !!
    !! - `=` : Assignment.
    !!   - Implicit casting from a symmetric tensor to a general tensor format.

    use, intrinsic :: iso_fortran_env
    use mod_ten_3D2Osym
    use mod_ten_3D2O
    implicit none
    private

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_3D2O_3D2Osym
    end interface

contains

    ! =========================================================================
    ! ASSIGNMENT
    ! =========================================================================

    pure subroutine assign_3D2O_3D2Osym(a, b)
        !! Explicit assignment from a symmetric tensor to a general tensor.
        !! 
        !! Unpacks the 6-component array of the symmetric tensor into the 
        !! 9-component column-major array of the general tensor.
        implicit none
        type(ten_3D2O), intent(out) :: a
            !! The target general tensor to be overwritten.
        type(ten_3D2Osym), intent(in) :: b  
            !! The source symmetric tensor.
        
        a%vals(1) = b%vals(1) ! xx
        a%vals(2) = b%vals(4) ! yx <- xy
        a%vals(3) = b%vals(6) ! zx <- xz
        a%vals(4) = b%vals(4) ! xy
        a%vals(5) = b%vals(2) ! yy
        a%vals(6) = b%vals(5) ! zy <- yz
        a%vals(7) = b%vals(6) ! xz
        a%vals(8) = b%vals(5) ! yz
        a%vals(9) = b%vals(3) ! zz
    end subroutine assign_3D2O_3D2Osym

end module mod_operator_3D2Osym_3D2O