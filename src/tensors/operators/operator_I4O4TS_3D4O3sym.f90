module mod_operator_I4O4TS_3D4O3sym
    !! Module mod_operator_I4O4TS_3D4O3sym
    !! ===================================
    !!
    !! Defines mixed algebraic operations between the scaled 3D fourth-order 
    !! symmetric identity tensor (`iden_4O4TS`, \(c\mathbb{I}^S\)) and 
    !! 3D fully symmetric fourth-order tensors (`ten_3D4O3sym`).
    !!
    !! The `ten_3D4O3sym` type uses a compressed 21-component storage.
    !! In this format, the scaled identity \(c\mathbb{I}^S\) affects only the 
    !! first six components (the diagonal of the 6x6 Voigt matrix):
    !! - Components 1, 2, 3 (Normal): + c
    !! - Components 4, 5, 6 (Shear):  + 0.5 * c
    !!
    !! This module provides overloads for addition and subtraction.

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4TS
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: assignment (=)
    interface assignment (=)
        module procedure assign_3D4O3sym_I4O4TS
    end interface

contains

    pure subroutine assign_3D4O3sym_I4O4TS(a, b)
        implicit none
        type(ten_3D4O3sym), intent(out) :: a
            !! The target general tensor to be overwritten.
        type(iden_4O4TS), intent(in) :: b  
            !! The source symmetric tensor.
        real(real64) :: c, half_c

        c = b%val
        half_c = 0.5D0 * c

        a%vals = 0D0
        ! Add scaled identity components to the Voigt diagonal
        a%vals(1:3) = c
        a%vals(4:6) = half_c
    end subroutine assign_3D4O3sym_I4O4TS
    
end module mod_operator_I4O4TS_3D4O3sym