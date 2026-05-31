module mod_operator_I4O3TS_3D4O3sym
    !! Module mod_operator_I4O3TS_3D4O3sym
    !! ======================================
    !!
    !! Defines mixed algebraic operations between the scaled 3D fourth-order 
    !! identity tensor of type 3 (`iden_4O3TS`, \(c \cdot \delta_{ij}\delta_{kl}\)) 
    !! and 3D fully symmetric fourth-order tensors (`ten_3D4O3sym`).
    !!
    !! In the 21-component compressed Voigt storage, this identity tensor affects 
    !! the normal-normal interaction components: indices 1, 2, 3 (purely normal) 
    !! and 7, 8, 12 (normal coupling).
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of a tensor and the scaled type-3 4th-order identity.
    !! - `-` : Subtraction between a tensor and the scaled type-3 4th-order identity.
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3TS
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: assignment (=)
    interface assignment (=)
        module procedure assign_3D4O3sym_I4O3TS
    end interface
    
contains
    pure subroutine assign_3D4O3sym_I4O3TS(a, b)
        implicit none
        type(ten_3D4O3sym), intent(out) :: a
            !! The target general tensor to be overwritten.
        type(iden_4O3TS), intent(in) :: b  
            !! The source symmetric tensor.
        
        a%vals = 0D0
        a%vals(1:3) = b%val
        a%vals(7:8) = b%val
        a%vals(12)  = b%val
    end subroutine assign_3D4O3sym_I4O3TS

end module mod_operator_I4O3TS_3D4O3sym