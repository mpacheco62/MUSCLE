module mod_operator_I4O3T_3D4O3sym
    !! Module mod_operator_I4O3T_3D4O3sym
    !! =====================================
    !!
    !! Defines mixed algebraic operations between the standard 3D fourth-order 
    !! identity tensor of type 3 (`iden_4O3T`, \(\delta_{ij}\delta_{kl}\)) and 
    !! 3D fully symmetric fourth-order tensors (`ten_3D4O3sym`).
    !!
    !! In the 21-component compressed Voigt storage, this identity tensor 
    !! increments the diagonal terms of the normal-normal interaction block:
    !! indices 1, 2, 3 (purely normal) and 7, 8, 12 (normal coupling).
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of a tensor and the type-3 fourth-order identity.
    !! - `-` : Subtraction between a tensor and the type-3 fourth-order identity.
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3T
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O3sym_I4O3T
        module procedure sum_I4O3T_3D4O3sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O3sym_I4O3T
        module procedure sub_I4O3T_3D4O3sym
    end interface
    
contains

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

end module mod_operator_I4O3T_3D4O3sym