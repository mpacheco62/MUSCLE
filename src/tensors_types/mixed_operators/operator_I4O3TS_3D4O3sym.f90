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
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O3sym_I4O3TS
        module procedure sum_I4O3TS_3D4O3sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O3sym_I4O3TS
        module procedure sub_I4O3TS_3D4O3sym
    end interface
    
contains

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

end module mod_operator_I4O3TS_3D4O3sym