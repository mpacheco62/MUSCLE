module mod_operator_I4O3TS_3D4O2sym
    !! Module mod_operator_I4O3TS_3D4O2sym
    !! ======================================
    !!
    !! Defines mixed algebraic operations between the scaled 3D fourth-order 
    !! identity tensor of type 3 (`iden_4O3TS`, \(c \cdot \delta_{ij}\delta_{kl}\)) 
    !! and 3D fourth-order tensors with minor symmetries (`ten_3D4O2sym`).
    !!
    !! In the 6x6 Voigt matrix representation, this identity tensor affects 
    !! only the upper-left 3x3 block, which corresponds to normal-normal 
    !! stress/strain components.
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
    use mod_ten_3D4O2sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O2sym_I4O3TS
        module procedure sum_I4O3TS_3D4O2sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O2sym_I4O3TS
        module procedure sub_I4O3TS_3D4O2sym
    end interface
    
contains

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

end module mod_operator_I4O3TS_3D4O2sym