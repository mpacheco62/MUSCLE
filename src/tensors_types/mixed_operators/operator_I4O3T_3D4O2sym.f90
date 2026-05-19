! --- INICIO DEL ARCHIVO CORREGIDO: operator_I4O3T_3D4O2sym.f90 ---

module mod_operator_I4O3T_3D4O2sym
    !! author: MPacheco
    !! version: 1.1 - Added FORD documentation and optimized argument types.
    !!
    !! Module mod_operator_I4O3T_3D4O2sym
    !! =====================================
    !!
    !! Defines mixed algebraic operations between the standard 3D fourth-order 
    !! identity tensor of type 3 (`iden_4O3T`, \(\delta_{ij}\delta_{kl}\)) and 
    !! 3D fourth-order tensors with minor symmetries (`ten_3D4O2sym`).
    !!
    !! In the 6x6 Voigt matrix representation, this identity tensor corresponds 
    !! to a matrix with 1.0 in the upper-left 3x3 block and 0.0 elsewhere.
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
    use mod_ten_3D4O2sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O2sym_I4O3T
        module procedure sum_I4O3T_3D4O2sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O2sym_I4O3T
        module procedure sub_I4O3T_3D4O2sym
    end interface
    
contains

    pure function sum_I4O3T_3D4O2sym(I4, a) result(res)
        !! Computes the sum \(\mathbb{res} = \mathbb{I} + \mathbb{A}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
            !! Standard 4th-order identity tensor (\(\delta_{ij}\delta_{kl}\)).
        type(ten_3D4O2sym), intent(in) :: a
            !! 4th-order tensor with minor symmetries.
        type(ten_3D4O2sym) :: res
            !! Resulting 4th-order tensor.
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1.0D0
    end function sum_I4O3T_3D4O2sym

    pure function sum_3D4O2sym_I4O3T(a, I4) result(res)
        !! Computes the sum \(\mathbb{res} = \mathbb{A} + \mathbb{I}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1.0D0
    end function sum_3D4O2sym_I4O3T

    pure function sub_I4O3T_3D4O2sym(I4, a) result(res)
        !! Computes the subtraction \(\mathbb{res} = \mathbb{I} - \mathbb{A}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = -a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1.0D0
    end function sub_I4O3T_3D4O2sym

    pure function sub_3D4O2sym_I4O3T(a, I4) result(res)
        !! Computes the subtraction \(\mathbb{res} = \mathbb{A} - \mathbb{I}\).
        implicit none
        type(iden_4O3T), intent(in) :: I4
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = a%vals
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) - 1.0D0
    end function sub_3D4O2sym_I4O3T

end module mod_operator_I4O3T_3D4O2sym