module mod_operator_I2O_3D2O
    !! Module mod_operator_I2O_3D2O
    !! ============================
    !!
    !! Defines mixed algebraic operations between the standard 3D identity tensor
    !! (`iden_2O`) and general (non-symmetric) 3D second-order tensors (`ten_3D2O`).
    !!
    !! This module handles the 9-component column-major storage format where
    !! diagonal elements are located at indices 1, 5, and 9.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of a tensor and the identity.
    !! - `-` : Subtraction of a tensor and the identity.
    !! - `*` : Single contraction (matrix product) with the identity.
    !! - `.ddot.` : Double contraction (Trace of the tensor).
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_2O
    use mod_ten_3D2O
    implicit none
    private
        
    public :: operator(*)
    interface operator (*)
        module procedure mul_3D2O_I2O
        module procedure mul_I2O_3D2O
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2O_3D2O
        module procedure ddot_3D2O_I2O
    end interface

contains

    pure function ddot_I2O_3D2O(I2, b) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{I} : \mathbf{B}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(5) + b%vals(9)
    end function ddot_I2O_3D2O

    pure function ddot_3D2O_I2O(b, I2) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{B} : \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(5) + b%vals(9)
    end function ddot_3D2O_I2O

    pure function mul_I2O_3D2O(I2, a) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{I} \cdot \mathbf{A} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
    end function mul_I2O_3D2O

    pure function mul_3D2O_I2O(a, I2) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{A} \cdot \mathbf{I} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
    end function mul_3D2O_I2O

end module mod_operator_I2O_3D2O