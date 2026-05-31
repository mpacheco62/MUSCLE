module mod_operator_I2OS_3D2O
    !! Module mod_operator_I2OS_3D2O
    !! =====================================
    !!
    !! Defines mixed algebraic operations between the scaled 2nd-order identity 
    !! tensor (`iden_2OS`, \(c\mathbf{I}\)) and general 3D second-order 
    !! tensors (`ten_3D2O`).
    !!
    !! This module handles the 9-component column-major storage format where
    !! diagonal elements are located at indices 1, 5, and 9.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of a tensor and a scaled identity.
    !! - `-` : Subtraction between a tensor and a scaled identity.
    !! - `*` : Scalar-tensor multiplication (contraction with scaled identity).
    !! - `.ddot.` : Double contraction, computing \(c \cdot \text{tr}(\mathbf{B})\).
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_2OS
    use mod_ten_3D2O
    implicit none
    private
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_3D2O_I2OS
        module procedure mul_I2OS_3D2O
    end interface

contains

    pure function mul_I2OS_3D2O(I2, a) result(res)
        !! Computes the scaled contraction \(\mathbf{res} = c\mathbf{I} \cdot \mathbf{A} = c\mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = I2%val * a%vals
    end function mul_I2OS_3D2O

    pure function mul_3D2O_I2OS(a, I2) result(res)
        !! Computes the scaled contraction \(\mathbf{res} = \mathbf{A} \cdot c\mathbf{I} = c\mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = I2%val * a%vals
    end function mul_3D2O_I2OS

end module mod_operator_I2OS_3D2O