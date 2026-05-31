module mod_operator_I2O_2D2Osym
    !! Module mod_operator_I2O_2D2Osym
    !! ====================================
    !!
    !! Defines mixed algebraic operations between the standard 2D identity tensor
    !! (`iden_2O`) and 2D symmetric second-order tensors (`ten_2D2Osym`).
    !!
    !! This module provides overloads for basic arithmetic, allowing for intuitive,
    !! readable tensor expressions like `stress - I` or `I : strain`.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of a tensor and the identity.
    !! - `-` : Subtraction of a tensor and the identity.
    !! - `*` : Single contraction (matrix product) with the identity.
    !! - `.ddot.` : Double contraction, which computes the trace of the tensor.
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_2O
    use mod_ten_2D2Osym
    implicit none
    private
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_2D2Osym_I2O
        module procedure mul_I2O_2D2Osym
    end interface

contains

    pure function mul_I2O_2D2Osym(I2, a) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{I} \cdot \mathbf{A} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
    end function mul_I2O_2D2Osym

    pure function mul_2D2Osym_I2O(a, I2) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{A} \cdot \mathbf{I} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
    end function mul_2D2Osym_I2O
    
end module mod_operator_I2O_2D2Osym