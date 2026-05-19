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
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_2D2Osym_I2O
        module procedure sum_I2O_2D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D2Osym_I2O
        module procedure sub_I2O_2D2Osym
    end interface
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_2D2Osym_I2O
        module procedure mul_I2O_2D2Osym
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2O_2D2Osym
        module procedure ddot_2D2Osym_I2O
    end interface

contains

    pure function ddot_I2O_2D2Osym(I2, b) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{I} : \mathbf{B}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: b
            !! The symmetric 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! The resulting scalar (trace of \(\mathbf{B}\)).
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_I2O_2D2Osym

    pure function ddot_2D2Osym_I2O(b, I2) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{B} : \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: b
            !! The symmetric 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! The resulting scalar (trace of \(\mathbf{B}\)).
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_2D2Osym_I2O

    pure function sum_I2O_2D2Osym(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
    end function sum_I2O_2D2Osym

    pure function sum_2D2Osym_I2O(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
    end function sum_2D2Osym_I2O

    pure function sub_I2O_2D2Osym(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = -a%vals
        res%vals(1:3) = 1.0D0 + res%vals(1:3)
    end function sub_I2O_2D2Osym

    pure function sub_2D2Osym_I2O(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
            !! The standard 2nd-order identity tensor \(\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor \(\mathbf{A}\).
        type(ten_2D2Osym) :: res
            !! The resulting symmetric 2nd-order tensor \(\mathbf{res}\).
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - 1.0D0
    end function sub_2D2Osym_I2O

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