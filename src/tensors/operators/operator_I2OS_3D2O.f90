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
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_3D2O_I2OS
        module procedure sum_I2OS_3D2O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D2O_I2OS
        module procedure sub_I2OS_3D2O
    end interface
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_3D2O_I2OS
        module procedure mul_I2OS_3D2O
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2OS_3D2O
        module procedure ddot_3D2O_I2OS
    end interface

contains

    pure function ddot_I2OS_3D2O(I2, b) result(res)
        !! Computes the scaled double contraction \(c\mathbf{I} : \mathbf{B} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
            !! Scaled identity tensor \(c\mathbf{I}\).
        type(ten_3D2O), intent(in) :: b
            !! General 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! Scalar result.
        res = I2%val * (b%vals(1) + b%vals(5) + b%vals(9))
    end function ddot_I2OS_3D2O

    pure function ddot_3D2O_I2OS(b, I2) result(res)
        !! Computes the scaled double contraction \(\mathbf{B} : c\mathbf{I} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: b
        real(real64) :: res
        res = I2%val * (b%vals(1) + b%vals(5) + b%vals(9))
    end function ddot_3D2O_I2OS

    pure function sum_I2OS_3D2O(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = c\mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) + I2%val
        res%vals(5) = res%vals(5) + I2%val
        res%vals(9) = res%vals(9) + I2%val
    end function sum_I2OS_3D2O

    pure function sum_3D2O_I2OS(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) + I2%val
        res%vals(5) = res%vals(5) + I2%val
        res%vals(9) = res%vals(9) + I2%val
    end function sum_3D2O_I2OS

    pure function sub_I2OS_3D2O(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = c\mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = -a%vals
        res%vals(1) = I2%val + res%vals(1)
        res%vals(5) = I2%val + res%vals(5)
        res%vals(9) = I2%val + res%vals(9)
    end function sub_I2OS_3D2O

    pure function sub_3D2O_I2OS(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res
        res%vals = a%vals
        res%vals(1) = res%vals(1) - I2%val
        res%vals(5) = res%vals(5) - I2%val
        res%vals(9) = res%vals(9) - I2%val
    end function sub_3D2O_I2OS

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