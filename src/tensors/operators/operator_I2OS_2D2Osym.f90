! --- INICIO DEL ARCHIVO CORREGIDO: operator_I2OS_2D2Osym.f90 ---

module mod_operator_I2OS_2D2Osym
    !! author: MPacheco
    !! version: 1.1 - Added FORD documentation and micro-optimizations.
    !!
    !! Module mod_operator_I2OS_2D2Osym
    !! =====================================
    !!
    !! Defines mixed algebraic operations between the scaled 2nd-order identity 
    !! tensor (`iden_2OS`, \(c\mathbf{I}\)) and 2D symmetric second-order 
    !! tensors (`ten_2D2Osym`).
    !!
    !! This module enables operations like `(2.0*I) + stress` or `(1.0*I) : strain`.
    !! It specifically targets the 4-component symmetric storage (xx, yy, zz, xy).
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of a tensor and a scaled identity.
    !! - `-` : Subtraction between a tensor and a scaled identity.
    !! - `.ddot.` : Double contraction, computing \(c \cdot \text{tr}(\mathbf{B})\).
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_2OS
    use mod_ten_2D2Osym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_2D2Osym_I2OS
        module procedure sum_I2OS_2D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D2Osym_I2OS
        module procedure sub_I2OS_2D2Osym
    end interface
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2OS_2D2Osym
        module procedure ddot_2D2Osym_I2OS
    end interface

contains

    pure function ddot_I2OS_2D2Osym(I2, b) result(res)
        !! Computes the scaled double contraction \(c\mathbf{I} : \mathbf{B} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
            !! Scaled identity tensor \(c\mathbf{I}\).
        type(ten_2D2Osym), intent(in) :: b
            !! Symmetric 2nd-order tensor \(\mathbf{B}\).
        real(real64) :: res
            !! Scalar result.
        res = I2%val * (b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_I2OS_2D2Osym

    pure function ddot_2D2Osym_I2OS(b, I2) result(res)
        !! Computes the scaled double contraction \(\mathbf{B} : c\mathbf{I} = c \cdot \text{tr}(\mathbf{B})\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: b
        real(real64) :: res
        res = I2%val * (b%vals(1) + b%vals(2) + b%vals(3))
    end function ddot_2D2Osym_I2OS

    pure function sum_I2OS_2D2Osym(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = c\mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_I2OS_2D2Osym

    pure function sum_2D2Osym_I2OS(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + I2%val
    end function sum_2D2Osym_I2OS

    pure function sub_I2OS_2D2Osym(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = c\mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = -a%vals
        res%vals(1:3) = I2%val + res%vals(1:3)
    end function sub_I2OS_2D2Osym

    pure function sub_2D2Osym_I2OS(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - c\mathbf{I}\).
        implicit none
        type(iden_2OS), intent(in) :: I2
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - I2%val
    end function sub_2D2Osym_I2OS

end module mod_operator_I2OS_2D2Osym