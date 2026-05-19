module mod_operator_I2O_3D2Osym
    !! Module mod_operator_I2O_3D2Osym
    !! ====================================
    !!
    !! Defines mixed algebraic operations between the standard 3D identity tensor
    !! (`iden_2O`) and 3D symmetric second-order tensors (`ten_3D2Osym`).
    !!
    !! This module handles the 6-component Voigt storage format where diagonal 
    !! elements occupy the first three positions. It also includes the dyadic 
    !! product with the identity, resulting in 4th-order tensors.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of a tensor and the identity.
    !! - `-` : Subtraction of a tensor and the identity.
    !! - `*` : Single contraction (matrix product) with the identity.
    !! - `.ddot.` : Double contraction, computing the trace of the tensor.
    !! - `.tdot.` : Dyadic product (\(\mathbf{a} \otimes \mathbf{I}\) or \(\mathbf{I} \otimes \mathbf{a}\)).
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_2O
    use mod_ten_3D2Osym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_3D2Osym_I2O
        module procedure sum_I2O_3D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D2Osym_I2O
        module procedure sub_I2O_3D2Osym
    end interface
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_3D2Osym_I2O
        module procedure mul_I2O_3D2Osym
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2O_3D2Osym
        module procedure ddot_3D2Osym_I2O
    end interface

    public :: operator(.tdot.)
    interface operator (.tdot.)
        module procedure tdot_3D2Osym_I2O
        module procedure tdot_I2O_3D2Osym
    end interface

contains

    pure function ddot_I2O_3D2Osym(I2, b) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{I} : \mathbf{B}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_I2O_3D2Osym

    pure function ddot_3D2Osym_I2O(b, I2) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{B} : \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_3D2Osym_I2O

    pure function sum_I2O_3D2Osym(I2, a) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{I} + \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
    end function sum_I2O_3D2Osym

    pure function sum_3D2Osym_I2O(a, I2) result(res)
        !! Computes the sum \(\mathbf{res} = \mathbf{A} + \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) + 1.0D0
    end function sum_3D2Osym_I2O

    pure function sub_I2O_3D2Osym(I2, a) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{I} - \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = -a%vals
        res%vals(1:3) = 1.0D0 + res%vals(1:3)
    end function sub_I2O_3D2Osym

    pure function sub_3D2Osym_I2O(a, I2) result(res)
        !! Computes the subtraction \(\mathbf{res} = \mathbf{A} - \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
        res%vals(1:3) = res%vals(1:3) - 1.0D0
    end function sub_3D2Osym_I2O

    pure function mul_I2O_3D2Osym(I2, a) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{I} \cdot \mathbf{A} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
    end function mul_I2O_3D2Osym

    pure function mul_3D2Osym_I2O(a, I2) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{A} \cdot \mathbf{I} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
    end function mul_3D2Osym_I2O

    pure function tdot_3D2Osym_I2O(a, I2) result(res)
        !! Computes the dyadic product \(\mathbb{C} = \mathbf{a} \otimes \mathbf{I}\).
        !! In 6x6 Voigt notation, this fills the first three columns with the vector \(\mathbf{a}\).
        use mod_ten_3D4O2sym
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = 0.0D0
        res%vals(:,1) = a%vals
        res%vals(:,2) = a%vals
        res%vals(:,3) = a%vals
    end function tdot_3D2Osym_I2O

    pure function tdot_I2O_3D2Osym(I2, a) result(res)
        !! Computes the dyadic product \(\mathbb{C} = \mathbf{I} \otimes \mathbf{a}\).
        !! In 6x6 Voigt notation, this fills the first three rows with the vector \(\mathbf{a}\).
        use mod_ten_3D4O2sym
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = 0.0D0
        res%vals(1,:) = a%vals
        res%vals(2,:) = a%vals
        res%vals(3,:) = a%vals
    end function tdot_I2O_3D2Osym
    
end module mod_operator_I2O_3D2Osym