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
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_3D2Osym_I2O
        module procedure mul_I2O_3D2Osym
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_3D2Osym_I2O
    end interface

contains

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
    
    pure subroutine assign_3D2Osym_I2O(a, b)
        !! Explicit assignment from a second order identity to a 3D symmetric tensor.
        !! 
        !! This is necessary to allow statements like `a = I2` where `a` is a `ten_3D2Osym` and `I2` is an `iden_2O`.
        !! The resulting tensor `a` will have its diagonal components set to 1 and off-diagonal components set to 0.
        implicit none
        type(ten_3D2Osym), intent(out) :: a
            !! The target symmetric tensor to be overwritten.
        type(iden_2O), intent(in) :: b
            !! The source second-order identity tensor.

        a%vals(1) = 1.0D0 ! xx
        a%vals(2) = 1.0D0 ! yy
        a%vals(3) = 1.0D0 ! zz
        a%vals(4) = 0.0D0 ! xy
        a%vals(5) = 0.0D0 ! yz
        a%vals(6) = 0.0D0 ! xz
    end subroutine assign_3D2Osym_I2O
end module mod_operator_I2O_3D2Osym