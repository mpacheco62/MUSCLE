module mod_operator_I2OS_3D2Osym
    !! Module mod_operator_I2OS_3D2Osym
    !! ======================================
    !!
    !! Defines mixed algebraic operations between the scaled 2nd-order identity 
    !! tensor (`iden_2OS`, \(c\mathbf{I}\)) and 3D symmetric second-order 
    !! tensors (`ten_3D2Osym`).
    !!
    !! This module handles the 6-component Voigt storage format where diagonal 
    !! elements occupy the first three positions (1, 2, and 3).
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
    use mod_ten_3D2Osym
    implicit none
    private

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_3D2Osym_I2OS
    end interface

contains

    pure subroutine assign_3D2Osym_I2OS(a, b)
        !! Explicit assignment from a scaled second order identity to a 3D symmetric tensor.
        !! 
        !! This is necessary to allow statements like `a = I2` where `a` is a `ten_3D2Osym` and `I2` is an `iden_2OS`.
        !! The resulting tensor `a` will have its diagonal components set to the scale factor and off-diagonal components set to 0.
        implicit none
        type(ten_3D2Osym), intent(out) :: a
            !! The target symmetric tensor to be overwritten.
        type(iden_2OS), intent(in) :: b
            !! The source scaled second-order identity tensor.

        a%vals(1) = b%val ! xx
        a%vals(2) = b%val ! yy
        a%vals(3) = b%val ! zz
        a%vals(4) = 0.0D0 ! xy
        a%vals(5) = 0.0D0 ! yz
        a%vals(6) = 0.0D0 ! xz
    end subroutine assign_3D2Osym_I2OS

end module mod_operator_I2OS_3D2Osym