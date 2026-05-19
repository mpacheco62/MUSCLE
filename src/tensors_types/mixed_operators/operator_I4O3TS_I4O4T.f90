module mod_operator_I4O3TS_I4O4T
    !! Module mod_operator_I4O3TS_I4O4T
    !! ======================================
    !!
    !! Defines mixed algebraic operations between two different types of 3D 
    !! fourth-order identity tensors:
    !! 1. Scaled Type 3 (`iden_4O3TS`): \(\mathbb{I}_{3TS} = c \cdot \delta_{ij}\delta_{kl}\)
    !! 2. Standard Type 4 (`iden_4O4T`): \(\mathbb{I}_{4T} = \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\)
    !!
    !! The result is a fully symmetric 4th-order tensor (`ten_3D4O3sym`) stored 
    !! in the 21-component compressed Voigt format.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of Scaled Type-3 and Standard Type-4 identity tensors.
    !! - `-` : Subtraction between Scaled Type-3 and Standard Type-4 identity tensors.
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3TS
    use mod_iden_4O4T
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O3TS_I4O4T
        module procedure sum_I4O4T_I4O3TS
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O3TS_I4O4T
        module procedure sub_I4O4T_I4O3TS
    end interface
    
contains

    pure function sum_I4O3TS_I4O4T(I3S, I4) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}_{3T} + \mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I3S
            !! Scaled Type-3 identity (\(c \cdot \delta_{ij}\delta_{kl}\)).
        type(iden_4O4T), intent(in) :: I4
            !! Standard Type-4 identity (Symmetric \(\mathbf{I}^S\)).
        type(ten_3D4O3sym) :: res
            !! Resulting 21-component symmetric tensor.
        
        res%vals = (/1.0D0 + I3S%val, 1.0D0 + I3S%val, 1.0D0 + I3S%val, & ! Diagonals
                     0.5D0, 0.5D0, 0.5D0,                               & ! Shear diagonals
                     I3S%val, I3S%val, 0.0D0, 0.0D0, 0.0D0, I3S%val,    & ! Coupling
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          & ! Others
                     0.0D0, 0.0D0, 0.0D0/)
    end function sum_I4O3TS_I4O4T

    pure function sum_I4O4T_I4O3TS(I4, I3S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{4T} + c\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4T), intent(in) :: I4
        type(iden_4O3TS), intent(in) :: I3S
        type(ten_3D4O3sym) :: res
        res = sum_I4O3TS_I4O4T(I3S, I4)
    end function sum_I4O4T_I4O3TS

    pure function sub_I4O3TS_I4O4T(I3S, I4) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}_{3T} - \mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I3S
        type(iden_4O4T), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/I3S%val - 1.0D0, I3S%val - 1.0D0, I3S%val - 1.0D0, &
                     -0.5D0, -0.5D0, -0.5D0,                            &
                     I3S%val, I3S%val, 0.0D0, 0.0D0, 0.0D0, I3S%val,    &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O3TS_I4O4T

    pure function sub_I4O4T_I4O3TS(I4, I3S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{4T} - c\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4T), intent(in) :: I4
        type(iden_4O3TS), intent(in) :: I3S
        type(ten_3D4O3sym) :: res
        res%vals = (/1.0D0 - I3S%val, 1.0D0 - I3S%val, 1.0D0 - I3S%val, &
                     0.5D0, 0.5D0, 0.5D0,                               &
                     -I3S%val, -I3S%val, 0.0D0, 0.0D0, 0.0D0, -I3S%val, &
                      0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,         &
                      0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O4T_I4O3TS

end module mod_operator_I4O3TS_I4O4T