module mod_operator_I4O3TS_I4O4TS
    !! Module mod_operator_I4O3TS_I4O4TS
    !! =======================================
    !!
    !! Defines mixed algebraic operations between two scaled 3D fourth-order 
    !! identity tensors:
    !! 1. Scaled Type 3 (`iden_4O3TS`): \(\mathbb{I}_{3TS} = c_1 \cdot \delta_{ij}\delta_{kl}\)
    !! 2. Scaled Type 4 (`iden_4O4TS`): \(\mathbb{I}_{4TS} = c_2 \cdot \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\)
    !!
    !! The result is a fully symmetric 4th-order tensor (`ten_3D4O3sym`) stored 
    !! in the 21-component compressed Voigt format.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of two different scaled 4th-order identity tensors.
    !! - `-` : Subtraction between two different scaled 4th-order identity tensors.
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3TS
    use mod_iden_4O4TS
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O3TS_I4O4TS
        module procedure sum_I4O4TS_I4O3TS
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O3TS_I4O4TS
        module procedure sub_I4O4TS_I4O3TS
    end interface
    
contains

    pure function sum_I4O3TS_I4O4TS(I3S, I4S) result(res)
        !! Computes \(\mathbb{res} = c_1\mathbb{I}_{3T} + c_2\mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I3S
            !! Scaled Type-3 identity (\(c_1 \cdot \delta_{ij}\delta_{kl}\)).
        type(iden_4O4TS), intent(in) :: I4S
            !! Scaled Type-4 identity (\(c_2 \cdot \mathbf{I}^S\)).
        type(ten_3D4O3sym) :: res
            !! Resulting 21-component symmetric tensor.
        
        res%vals = (/I4S%val + I3S%val, I4S%val + I3S%val, I4S%val + I3S%val, & ! Diagonals
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val,       & ! Shear diagonals
                     I3S%val, I3S%val, 0.0D0, 0.0D0, 0.0D0, I3S%val,          & ! Normal coupling
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,                & ! Others
                     0.0D0, 0.0D0, 0.0D0/)
    end function sum_I4O3TS_I4O4TS

    pure function sum_I4O4TS_I4O3TS(I4S, I3S) result(res)
        !! Computes \(\mathbb{res} = c_2\mathbb{I}_{4T} + c_1\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
        type(iden_4O3TS), intent(in) :: I3S
        type(ten_3D4O3sym) :: res
        res = sum_I4O3TS_I4O4TS(I3S, I4S)
    end function sum_I4O4TS_I4O3TS

    pure function sub_I4O3TS_I4O4TS(I3S, I4S) result(res)
        !! Computes \(\mathbb{res} = c_1\mathbb{I}_{3T} - c_2\mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3TS), intent(in) :: I3S
        type(iden_4O4TS), intent(in) :: I4S
        type(ten_3D4O3sym) :: res
        res%vals = (/I3S%val - I4S%val, I3S%val - I4S%val, I3S%val - I4S%val, &
                     -0.5D0 * I4S%val, -0.5D0 * I4S%val, -0.5D0 * I4S%val,    &
                     I3S%val, I3S%val, 0.0D0, 0.0D0, 0.0D0, I3S%val,          &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,                &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O3TS_I4O4TS

    pure function sub_I4O4TS_I4O3TS(I4S, I3S) result(res)
        !! Computes \(\mathbb{res} = c_2\mathbb{I}_{4T} - c_1\mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
        type(iden_4O3TS), intent(in) :: I3S
        type(ten_3D4O3sym) :: res
        res%vals = (/I4S%val - I3S%val, I4S%val - I3S%val, I4S%val - I3S%val, &
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val,       &
                     -I3S%val, -I3S%val, 0.0D0, 0.0D0, 0.0D0, -I3S%val,       &
                      0.0D0,  0.0D0, 0.0D0, 0.0D0, 0.0D0,  0.0D0,             &
                      0.0D0,  0.0D0, 0.0D0/)
    end function sub_I4O4TS_I4O3TS

end module mod_operator_I4O3TS_I4O4TS