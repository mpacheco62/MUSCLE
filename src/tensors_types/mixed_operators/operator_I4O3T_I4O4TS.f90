module mod_operator_I4O3T_I4O4TS
    !! Module mod_operator_I4O3T_I4O4TS
    !! ======================================
    !!
    !! Defines mixed algebraic operations between the standard 3D fourth-order 
    !! identity tensor of type 3 (`iden_4O3T`, \(\delta_{ij}\delta_{kl}\)) and 
    !! the scaled 3D fourth-order symmetric identity tensor (`iden_4O4TS`, \(c \cdot \mathbf{I}^S\)).
    !!
    !! The result is a fully symmetric 4th-order tensor (`ten_3D4O3sym`) stored 
    !! in the 21-component compressed Voigt format.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of Type-3 and Scaled Type-4 identity tensors.
    !! - `-` : Subtraction between Type-3 and Scaled Type-4 identity tensors.
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3T
    use mod_iden_4O4TS
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O3T_I4O4TS
        module procedure sum_I4O4TS_I4O3T
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O3T_I4O4TS
        module procedure sub_I4O4TS_I4O3T
    end interface
    
contains

    pure function sum_I4O3T_I4O4TS(I3, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{3T} + c\mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3T), intent(in) :: I3
            !! Standard Type-3 identity (\(\delta_{ij}\delta_{kl}\)).
        type(iden_4O4TS), intent(in) :: I4S
            !! Scaled Type-4 identity (\(c\mathbf{I}^S\)).
        type(ten_3D4O3sym) :: res
            !! Resulting 21-component symmetric tensor.
        
        res%vals = (/I4S%val + 1.0D0, I4S%val + 1.0D0, I4S%val + 1.0D0, & ! Diagonals
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val, & ! Shears
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0,          & ! Normal coupling
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          & ! Mixed
                     0.0D0, 0.0D0, 0.0D0/)                                ! Mixed
    end function sum_I4O3T_I4O4TS

    pure function sum_I4O4TS_I4O3T(I4S, I3) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}_{4T} + \mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
        type(iden_4O3T), intent(in) :: I3
        type(ten_3D4O3sym) :: res
        res%vals = (/I4S%val + 1.0D0, I4S%val + 1.0D0, I4S%val + 1.0D0, &
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val, &
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0,          &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sum_I4O4TS_I4O3T

    pure function sub_I4O3T_I4O4TS(I3, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{3T} - c\mathbb{I}_{4T}\).
        implicit none
        type(iden_4O3T), intent(in) :: I3
        type(iden_4O4TS), intent(in) :: I4S
        type(ten_3D4O3sym) :: res
        res%vals = (/1.0D0 - I4S%val, 1.0D0 - I4S%val, 1.0D0 - I4S%val, &
                     -0.5D0 * I4S%val, -0.5D0 * I4S%val, -0.5D0 * I4S%val, &
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0,          &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,          &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O3T_I4O4TS

    pure function sub_I4O4TS_I4O3T(I4S, I3) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}_{4T} - \mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
        type(iden_4O3T), intent(in) :: I3
        type(ten_3D4O3sym) :: res
        res%vals = (/I4S%val - 1.0D0, I4S%val - 1.0D0, I4S%val - 1.0D0, &
                     0.5D0 * I4S%val, 0.5D0 * I4S%val, 0.5D0 * I4S%val, &
                     -1.0D0, -1.0D0, 0.0D0, 0.0D0, 0.0D0, -1.0D0,       &
                      0.0D0,  0.0D0, 0.0D0, 0.0D0, 0.0D0,  0.0D0,       &
                      0.0D0,  0.0D0, 0.0D0/)
    end function sub_I4O4TS_I4O3T

end module mod_operator_I4O3T_I4O4TS