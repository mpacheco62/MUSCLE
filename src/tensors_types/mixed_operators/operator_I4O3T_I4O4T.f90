module mod_operator_I4O3T_I4O4T
    !! Module mod_operator_I4O3T_I4O4T
    !! =====================================
    !!
    !! Defines mixed algebraic operations between two different types of 3D 
    !! fourth-order identity tensors:
    !! 1. Type 3 (`iden_4O3T`): \(\mathbb{I}_{3T} = \delta_{ij}\delta_{kl}\)
    !! 2. Type 4 (`iden_4O4T`): \(\mathbb{I}_{4T} = \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\)
    !!
    !! The result of these operations is a fully symmetric 4th-order tensor 
    !! (`ten_3D4O3sym`) stored in the 21-component compressed Voigt format.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of Type-3 and Type-4 identity tensors.
    !! - `-` : Subtraction between Type-3 and Type-4 identity tensors.
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3T
    use mod_iden_4O4T
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O3T_I4O4T
        module procedure sum_I4O4T_I4O3T
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O3T_I4O4T
        module procedure sub_I4O4T_I4O3T
    end interface
    
contains

    pure function sum_I4O3T_I4O4T(I3, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{3T} + \mathbb{I}_{4T}\).
        !! Result indices in Voigt: Diagonals = 2.0, Shears = 0.5, Coupling = 1.0.
        implicit none
        type(iden_4O3T), intent(in) :: I3
            !! Standard Type-3 identity (\(\delta_{ij}\delta_{kl}\)).
        type(iden_4O4T), intent(in) :: I4
            !! Standard Type-4 identity (Symmetric \(\mathbf{I}^S\)).
        type(ten_3D4O3sym) :: res
            !! Resulting 21-component symmetric tensor.
        
        res%vals = (/2.0D0, 2.0D0, 2.0D0, 0.5D0, 0.5D0, 0.5D0, &
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                     0.0D0, 0.0D0, 0.0D0/)
        ! Mapping: (1-3: diag, 4-6: shear, 7-8: coupling, 12: coupling, others: 0)
    end function sum_I4O3T_I4O4T

    pure function sum_I4O4T_I4O3T(I4, I3) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{4T} + \mathbb{I}_{3T}\).
        implicit none
        type(iden_4O4T), intent(in) :: I4
        type(iden_4O3T), intent(in) :: I3
        type(ten_3D4O3sym) :: res
        res%vals = (/2.0D0, 2.0D0, 2.0D0, 0.5D0, 0.5D0, 0.5D0, &
                     1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, &
                     0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sum_I4O4T_I4O3T

    pure function sub_I4O3T_I4O4T(I3, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{3T} - \mathbb{I}_{4T}\).
        !! Result indices in Voigt: Diagonals = 0.0, Shears = -0.5, Coupling = 1.0.
        implicit none
        type(iden_4O3T), intent(in) :: I3
        type(iden_4O4T), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/0.0D0, 0.0D0, 0.0D0, -0.5D0, -0.5D0, -0.5D0, &
                     1.0D0, 1.0D0, 0.0D0,  0.0D0,  0.0D0,  1.0D0, &
                     0.0D0, 0.0D0, 0.0D0,  0.0D0,  0.0D0,  0.0D0, &
                     0.0D0, 0.0D0, 0.0D0/)
    end function sub_I4O3T_I4O4T

    pure function sub_I4O4T_I4O3T(I4, I3) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}_{4T} - \mathbb{I}_{3T}\).
        !! Result indices in Voigt: Diagonals = 0.0, Shears = 0.5, Coupling = -1.0.
        implicit none
        type(iden_4O4T), intent(in) :: I4
        type(iden_4O3T), intent(in) :: I3
        type(ten_3D4O3sym) :: res
        res%vals = (/ 0.0D0,  0.0D0, 0.0D0, 0.5D0, 0.5D0,  0.5D0, &
                     -1.0D0, -1.0D0, 0.0D0, 0.0D0, 0.0D0, -1.0D0, &
                      0.0D0,  0.0D0, 0.0D0, 0.0D0, 0.0D0,  0.0D0, &
                      0.0D0,  0.0D0, 0.0D0/)
    end function sub_I4O4T_I4O3T

end module mod_operator_I4O3T_I4O4T