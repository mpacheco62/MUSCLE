module mod_operator_I4O4T_3D4O3sym
    !! Module mod_operator_I4O4T_3D4O3sym
    !! ==================================
    !!
    !! Defines mixed algebraic operations between the standard 3D fourth-order 
    !! symmetric identity tensor (`iden_4O4T`, \(\mathbb{I}^S\)) and 3D 
    !! fully symmetric fourth-order tensors (`ten_3D4O3sym`).
    !!
    !! The `ten_3D4O3sym` type uses a compressed 21-component storage.
    !! In this format, the identity \(\mathbb{I}^S\) affects only the first 
    !! six components (the diagonal of the 6x6 Voigt matrix):
    !! - Components 1, 2, 3 (Normal): + 1.0
    !! - Components 4, 5, 6 (Shear):  + 0.5
    !!
    !! This module provides overloads for addition and subtraction.

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4T
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O4T_3D4O3sym
        module procedure sum_3D4O3sym_I4O4T
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O4T_3D4O3sym
        module procedure sub_3D4O3sym_I4O4T
    end interface

contains

    ! =========================================================================
    ! ADDITION
    ! =========================================================================

    pure function sum_I4O4T_3D4O3sym(I4, a) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}^S + \mathbf{A}\).
        implicit none
        type(iden_4O4T),    intent(in) :: I4
            !! Standard fourth-order symmetric identity tensor.
        type(ten_3D4O3sym), intent(in) :: a
            !! Fully symmetric fourth-order tensor (21 components).
        type(ten_3D4O3sym)             :: res

        res%vals = a%vals
        ! Add 1.0 to normal diagonals (11, 22, 33)
        res%vals(1:3) = res%vals(1:3) + 1.0D0
        ! Add 0.5 to shear diagonals (44, 55, 66 in Voigt-notation indices)
        res%vals(4:6) = res%vals(4:6) + 0.5D0
    end function sum_I4O4T_3D4O3sym

    pure function sum_3D4O3sym_I4O4T(a, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} + \mathbb{I}^S\).
        implicit none
        type(ten_3D4O3sym), intent(in) :: a
        type(iden_4O4T),    intent(in) :: I4
        type(ten_3D4O3sym)             :: res
        
        ! Delegate to ensure identical behavior and enable compiler inlining
        res = sum_I4O4T_3D4O3sym(I4, a)
    end function sum_3D4O3sym_I4O4T

    ! =========================================================================
    ! SUBTRACTION
    ! =========================================================================

    pure function sub_I4O4T_3D4O3sym(I4, a) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}^S - \mathbf{A}\).
        implicit none
        type(iden_4O4T),    intent(in) :: I4
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym)             :: res

        ! Compute -A
        res%vals = -a%vals
        ! Add Identity components: res = I + (-A)
        res%vals(1:3) = res%vals(1:3) + 1.0D0
        res%vals(4:6) = res%vals(4:6) + 0.5D0
    end function sub_I4O4T_3D4O3sym

    pure function sub_3D4O3sym_I4O4T(a, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} - \mathbb{I}^S\).
        implicit none
        type(ten_3D4O3sym), intent(in) :: a
        type(iden_4O4T),    intent(in) :: I4
        type(ten_3D4O3sym)             :: res

        res%vals = a%vals
        ! Subtract identity components from diagonal
        res%vals(1:3) = res%vals(1:3) - 1.0D0
        res%vals(4:6) = res%vals(4:6) - 0.5D0
    end function sub_3D4O3sym_I4O4T

end module mod_operator_I4O4T_3D4O3sym