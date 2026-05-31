module mod_operator_I4O4T_3D4O2sym
    !! Module mod_operator_I4O4T_3D4O2sym
    !! ==================================
    !!
    !! Defines mixed algebraic operations between the standard 3D fourth-order 
    !! symmetric identity tensor (`iden_4O4T`, \(\mathbb{I}^S\)) and 3D 
    !! fourth-order tensors with minor symmetries (`ten_3D4O2sym`).
    !!
    !! In Voigt notation (6x6 matrix), \(\mathbb{I}^S\) is represented as:
    !! diag(1, 1, 1, 0.5, 0.5, 0.5).
    !!
    !! This module provides overloads for addition, subtraction, and 
    !! explicit assignment.

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4T
    use mod_ten_3D4O2sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O4T_3D4O2sym
        module procedure sum_3D4O2sym_I4O4T
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O4T_3D4O2sym
        module procedure sub_3D4O2sym_I4O4T
    end interface
    
    public :: assignment(=)
    interface assignment(=)
        module procedure assign_3D4O2sym_I4O4T
    end interface

contains

    ! =========================================================================
    ! ADDITION
    ! =========================================================================

    pure function sum_I4O4T_3D4O2sym(I4, a) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}^S + \mathbf{A}\).
        implicit none
        type(iden_4O4T),   intent(in) :: I4
            !! Standard fourth-order symmetric identity tensor.
        type(ten_3D4O2sym), intent(in) :: a
            !! Fourth-order tensor with minor symmetries.
        type(ten_3D4O2sym)             :: res

        res%vals = a%vals
        ! Update Voigt diagonal components
        res%vals(1,1) = res%vals(1,1) + 1.0D0
        res%vals(2,2) = res%vals(2,2) + 1.0D0
        res%vals(3,3) = res%vals(3,3) + 1.0D0
        res%vals(4,4) = res%vals(4,4) + 0.5D0
        res%vals(5,5) = res%vals(5,5) + 0.5D0
        res%vals(6,6) = res%vals(6,6) + 0.5D0
    end function sum_I4O4T_3D4O2sym

    pure function sum_3D4O2sym_I4O4T(a, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} + \mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
        type(iden_4O4T),   intent(in) :: I4
        type(ten_3D4O2sym)             :: res
        
        ! Commutative property: delegate to the other implementation
        res = sum_I4O4T_3D4O2sym(I4, a)
    end function sum_3D4O2sym_I4O4T

    ! =========================================================================
    ! SUBTRACTION
    ! =========================================================================

    pure function sub_I4O4T_3D4O2sym(I4, a) result(res)
        !! Computes \(\mathbb{res} = \mathbb{I}^S - \mathbf{A}\).
        implicit none
        type(iden_4O4T),   intent(in) :: I4
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym)             :: res

        res%vals = -a%vals
        res%vals(1,1) = res%vals(1,1) + 1.0D0
        res%vals(2,2) = res%vals(2,2) + 1.0D0
        res%vals(3,3) = res%vals(3,3) + 1.0D0
        res%vals(4,4) = res%vals(4,4) + 0.5D0
        res%vals(5,5) = res%vals(5,5) + 0.5D0
        res%vals(6,6) = res%vals(6,6) + 0.5D0
    end function sub_I4O4T_3D4O2sym

    pure function sub_3D4O2sym_I4O4T(a, I4) result(res)
        !! Computes \(\mathbb{res} = \mathbf{A} - \mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
        type(iden_4O4T),   intent(in) :: I4
        type(ten_3D4O2sym)             :: res

        res%vals = a%vals
        res%vals(1,1) = res%vals(1,1) - 1.0D0
        res%vals(2,2) = res%vals(2,2) - 1.0D0
        res%vals(3,3) = res%vals(3,3) - 1.0D0
        res%vals(4,4) = res%vals(4,4) - 0.5D0
        res%vals(5,5) = res%vals(5,5) - 0.5D0
        res%vals(6,6) = res%vals(6,6) - 0.5D0
    end function sub_3D4O2sym_I4O4T

    ! =========================================================================
    ! ASSIGNMENT
    ! =========================================================================

    pure subroutine assign_3D4O2sym_I4O4T(a, I4)
        !! Assigns the symmetric identity tensor directly to a fourth-order tensor: \(\mathbf{A} = \mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(out) :: a
        type(iden_4O4T),   intent(in)  :: I4
        
        a%vals = 0.0D0
        a%vals(1,1) = 1.0D0
        a%vals(2,2) = 1.0D0
        a%vals(3,3) = 1.0D0
        a%vals(4,4) = 0.5D0
        a%vals(5,5) = 0.5D0
        a%vals(6,6) = 0.5D0
    end subroutine assign_3D4O2sym_I4O4T

end module mod_operator_I4O4T_3D4O2sym