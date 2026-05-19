module mod_operator_I4O4TS_3D4O2sym
    !! Module mod_operator_I4O4TS_3D4O2sym
    !! ===================================
    !!
    !! Defines mixed algebraic operations between the scaled 3D fourth-order 
    !! symmetric identity tensor (`iden_4O4TS`, \(c\mathbb{I}^S\)) and 
    !! 3D fourth-order tensors with minor symmetries (`ten_3D4O2sym`).
    !!
    !! In the 6x6 Voigt matrix representation, the scaled symmetric identity 
    !! \(c\mathbb{I}^S\) corresponds to:
    !! \[ \text{diag}(c, c, c, 0.5c, 0.5c, 0.5c) \]
    !!
    !! This module provides overloads for addition, subtraction, and 
    !! explicit assignment.

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4TS
    use mod_ten_3D4O2sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O4TS_3D4O2sym
        module procedure sum_3D4O2sym_I4O4TS
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O4TS_3D4O2sym
        module procedure sub_3D4O2sym_I4O4TS
    end interface
    
    public :: assignment(=)
    interface assignment(=)
        module procedure assign_3D4O2sym_I4O4TS
    end interface

contains

    ! =========================================================================
    ! ADDITION
    ! =========================================================================

    pure function sum_I4O4TS_3D4O2sym(I4S, a) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}^S + \mathbb{A}\).
        implicit none
        type(iden_4O4TS),  intent(in) :: I4S
            !! Scaled fourth-order symmetric identity tensor.
        type(ten_3D4O2sym), intent(in) :: a
            !! Fourth-order tensor with minor symmetries (6x6 matrix).
        type(ten_3D4O2sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = a%vals
        ! Add components to the Voigt diagonal
        res%vals(1,1) = res%vals(1,1) + c
        res%vals(2,2) = res%vals(2,2) + c
        res%vals(3,3) = res%vals(3,3) + c
        res%vals(4,4) = res%vals(4,4) + half_c
        res%vals(5,5) = res%vals(5,5) + half_c
        res%vals(6,6) = res%vals(6,6) + half_c
    end function sum_I4O4TS_3D4O2sym

    pure function sum_3D4O2sym_I4O4TS(a, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{A} + c\mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
        type(iden_4O4TS),  intent(in) :: I4S
        type(ten_3D4O2sym)             :: res
        
        ! Commutative property: delegate to the primary implementation
        res = sum_I4O4TS_3D4O2sym(I4S, a)
    end function sum_3D4O2sym_I4O4TS

    ! =========================================================================
    ! SUBTRACTION
    ! =========================================================================

    pure function sub_I4O4TS_3D4O2sym(I4S, a) result(res)
        !! Computes \(\mathbb{res} = c\mathbb{I}^S - \mathbb{A}\).
        implicit none
        type(iden_4O4TS),  intent(in) :: I4S
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = -a%vals
        res%vals(1,1) = res%vals(1,1) + c
        res%vals(2,2) = res%vals(2,2) + c
        res%vals(3,3) = res%vals(3,3) + c
        res%vals(4,4) = res%vals(4,4) + half_c
        res%vals(5,5) = res%vals(5,5) + half_c
        res%vals(6,6) = res%vals(6,6) + half_c
    end function sub_I4O4TS_3D4O2sym

    pure function sub_3D4O2sym_I4O4TS(a, I4S) result(res)
        !! Computes \(\mathbb{res} = \mathbb{A} - c\mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(in) :: a
        type(iden_4O4TS),  intent(in) :: I4S
        type(ten_3D4O2sym)             :: res
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        res%vals = a%vals
        res%vals(1,1) = res%vals(1,1) - c
        res%vals(2,2) = res%vals(2,2) - c
        res%vals(3,3) = res%vals(3,3) - c
        res%vals(4,4) = res%vals(4,4) - half_c
        res%vals(5,5) = res%vals(5,5) - half_c
        res%vals(6,6) = res%vals(6,6) - half_c
    end function sub_3D4O2sym_I4O4TS

    ! =========================================================================
    ! ASSIGNMENT
    ! =========================================================================

    pure subroutine assign_3D4O2sym_I4O4TS(a, I4S)
        !! Directly assigns a scaled symmetric identity to a 4th-order tensor: \(\mathbb{A} = c\mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(out) :: a
        type(iden_4O4TS),  intent(in)  :: I4S
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        a%vals = 0.0D0
        a%vals(1,1) = c
        a%vals(2,2) = c
        a%vals(3,3) = c
        a%vals(4,4) = half_c
        a%vals(5,5) = half_c
        a%vals(6,6) = half_c
    end subroutine assign_3D4O2sym_I4O4TS

end module mod_operator_I4O4TS_3D4O2sym