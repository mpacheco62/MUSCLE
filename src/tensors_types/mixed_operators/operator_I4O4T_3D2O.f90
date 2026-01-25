module mod_operator_I4O4T_3D2O
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4T
    use mod_ten_3D2Osym
    use mod_ten_3D2O
    implicit none
    private
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D2O_I4O4T
        module procedure ddot_I4O4T_3D2O
    end interface

    contains


    pure function ddot_3D2O_I4O4T(I2, a) result(res)
        implicit none
        class(iden_4O4T), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals(1) = a%vals(1)
        res%vals(2) = a%vals(5)
        res%vals(3) = a%vals(9)
        res%vals(4) = (a%vals(2)+a%vals(4))/2.0D0
        res%vals(5) = (a%vals(6)+a%vals(8))/2.0D0
        res%vals(6) = (a%vals(3)+a%vals(7))/2.0D0
    end function ddot_3D2O_I4O4T

    pure function ddot_I4O4T_3D2O(a, I2) result(res)
        implicit none
        class(iden_4O4T), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals(1) = a%vals(1)
        res%vals(2) = a%vals(5)
        res%vals(3) = a%vals(9)
        res%vals(4) = (a%vals(2)+a%vals(4))/2.0D0
        res%vals(5) = (a%vals(6)+a%vals(8))/2.0D0
        res%vals(6) = (a%vals(3)+a%vals(7))/2.0D0
    end function ddot_I4O4T_3D2O

end module mod_operator_I4O4T_3D2O