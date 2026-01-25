module mod_operator_I4O4TS_3D2O
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4TS
    use mod_ten_3D2Osym
    use mod_ten_3D2O
    implicit none
    private
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D2O_I4O4TS
        module procedure ddot_I4O4TS_3D2O
    end interface

    contains


    pure function ddot_3D2O_I4O4TS(I2, a) result(res)
        implicit none
        class(iden_4O4TS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals(1) = a%vals(1)*I2%val
        res%vals(2) = a%vals(5)*I2%val
        res%vals(3) = a%vals(9)*I2%val
        res%vals(4) = (a%vals(2)+a%vals(4))/2.0D0*I2%val
        res%vals(5) = (a%vals(6)+a%vals(8))/2.0D0*I2%val
        res%vals(6) = (a%vals(3)+a%vals(7))/2.0D0*I2%val
    end function ddot_3D2O_I4O4TS

    pure function ddot_I4O4TS_3D2O(a, I2) result(res)
        implicit none
        class(iden_4O4TS), intent(in) :: I2
        type(ten_3D2O), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals(1) = a%vals(1)*I2%val
        res%vals(2) = a%vals(5)*I2%val
        res%vals(3) = a%vals(9)*I2%val
        res%vals(4) = (a%vals(2)+a%vals(4))/2.0D0*I2%val
        res%vals(5) = (a%vals(6)+a%vals(8))/2.0D0*I2%val
        res%vals(6) = (a%vals(3)+a%vals(7))/2.0D0*I2%val
    end function ddot_I4O4TS_3D2O
end module mod_operator_I4O4TS_3D2O