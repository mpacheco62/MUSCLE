module mod_operator_I4O4TS_3D2Osym
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4TS
    use mod_ten_3D2Osym
    implicit none
    private
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D2Osym_I4O4TS
        module procedure ddot_I4O4TS_3D2Osym
    end interface

    contains


    pure function ddot_3D2Osym_I4O4TS(I2, a) result(res)
        implicit none
        class(iden_4O4TS), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals*I2%val
    end function ddot_3D2Osym_I4O4TS

    pure function ddot_I4O4TS_3D2Osym(a, I2) result(res)
        implicit none
        class(iden_4O4TS), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals*I2%val
    end function ddot_I4O4TS_3D2Osym

end module mod_operator_I4O4TS_3D2Osym