module mod_operator_I4O4T_3D4O2sym
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4T
    use mod_ten_3D4O2sym
    implicit none
    private
    
    !
    !  | (1111) (1122) (1133) (1112) (1123) (1113) |
    !  | (2211) (2222) (2233) (2212) (2223) (2213) |
    !  | (3311) (3322) (3333) (3312) (3323) (3313) |
    !  | (1211) (1222) (1233) (1212) (1223) (1213) |
    !  | (2311) (2322) (2333) (2312) (2323) (2313) |
    !  | (1311) (1322) (1333) (1312) (1323) (1313) |

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O2sym_I4O4T
        module procedure sum_I4O4T_3D4O2sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O2sym_I4O4T
        module procedure sub_I4O4T_3D4O2sym
    end interface
    
    public :: assignment(=)
    interface assignment(=)
        module procedure assign_3D4O2sym_I4O4T
    end interface

    contains


    pure function sum_I4O4T_3D4O2sym(I2, a) result(res)
        implicit none
        class(iden_4O4T), intent(in) :: I2
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        call res%init(a%vals)
        res%vals(1,1) = res%vals(1,1) + 1D0
        res%vals(2,2) = res%vals(2,2) + 1D0
        res%vals(3,3) = res%vals(3,3) + 1D0
        res%vals(4,4) = res%vals(4,4) + 0.5D0
        res%vals(5,5) = res%vals(5,5) + 0.5D0
        res%vals(6,6) = res%vals(6,6) + 0.5D0
    end function sum_I4O4T_3D4O2sym

    pure function sum_3D4O2sym_I4O4T(a, I2) result(res)
        implicit none
        class(iden_4O4T), intent(in) :: I2
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        call res%init(a%vals)
        res%vals(1,1) = res%vals(1,1) + 1D0
        res%vals(2,2) = res%vals(2,2) + 1D0
        res%vals(3,3) = res%vals(3,3) + 1D0
        res%vals(4,4) = res%vals(4,4) + 0.5D0
        res%vals(5,5) = res%vals(5,5) + 0.5D0
        res%vals(6,6) = res%vals(6,6) + 0.5D0
    end function sum_3D4O2sym_I4O4T

    pure function sub_I4O4T_3D4O2sym(I2, a) result(res)
        implicit none
        class(iden_4O4T), intent(in) :: I2
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        call res%init(-a%vals)
        res%vals(1,1) = res%vals(1,1) + 1D0
        res%vals(2,2) = res%vals(2,2) + 1D0
        res%vals(3,3) = res%vals(3,3) + 1D0
        res%vals(4,4) = res%vals(4,4) + 0.5D0
        res%vals(5,5) = res%vals(5,5) + 0.5D0
        res%vals(6,6) = res%vals(6,6) + 0.5D0
    end function sub_I4O4T_3D4O2sym

    pure function sub_3D4O2sym_I4O4T(a, I2) result(res)
        implicit none
        class(iden_4O4T), intent(in) :: I2
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        call res%init(a%vals)
        res%vals(1,1) = res%vals(1,1) - 1D0
        res%vals(2,2) = res%vals(2,2) - 1D0
        res%vals(3,3) = res%vals(3,3) - 1D0
        res%vals(4,4) = res%vals(4,4) - 0.5D0
        res%vals(5,5) = res%vals(5,5) - 0.5D0
        res%vals(6,6) = res%vals(6,6) - 0.5D0
    end function sub_3D4O2sym_I4O4T

    pure subroutine assign_3D4O2sym_I4O4T(a, I2)
        implicit none
        class(iden_4O4T), intent(in) :: I2
        type(ten_3D4O2sym), intent(out) :: a
        a%vals = 0D0
        a%vals(1,1) = 1D0
        a%vals(2,2) = 1D0
        a%vals(3,3) = 1D0
        a%vals(4,4) = 0.5D0
        a%vals(5,5) = 0.5D0
        a%vals(6,6) = 0.5D0
    end subroutine assign_3D4O2sym_I4O4T
end module mod_operator_I4O4T_3D4O2sym