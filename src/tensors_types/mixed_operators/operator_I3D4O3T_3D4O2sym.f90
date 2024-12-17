module mod_operator_I3D4O3T_3D4O2sym
    use, intrinsic :: iso_fortran_env
    use mod_iden_3D4O3T
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
        module procedure sum_3D4O2sym_I3D4O3T
        module procedure sum_I3D4O3T_3D4O2sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O2sym_I3D4O3T
        module procedure sub_I3D4O3T_3D4O2sym
    end interface
    
    contains

    pure module function sum_I3D4O3T_3D4O2sym(I2, a) result(res)
        implicit none
        class(iden_3D4O3T), intent(in) :: I2
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        call res%init(a%vals)
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1D0
    end function sum_I3D4O3T_3D4O2sym

    pure module function sum_3D4O2sym_I3D4O3T(a, I2) result(res)
        implicit none
        class(iden_3D4O3T), intent(in) :: I2
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        call res%init(a%vals)
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1D0
    end function sum_3D4O2sym_I3D4O3T

    pure module function sub_I3D4O3T_3D4O2sym(I2, a) result(res)
        implicit none
        class(iden_3D4O3T), intent(in) :: I2
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        call res%init(-a%vals)
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) + 1D0
    end function sub_I3D4O3T_3D4O2sym

    pure module function sub_3D4O2sym_I3D4O3T(a, I2) result(res)
        implicit none
        class(iden_3D4O3T), intent(in) :: I2
        type(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        call res%init(a%vals)
        res%vals(1:3, 1:3) = res%vals(1:3, 1:3) - 1D0
    end function sub_3D4O2sym_I3D4O3T
end module mod_operator_I3D4O3T_3D4O2sym