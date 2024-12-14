module mod_operator_I3D4O3TMod_3D4O3sym
    use, intrinsic :: iso_fortran_env
    use mod_iden_3D4O3TMod
    use mod_ten_3D4O3sym
    implicit none
    private
    
    !
    !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
    !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
    !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
    !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
    !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
    !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O3sym_I3D4O3TMod
        module procedure sum_I3D4O3TMod_3D4O3sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O3sym_I3D4O3TMod
        module procedure sub_I3D4O3TMod_3D4O3sym
    end interface
    
    ! public :: operator(.ddot.)
    ! interface operator (.ddot.)
    !     module procedure ddot_I3D4O3TMod_3D4O3sym
    !     module procedure ddot_3D4O3sym_I3D4O3TMod
    ! end interface

    contains

    ! pure module function ddot_I3D4O3TMod_3D4O3sym(I2, b) result(res)
    !     implicit none
    !     class(iden_3D4O3TMod), intent(in) :: I2
    !     class(ten_3D4O3sym), intent(in) :: b
    !     real(real64) :: res
    !     res = b%vals(1) + b%vals(2) + b%vals(3)
    ! end function ddot_I3D4O3TMod_3D4O3sym

    ! pure module function ddot_3D4O3sym_I3D4O3TMod(b, I2) result(res)
    !     implicit none
    !     class(iden_3D4O3TMod), intent(in) :: I2
    !     class(ten_3D4O3sym), intent(in) :: b
    !     real(real64) :: res
    !     res = b%vals(1) + b%vals(2) + b%vals(3)
    ! end function ddot_3D4O3sym_I3D4O3TMod

    pure module function sum_I3D4O3TMod_3D4O3sym(I2, a) result(res)
        implicit none
        class(iden_3D4O3TMod), intent(in) :: I2
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        res%vals( 1: 3) = a%vals( 1: 3) + I2%val
        res%vals( 4: 6) = a%vals( 4: 6)
        res%vals( 7: 8) = a%vals( 7: 8) + I2%val
        res%vals( 9:11) = a%vals( 9:11)
        res%vals(   12) = a%vals(   12) + I2%val
        res%vals(13:21) = a%vals(13:21)
    end function sum_I3D4O3TMod_3D4O3sym

    pure module function sum_3D4O3sym_I3D4O3TMod(a, I2) result(res)
        implicit none
        class(iden_3D4O3TMod), intent(in) :: I2
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        res%vals( 1: 3) = a%vals( 1: 3) + I2%val
        res%vals( 4: 6) = a%vals( 4: 6)
        res%vals( 7: 8) = a%vals( 7: 8) + I2%val
        res%vals( 9:11) = a%vals( 9:11)
        res%vals(   12) = a%vals(   12) + I2%val
        res%vals(13:21) = a%vals(13:21)
    end function sum_3D4O3sym_I3D4O3TMod

    pure module function sub_I3D4O3TMod_3D4O3sym(I2, a) result(res)
        implicit none
        class(iden_3D4O3TMod), intent(in) :: I2
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        res%vals( 1: 3) = -a%vals( 1: 3) + I2%val
        res%vals( 4: 6) = -a%vals( 4: 6)
        res%vals( 7: 8) = -a%vals( 7: 8) + I2%val
        res%vals( 9:11) = -a%vals( 9:11)
        res%vals(   12) = -a%vals(   12) + I2%val
        res%vals(13:21) = -a%vals(13:21)
    end function sub_I3D4O3TMod_3D4O3sym

    pure module function sub_3D4O3sym_I3D4O3TMod(a, I2) result(res)
        implicit none
        class(iden_3D4O3TMod), intent(in) :: I2
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        res%vals( 1: 3) = a%vals( 1: 3) - I2%val
        res%vals( 4: 6) = a%vals( 4: 6)
        res%vals( 7: 8) = a%vals( 7: 8) - I2%val
        res%vals( 9:11) = a%vals( 9:11)
        res%vals(   12) = a%vals(   12) - I2%val
        res%vals(13:21) = a%vals(13:21)
    end function sub_3D4O3sym_I3D4O3TMod
end module mod_operator_I3D4O3TMod_3D4O3sym