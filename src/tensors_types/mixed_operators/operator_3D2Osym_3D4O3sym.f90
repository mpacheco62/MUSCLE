module mod_operator_3D2Osym_3D4O3sym
    use, intrinsic :: iso_fortran_env
    use mod_ten_3D2Osym
    use mod_ten_3D4O3sym
    implicit none
    private


    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D4O3sym_3D2Osym
        module procedure ddot_3D2Osym_3D4O3sym
    end interface

    public :: operator(.tdot.)
    interface operator (.tdot.)
        module procedure tdot_3D2Osym_3D2Osym
    end interface

    contains

    pure module function ddot_3D4O3sym_3D2Osym(a, b) result(res)
        !
        !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
        !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
        !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
        !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
        !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
        !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |
        implicit none
        class(ten_3D4O3sym), intent(in) :: a
        class(ten_3D2Osym), intent(in) :: b
        type(ten_3D2Osym) :: res
        
        real(real64) :: a1111b11, a1122b22, a1133b33, a1112b12, a1123b23, a1113b13
        real(real64) :: a2211b11, a2222b22, a2233b33, a2212b12, a2223b23, a2213b13
        real(real64) :: a3311b11, a3322b22, a3333b33, a3312b12, a3323b23, a3313b13
        real(real64) :: a1112b11, a2212b22, a3312b33, a1212b12, a1223b23, a1213b13
        real(real64) :: a1123b11, a2223b22, a3323b33, a2312b12, a2323b23, a2313b13
        real(real64) :: a1113b11, a2213b22, a3313b33, a1312b12, a1323b23, a1313b13

        a1111b11 = a%vals( 1)*b%vals(1)
        a1122b22 = a%vals( 7)*b%vals(2)
        a1133b33 = a%vals(12)*b%vals(3)
        a1112b12 = a%vals(16)*b%vals(4)
        a1123b23 = a%vals(19)*b%vals(5)
        a1113b13 = a%vals(21)*b%vals(6)

        a2211b11 = a%vals( 7)*b%vals(1)
        a2222b22 = a%vals( 2)*b%vals(2)
        a2233b33 = a%vals( 8)*b%vals(3)
        a2212b12 = a%vals(13)*b%vals(4)
        a2223b23 = a%vals(17)*b%vals(5)
        a2213b13 = a%vals(20)*b%vals(6)

        a3311b11 = a%vals(12)*b%vals(1)
        a3322b22 = a%vals( 8)*b%vals(2)
        a3333b33 = a%vals( 3)*b%vals(3)
        a3312b12 = a%vals( 9)*b%vals(4)
        a3323b23 = a%vals(14)*b%vals(5)
        a3313b13 = a%vals(18)*b%vals(6)

        a1112b11 = a%vals(16)*b%vals(1)
        a2212b22 = a%vals(13)*b%vals(2)
        a3312b33 = a%vals( 9)*b%vals(3)
        a1212b12 = a%vals( 4)*b%vals(4)
        a1223b23 = a%vals(10)*b%vals(5)
        a1213b13 = a%vals(15)*b%vals(6)

        a1123b11 = a%vals(19)*b%vals(1)
        a2223b22 = a%vals(17)*b%vals(2)
        a3323b33 = a%vals(14)*b%vals(3)
        a2312b12 = a%vals(10)*b%vals(4)
        a2323b23 = a%vals( 5)*b%vals(5)
        a2313b13 = a%vals(11)*b%vals(6)

        a1113b11 = a%vals(21)*b%vals(1)
        a2213b22 = a%vals(20)*b%vals(2)
        a3313b33 = a%vals(18)*b%vals(3)
        a1312b12 = a%vals(15)*b%vals(4)
        a1323b23 = a%vals(11)*b%vals(5)
        a1313b13 = a%vals( 6)*b%vals(6)
        
        res%vals(1) = a1111b11 + a1122b22 + a1133b33 + 2*a1112b12 + 2*a1123b23 + 2*a1113b13
        res%vals(2) = a2211b11 + a2222b22 + a2233b33 + 2*a2212b12 + 2*a2223b23 + 2*a2213b13
        res%vals(3) = a3311b11 + a3322b22 + a3333b33 + 2*a3312b12 + 2*a3323b23 + 2*a3313b13
        res%vals(4) = a1112b11 + a2212b22 + a3312b33 + 2*a1212b12 + 2*a1223b23 + 2*a1213b13
        res%vals(5) = a1123b11 + a2223b22 + a3323b33 + 2*a2312b12 + 2*a2323b23 + 2*a2313b13
        res%vals(6) = a1113b11 + a2213b22 + a3313b33 + 2*a1312b12 + 2*a1323b23 + 2*a1313b13
        
    end function ddot_3D4O3sym_3D2Osym

    pure module function ddot_3D2Osym_3D4O3sym(b, a) result(res)
        !
        !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
        !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
        !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
        !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
        !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
        !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |
        implicit none
        class(ten_3D2Osym), intent(in) :: b
        class(ten_3D4O3sym), intent(in) :: a
        type(ten_3D2Osym) :: res
        
        real(real64) :: a1111b11, a2211b22, a3311b33, a1211b12, a2311b23, a1311b13
        real(real64) :: a1122b11, a2222b22, a3322b33, a1222b12, a2322b23, a1322b13
        real(real64) :: a1133b11, a2233b22, a3333b33, a1233b12, a2333b23, a1333b13
        real(real64) :: a1112b11, a2212b22, a3312b33, a1212b12, a2312b23, a1312b13
        real(real64) :: a1123b11, a2223b22, a3323b33, a1223b12, a2323b23, a1323b13
        real(real64) :: a1113b11, a2213b22, a3313b33, a1213b12, a2313b23, a1313b13

        a1111b11 = a%vals( 1)*b%vals(1)
        a2211b22 = a%vals( 7)*b%vals(2)
        a3311b33 = a%vals(12)*b%vals(3)
        a1211b12 = a%vals(16)*b%vals(4)
        a2311b23 = a%vals(19)*b%vals(5)
        a1311b13 = a%vals(21)*b%vals(6)

        a1122b11 = a%vals( 7)*b%vals(1)
        a2222b22 = a%vals( 2)*b%vals(2)
        a3322b33 = a%vals( 8)*b%vals(3)
        a1222b12 = a%vals(13)*b%vals(4)
        a2322b23 = a%vals(17)*b%vals(5)
        a1322b13 = a%vals(20)*b%vals(6)

        a1133b11 = a%vals(12)*b%vals(1)
        a2233b22 = a%vals( 8)*b%vals(2)
        a3333b33 = a%vals( 3)*b%vals(3)
        a1233b12 = a%vals( 9)*b%vals(4)
        a2333b23 = a%vals(14)*b%vals(5)
        a1333b13 = a%vals(18)*b%vals(6)

        a1112b11 = a%vals(16)*b%vals(1)
        a2212b22 = a%vals(13)*b%vals(2)
        a3312b33 = a%vals( 9)*b%vals(3)
        a1212b12 = a%vals( 4)*b%vals(4)
        a2312b23 = a%vals(10)*b%vals(5)
        a1312b13 = a%vals(15)*b%vals(6)

        a1123b11 = a%vals(19)*b%vals(1)
        a2223b22 = a%vals(17)*b%vals(2)
        a3323b33 = a%vals(14)*b%vals(3)
        a1223b12 = a%vals(10)*b%vals(4)
        a2323b23 = a%vals( 5)*b%vals(5)
        a1323b13 = a%vals(11)*b%vals(6)

        a1113b11 = a%vals(21)*b%vals(1)
        a2213b22 = a%vals(20)*b%vals(2)
        a3313b33 = a%vals(18)*b%vals(3)
        a1213b12 = a%vals(15)*b%vals(4)
        a2313b23 = a%vals(11)*b%vals(5)
        a1313b13 = a%vals( 6)*b%vals(6)
        
        res%vals(1) = a1111b11 + a2211b22 + a3311b33 + 2*a1211b12 + 2*a2311b23 + 2*a1311b13
        res%vals(2) = a1122b11 + a2222b22 + a3322b33 + 2*a1222b12 + 2*a2322b23 + 2*a1322b13
        res%vals(3) = a1133b11 + a2233b22 + a3333b33 + 2*a1233b12 + 2*a2333b23 + 2*a1333b13
        res%vals(4) = a1112b11 + a2212b22 + a3312b33 + 2*a1212b12 + 2*a2312b23 + 2*a1312b13
        res%vals(5) = a1123b11 + a2223b22 + a3323b33 + 2*a1223b12 + 2*a2323b23 + 2*a1323b13
        res%vals(6) = a1113b11 + a2213b22 + a3313b33 + 2*a1213b12 + 2*a2313b23 + 2*a1313b13
        
    end function ddot_3D2Osym_3D4O3sym

    pure module function tdot_3D2Osym_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        type(ten_3D4O3sym) :: res
        res%vals( 1) = a%vals(1)*b%vals(1)
        res%vals( 7) = a%vals(1)*b%vals(2)
        res%vals(12) = a%vals(1)*b%vals(3)
        res%vals(16) = a%vals(1)*b%vals(4)
        res%vals(19) = a%vals(1)*b%vals(5)
        res%vals(21) = a%vals(1)*b%vals(6)

        res%vals( 2) = a%vals(2)*b%vals(2)
        res%vals( 8) = a%vals(2)*b%vals(3)
        res%vals(13) = a%vals(2)*b%vals(4)
        res%vals(17) = a%vals(2)*b%vals(5)
        res%vals(20) = a%vals(2)*b%vals(6)

        res%vals( 3) = a%vals(3)*b%vals(3)
        res%vals( 9) = a%vals(3)*b%vals(4)
        res%vals(14) = a%vals(3)*b%vals(5)
        res%vals(18) = a%vals(3)*b%vals(6)

        res%vals( 4) = a%vals(4)*b%vals(4)
        res%vals(10) = a%vals(4)*b%vals(5)
        res%vals(15) = a%vals(4)*b%vals(6)

        res%vals( 5) = a%vals(5)*b%vals(5)
        res%vals(11) = a%vals(5)*b%vals(6)
        
        res%vals( 6) = a%vals(6)*b%vals(6)
    end function tdot_3D2Osym_3D2Osym

end module mod_operator_3D2Osym_3D4O3sym