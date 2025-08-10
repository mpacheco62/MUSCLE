module mod_operator_2D2Osym_2D4O3sym
    use, intrinsic :: iso_fortran_env
    use mod_ten_2D2Osym
    use mod_ten_2D4O3sym
    implicit none
    private


    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_2D4O3sym_2D2Osym
        module procedure ddot_2D2Osym_2D4O3sym
    end interface

    public :: operator(.tdot.)
    interface operator (.tdot.)
        module procedure tdot_2D2Osym_2D2Osym
    end interface

    contains

    pure module function ddot_2D4O3sym_2D2Osym(a, b) result(res)
        !
        !  | ( 1:1111) ( 5:1122) ( 8:1133) (10:1112) |
        !  | ( 5:2211) ( 2:2222) ( 6:2233) ( 9:2212) |
        !  | ( 8:3311) ( 6:3322) ( 3:3333) ( 7:3312) |
        !  | (10:1211) ( 9:1222) ( 7:1233) ( 4:1212) |
        implicit none
        class(ten_2D4O3sym), intent(in) :: a
        class(ten_2D2Osym), intent(in) :: b
        type(ten_2D2Osym) :: res
        
        real(real64) :: a1111b11, a1122b22, a1133b33, a1112b12
        real(real64) :: a2211b11, a2222b22, a2233b33, a2212b12
        real(real64) :: a3311b11, a3322b22, a3333b33, a3312b12
        real(real64) :: a1112b11, a2212b22, a3312b33, a1212b12

        a1111b11 = a%vals( 1)*b%vals(1)
        a1122b22 = a%vals( 5)*b%vals(2)
        a1133b33 = a%vals( 8)*b%vals(3)
        a1112b12 = a%vals(10)*b%vals(4)

        a2211b11 = a%vals( 5)*b%vals(1)
        a2222b22 = a%vals( 2)*b%vals(2)
        a2233b33 = a%vals( 6)*b%vals(3)
        a2212b12 = a%vals( 9)*b%vals(4)

        a3311b11 = a%vals( 8)*b%vals(1)
        a3322b22 = a%vals( 6)*b%vals(2)
        a3333b33 = a%vals( 3)*b%vals(3)
        a3312b12 = a%vals( 7)*b%vals(4)

        a1112b11 = a%vals(10)*b%vals(1)
        a2212b22 = a%vals( 9)*b%vals(2)
        a3312b33 = a%vals( 7)*b%vals(3)
        a1212b12 = a%vals( 4)*b%vals(4)
        
        res%vals(1) = a1111b11 + a1122b22 + a1133b33 + 2*a1112b12
        res%vals(2) = a2211b11 + a2222b22 + a2233b33 + 2*a2212b12
        res%vals(3) = a3311b11 + a3322b22 + a3333b33 + 2*a3312b12
        res%vals(4) = a1112b11 + a2212b22 + a3312b33 + 2*a1212b12
        
    end function ddot_2D4O3sym_2D2Osym

    pure module function ddot_2D2Osym_2D4O3sym(b, a) result(res)
        !
        !  | ( 1:1111) ( 5:1122) ( 8:1133) (10:1112) |
        !  | ( 5:2211) ( 2:2222) ( 6:2233) ( 9:2212) |
        !  | ( 8:3311) ( 6:3322) ( 3:3333) ( 7:3312) |
        !  | (10:1211) ( 9:1222) ( 7:1233) ( 4:1212) |
        implicit none
        class(ten_2D2Osym), intent(in) :: b
        class(ten_2D4O3sym), intent(in) :: a
        type(ten_2D2Osym) :: res
        
        real(real64) :: a1111b11, a2211b22, a3311b33, a1211b12
        real(real64) :: a1122b11, a2222b22, a3322b33, a1222b12
        real(real64) :: a1133b11, a2233b22, a3333b33, a1233b12
        real(real64) :: a1112b11, a2212b22, a3312b33, a1212b12

        a1111b11 = a%vals( 1)*b%vals(1)
        a2211b22 = a%vals( 5)*b%vals(2)
        a3311b33 = a%vals( 8)*b%vals(3)
        a1211b12 = a%vals(10)*b%vals(4)

        a1122b11 = a%vals( 5)*b%vals(1)
        a2222b22 = a%vals( 2)*b%vals(2)
        a3322b33 = a%vals( 6)*b%vals(3)
        a1222b12 = a%vals( 9)*b%vals(4)

        a1133b11 = a%vals( 8)*b%vals(1)
        a2233b22 = a%vals( 6)*b%vals(2)
        a3333b33 = a%vals( 3)*b%vals(3)
        a1233b12 = a%vals( 7)*b%vals(4)

        a1112b11 = a%vals(10)*b%vals(1)
        a2212b22 = a%vals( 9)*b%vals(2)
        a3312b33 = a%vals( 7)*b%vals(3)
        a1212b12 = a%vals( 4)*b%vals(4)
        
        res%vals(1) = a1111b11 + a2211b22 + a3311b33 + 2*a1211b12
        res%vals(2) = a1122b11 + a2222b22 + a3322b33 + 2*a1222b12
        res%vals(3) = a1133b11 + a2233b22 + a3333b33 + 2*a1233b12
        res%vals(4) = a1112b11 + a2212b22 + a3312b33 + 2*a1212b12
   
    end function ddot_2D2Osym_2D4O3sym

    pure module function tdot_2D2Osym_2D2Osym(a, b) result(res)
        !
        !  | ( 1:1111) ( 5:1122) ( 8:1133) (10:1112) |
        !  | ( 5:2211) ( 2:2222) ( 6:2233) ( 9:2212) |
        !  | ( 8:3311) ( 6:3322) ( 3:3333) ( 7:3312) |
        !  | (10:1211) ( 9:1222) ( 7:1233) ( 4:1212) |
        implicit none
        class(ten_2D2Osym), intent(in) :: a, b
        type(ten_2D4O3sym) :: res
        res%vals( 1) = a%vals(1)*b%vals(1)
        res%vals( 5) = a%vals(1)*b%vals(2)
        res%vals( 8) = a%vals(1)*b%vals(3)
        res%vals(10) = a%vals(1)*b%vals(4)

        res%vals( 2) = a%vals(2)*b%vals(2)
        res%vals( 6) = a%vals(2)*b%vals(3)
        res%vals( 9) = a%vals(2)*b%vals(4)

        res%vals( 3) = a%vals(3)*b%vals(3)
        res%vals( 7) = a%vals(3)*b%vals(4)

        res%vals( 4) = a%vals(4)*b%vals(4)
    end function tdot_2D2Osym_2D2Osym

end module mod_operator_2D2Osym_2D4O3sym