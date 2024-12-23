module mod_ten_3D4O3sym
    use, intrinsic :: iso_fortran_env
    implicit none
    private
 
    type, public :: ten_3D4O3sym
        real(real64), dimension(21) :: vals
        contains
            generic, public :: init => init_ten_3D4O3sym, init2_ten_3D4O3sym
            procedure, private :: init_ten_3D4O3sym, init2_ten_3D4O3sym 
    end type ten_3D4O3sym

    public :: operator(.isequal.)
    interface operator (.isequal.)
        module procedure isequal_3D4O3sym
    end interface

    public :: operator(.inv.)
    interface operator (.inv.)
        module procedure inv_3D4O3sym
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O3sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O3sym
        module procedure subU_3D4O3sym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_3D4O3sym
        module procedure mul_3D4O3sym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_3D4O3sym_real64
    end interface
contains
    
    pure module subroutine init_ten_3D4O3sym(self, vals)
        implicit none
        class(ten_3D4O3sym), intent(inout) :: self
        real(real64), intent(in) :: vals(21)
        self%vals = vals
    end subroutine init_ten_3D4O3sym

    pure module subroutine init2_ten_3D4O3sym(self,             & 
                                         xxxx, yyyy, zzzz, &
                                         xyxy, yzyz, xzxz, &
                                         xxyy, yyzz,       &
                                         zzxy, xyyz, yzxz, &
                                         xxzz,             &
                                         yyxy, zzyz, xyxz, &
                                         xxxy, yyyz, zzxz, &
                                         xxyz, yyxz, xxxz  &
                                         )
        !
        !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
        !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
        !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
        !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
        !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
        !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |
        implicit none
        class(ten_3D4O3sym), intent(inout) :: self
        real(real64), intent(in) :: xxxx, xxyy, xxzz, xxxy, xxyz, xxxz
        real(real64), intent(in) :: yyyy, yyzz, yyxy, yyyz, yyxz
        real(real64), intent(in) :: zzzz, zzxy, zzyz, zzxz
        real(real64), intent(in) :: xyxy, xyyz, xyxz
        real(real64), intent(in) :: yzyz, yzxz
        real(real64), intent(in) :: xzxz
        self%vals = (/xxxx, yyyy, zzzz, xyxy, yzyz, xzxz, &
                      xxyy, yyzz, zzxy, xyyz, yzxz, &
                      xxzz, yyxy, zzyz, xyxz, &
                      xxxy, yyyz, zzxz, &
                      xxyz, yyxz, &
                      xxxz &
                      /)
    end subroutine

    pure module function isequal_3D4O3sym(a, b) result(res)
        implicit none
        class(ten_3D4O3sym), intent(in) :: a, b
        logical :: res
        real(real64), parameter :: EPS=1e-7, EPS_ABS=1e-30
        real(real64) :: norm_a, norm_b, norm_max, norm

        norm_a =   abs(a%vals( 1)) +   abs(a%vals( 2)) +   abs(a%vals( 3)) + &
                   abs(a%vals( 4)) +   abs(a%vals( 5)) +   abs(a%vals( 6)) + &
                 2*abs(a%vals( 7)) + 2*abs(a%vals( 8)) + 2*abs(a%vals( 9)) + &
                 2*abs(a%vals(10)) + 2*abs(a%vals(11)) + 2*abs(a%vals(12)) + & 
                 2*abs(a%vals(13)) + 2*abs(a%vals(14)) + 2*abs(a%vals(15)) + &
                 2*abs(a%vals(16)) + 2*abs(a%vals(17)) + 2*abs(a%vals(18)) + &
                 2*abs(a%vals(19)) + 2*abs(a%vals(20)) + 2*abs(a%vals(21))
        norm_b =   abs(b%vals( 1)) +   abs(b%vals( 2)) +   abs(b%vals( 3)) + &
                   abs(b%vals( 4)) +   abs(b%vals( 5)) +   abs(b%vals( 6)) + &
                 2*abs(b%vals( 7)) + 2*abs(b%vals( 8)) + 2*abs(b%vals( 9)) + &
                 2*abs(b%vals(10)) + 2*abs(b%vals(11)) + 2*abs(b%vals(12)) + & 
                 2*abs(b%vals(13)) + 2*abs(b%vals(14)) + 2*abs(b%vals(15)) + &
                 2*abs(b%vals(16)) + 2*abs(b%vals(17)) + 2*abs(b%vals(18)) + &
                 2*abs(b%vals(19)) + 2*abs(b%vals(20)) + 2*abs(b%vals(21))  

        norm_max = max(max(norm_a, norm_b), EPS_ABS)

        norm =   abs(a%vals( 1)-b%vals( 1)) +   abs(a%vals( 2)-b%vals( 2)) +   abs(a%vals( 3)-b%vals( 3)) + &
                 abs(a%vals( 4)-b%vals( 4)) +   abs(a%vals( 5)-b%vals( 5)) +   abs(a%vals( 6)-b%vals( 6)) + &
               2*abs(a%vals( 7)-b%vals( 7)) + 2*abs(a%vals( 8)-b%vals( 8)) + 2*abs(a%vals( 9)-b%vals( 9)) + &
               2*abs(a%vals(10)-b%vals(10)) + 2*abs(a%vals(11)-b%vals(11)) + 2*abs(a%vals(12)-b%vals(12)) + &
               2*abs(a%vals(13)-b%vals(13)) + 2*abs(a%vals(14)-b%vals(14)) + 2*abs(a%vals(15)-b%vals(15)) + &
               2*abs(a%vals(16)-b%vals(16)) + 2*abs(a%vals(17)-b%vals(17)) + 2*abs(a%vals(18)-b%vals(18)) + &
               2*abs(a%vals(19)-b%vals(19)) + 2*abs(a%vals(20)-b%vals(20)) + 2*abs(a%vals(21)-b%vals(21))
        
        if (norm/norm_max .gt. EPS) res=.false.
        if (norm/norm_max .le. EPS) res=.true.
    end function isequal_3D4O3sym

    pure module function sum_3D4O3sym(a, b) result(res)
        implicit none
        class(ten_3D4O3sym), intent(in) :: a, b
        type(ten_3D4O3sym) :: res
        res%vals = a%vals + b%vals
    end function sum_3D4O3sym

    pure module function sub_3D4O3sym(a, b) result(res)
        implicit none
        class(ten_3D4O3sym), intent(in) :: a, b
        type(ten_3D4O3sym) :: res
        res%vals = a%vals - b%vals
    end function sub_3D4O3sym

    pure module function subU_3D4O3sym(a) result(res)
        implicit none
        class(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        res%vals = -a%vals
    end function subU_3D4O3sym

    pure module function mul_real64_3D4O3sym(a, b) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O3sym), intent(in) :: b
        type(ten_3D4O3sym) :: res
        res%vals = a * b%vals
    end function mul_real64_3D4O3sym

    pure module function mul_3D4O3sym_real64(b, a) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O3sym), intent(in) :: b
        type(ten_3D4O3sym) :: res
        res%vals = a * b%vals
    end function mul_3D4O3sym_real64

    pure module function div_3D4O3sym_real64(b, a) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O3sym), intent(in) :: b
        type(ten_3D4O3sym) :: res
        res%vals = b%vals/a
    end function div_3D4O3sym_real64

    pure module function inv_3D4O3sym(a) result(res)
        use, intrinsic :: iso_fortran_env
        use inverses_mat
        implicit none
        class(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        real(real64) :: mat_a(6,6), mat_b(6,6), v(21)
        logical :: ok
        integer :: iok
        v = a%vals
        mat_a = reshape((/  v(1),  v(7), v(12), v(16), v(19), v(21), &
                            v(7),  v(2),  v(8), v(13), v(17), v(20), &
                           v(12),  v(8),  v(3),  v(9), v(14), v(18), & 
                           v(16), v(13),  v(9),  v(4), v(10), v(15), &
                           v(19), v(17), v(14), v(10),  v(5), v(11), &
                           v(21), v(20), v(18), v(15), v(11),  v(6)  &
                        /), (/6,6/))

        call M66INV(mat_a, mat_b, ok)


        ! call FINDInv(mat_a, mat_b, 6, iok)

        ! v = (/ mat_b(1,1), mat_b(2,2), mat_b(3,3), mat_b(4,4), mat_b(5,5), mat_b(6,6),  &
        !        mat_b(1,2), mat_b(2,3), mat_b(3,4), mat_b(4,5), mat_b(5,6),              &
        !        mat_b(1,3), mat_b(2,4), mat_b(3,5), mat_b(4,6),                          &
        !        mat_b(1,4), mat_b(2,5), mat_b(3,6),                                      &
        !        mat_b(1,5), mat_b(2,6),                                                  &
        !        mat_b(1,6)                                                               &
        !       /)

        v = (/ mat_b(1,1), mat_b(2,2), mat_b(3,3), mat_b(4,4)/4D0, mat_b(5,5)/4D0, mat_b(6,6)/4D0,  &
               mat_b(1,2), mat_b(2,3), mat_b(3,4)/2D0, mat_b(4,5)/4D0, mat_b(5,6)/4D0,              &
               mat_b(1,3), mat_b(2,4)/2D0, mat_b(3,5)/2D0, mat_b(4,6)/4D0,                          &
               mat_b(1,4)/2D0, mat_b(2,5)/2D0, mat_b(3,6)/2D0,                                      &
               mat_b(1,5)/2D0, mat_b(2,6)/2D0,                                                  &
               mat_b(1,6)/2D0                                                               &
              /)


        ! mat_a = reshape((/  v(1),  v(7), v(12), v(16), v(19), v(21), &
        !                     v(7),  v(2),  v(8), v(13), v(17), v(20), &
        !                     v(12),  v(8),  v(3),  v(9), v(14), v(18), & 
        !                     v(16), v(13),  v(9),  v(4), v(10), v(15), &
        !                     v(19), v(17), v(14), v(10),  v(5), v(11), &
        !                     v(21), v(20), v(18), v(15), v(11),  v(6)  &
        !                     /), (/6,6/))
        ! call M66INV(mat_a, mat_b, ok)
        ! call FINDInv(mat_a, mat_b, 6, iok)
        
        call res%init(v)
        
            !
            !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
            !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
            !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
            !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
            !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
            !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |

    end function inv_3D4O3sym

end module mod_ten_3D4O3sym