module mod_operator_3D2Osym_3D4O2sym
    use, intrinsic :: iso_fortran_env
    use mod_ten_3D2Osym
    use mod_ten_3D4O2sym
    implicit none
    private


    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D4O2sym_3D2Osym
        module procedure ddot_3D2Osym_3D4O2sym
    end interface

    contains

    pure function ddot_3D4O2sym_3D2Osym(a, b) result(res)
        !
        !  | (1111) (1122) (1133) (1112) (1123) (1113) |
        !  | (2211) (2222) (2233) (2212) (2223) (2213) |
        !  | (3311) (3322) (3333) (3312) (3323) (3313) |
        !  | (1211) (1222) (1233) (1212) (1223) (1213) |
        !  | (2311) (2322) (2333) (2312) (2323) (2313) |
        !  | (1311) (1322) (1333) (1312) (1323) (1313) |
        implicit none
        class(ten_3D4O2sym), intent(in) :: a
        class(ten_3D2Osym), intent(in) :: b
        type(ten_3D2Osym) :: res
        
        
        res%vals(:) =   a%vals(:,1)*b%vals(1) + a%vals(:,2)*b%vals(2) + a%vals(:,3)*b%vals(3) +    &
                      2*a%vals(:,4)*b%vals(4) + 2*a%vals(:,5)*b%vals(5) + 2*a%vals(:,6)*b%vals(6)
        
    end function ddot_3D4O2sym_3D2Osym

    pure function ddot_3D2Osym_3D4O2sym(b, a) result(res)
        !
        !  | (1111) (1122) (1133) (1112) (1123) (1113) |
        !  | (2211) (2222) (2233) (2212) (2223) (2213) |
        !  | (3311) (3322) (3333) (3312) (3323) (3313) |
        !  | (1211) (1222) (1233) (1212) (1223) (1213) |
        !  | (2311) (2322) (2333) (2312) (2323) (2313) |
        !  | (1311) (1322) (1333) (1312) (1323) (1313) |
        implicit none
        class(ten_3D2Osym), intent(in) :: b
        class(ten_3D4O2sym), intent(in) :: a
        type(ten_3D2Osym) :: res

        res%vals(:) =   a%vals(:,1)*b%vals(1) + a%vals(:,2)*b%vals(2) + a%vals(:,3)*b%vals(3) +    &
                      2*a%vals(:,4)*b%vals(4) + 2*a%vals(:,5)*b%vals(5) + 2*a%vals(:,6)*b%vals(6)
        
    end function ddot_3D2Osym_3D4O2sym

end module mod_operator_3D2Osym_3D4O2sym