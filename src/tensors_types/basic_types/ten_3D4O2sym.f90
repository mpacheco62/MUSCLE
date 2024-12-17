module mod_ten_3D4O2sym
    use, intrinsic :: iso_fortran_env
    implicit none
    private
 
    type, public :: ten_3D4O2sym
        real(real64), dimension(6,6) :: vals
        contains
            generic, public :: init => init_ten_3D4O2sym, init2_ten_3D4O2sym
            procedure, private :: init_ten_3D4O2sym, init2_ten_3D4O2sym 
    end type ten_3D4O2sym

    public :: operator(.isequal.)
    interface operator (.isequal.)
        module procedure isequal_3D4O2sym
    end interface

    ! public :: operator(.inv.)
    ! interface operator (.inv.)
    !     module procedure inv_3D4O2sym
    ! end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O2sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O2sym
        module procedure subU_3D4O2sym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_3D4O2sym
        module procedure mul_3D4O2sym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_3D4O2sym_real64
    end interface
contains
    
    pure module subroutine init_ten_3D4O2sym(self, vals)
        implicit none
        class(ten_3D4O2sym), intent(inout) :: self
        real(real64), intent(in) :: vals(6,6)
        self%vals = vals
    end subroutine init_ten_3D4O2sym

    pure module subroutine init2_ten_3D4O2sym(self,                               & 
                                              xxxx, xxyy, xxzz, xxxy, xxyz, xxxz, &
                                              yyxx, yyyy, yyzz, yyxy, yyyz, yyxz, &
                                              zzxx, zzyy, zzzz, zzxy, zzyz, zzxz, &
                                              xyxx, xyyy, xyzz, xyxy, xyyz, xyxz, &
                                              yzxx, yzyy, yzzz, yzxy, yzyz, yzxz, &
                                              xzxx, xzyy, xzzz, xzxy, xzyz, xzxz  &
                                              )
        !
        !  | (:1111) (1122) (1133) (1112) (1123) (1113) |
        !  | (:2211) (2222) (2233) (2212) (2223) (2213) |
        !  | (:3311) (3322) (3333) (3312) (3323) (3313) |
        !  | (:1211) (1222) (1233) (1212) (1223) (1213) |
        !  | (:2311) (2322) (2333) (2312) (2323) (2313) |
        !  | (:1311) (1322) (1333) (1312) (1323) (1313) |
        implicit none
        class(ten_3D4O2sym), intent(inout) :: self
        real(real64), intent(in) :: xxxx, xxyy, xxzz, xxxy, xxyz, xxxz
        real(real64), intent(in) :: yyxx, yyyy, yyzz, yyxy, yyyz, yyxz
        real(real64), intent(in) :: zzxx, zzyy, zzzz, zzxy, zzyz, zzxz
        real(real64), intent(in) :: xyxx, xyyy, xyzz, xyxy, xyyz, xyxz
        real(real64), intent(in) :: yzxx, yzyy, yzzz, yzxy, yzyz, yzxz
        real(real64), intent(in) :: xzxx, xzyy, xzzz, xzxy, xzyz, xzxz
        self%vals(:,1) = (/ xxxx, yyxx, zzxx, xyxx, yzxx, xzxx /)
        self%vals(:,2) = (/ xxyy, yyyy, zzyy, xyyy, yzyy, xzyy /)
        self%vals(:,3) = (/ xxzz, yyzz, zzzz, xyzz, yzzz, xzzz /)
        self%vals(:,4) = (/ xxxy, yyxy, zzxy, xyxy, yzxy, xzxy /)
        self%vals(:,5) = (/ xxyz, yyyz, zzyz, xyyz, yzyz, xzyz /)
        self%vals(:,6) = (/ xxxz, yyxz, zzxz, xyxz, yzxz, xzxz /)
    end subroutine

    pure module function isequal_3D4O2sym(a, b) result(res)
        implicit none
        class(ten_3D4O2sym), intent(in) :: a, b
        type(ten_3D4O2sym) :: temp
        logical :: res
        real(real64), parameter :: EPS=1e-7, EPS_ABS=1e-30
        real(real64) :: norm_a, norm_b, norm_max, norm
        integer :: i, j

        norm_a = 0D0
        norm_b = 0D0
        do i=1,6
            do j=1,6
                norm_a = norm_a + abs(a%vals(i,j))
                norm_b = norm_b + abs(b%vals(i,j))
            end do
        end do

        norm_max = max(max(norm_a, norm_b), EPS_ABS)

        temp = a - b
        norm = 0D0
        do i=1,6
            do j=1,6
                norm = norm + abs(temp%vals(i,j))
            end do
        end do
        
        if (norm/norm_max .gt. EPS) res=.false.
        if (norm/norm_max .le. EPS) res=.true.
    end function isequal_3D4O2sym

    pure module function sum_3D4O2sym(a, b) result(res)
        implicit none
        class(ten_3D4O2sym), intent(in) :: a, b
        type(ten_3D4O2sym) :: res
        res%vals = a%vals + b%vals
    end function sum_3D4O2sym

    pure module function sub_3D4O2sym(a, b) result(res)
        implicit none
        class(ten_3D4O2sym), intent(in) :: a, b
        type(ten_3D4O2sym) :: res
        res%vals = a%vals - b%vals
    end function sub_3D4O2sym

    pure module function subU_3D4O2sym(a) result(res)
        implicit none
        class(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = -a%vals
    end function subU_3D4O2sym

    pure module function mul_real64_3D4O2sym(a, b) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O2sym), intent(in) :: b
        type(ten_3D4O2sym) :: res
        res%vals = a * b%vals
    end function mul_real64_3D4O2sym

    pure module function mul_3D4O2sym_real64(b, a) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O2sym), intent(in) :: b
        type(ten_3D4O2sym) :: res
        res%vals = a * b%vals
    end function mul_3D4O2sym_real64

    pure module function div_3D4O2sym_real64(b, a) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O2sym), intent(in) :: b
        type(ten_3D4O2sym) :: res
        res%vals = b%vals/a
    end function div_3D4O2sym_real64

    pure module function inv_3D4O2sym(a) result(res)
        use, intrinsic :: iso_fortran_env
        use inverses_mat
        implicit none
        class(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        real(real64) :: mat_a(6,6), mat_b(6,6)
        logical :: ok
        mat_a = a%vals

        call M66INV(mat_a, mat_b, ok)

        mat_b(1:3, 4:6) = mat_b(1:3, 4:6)/2D0
        mat_b(4:6, 1:3) = mat_b(4:6, 1:3)/2D0
        mat_b(4:6, 4:6) = mat_b(4:6, 4:6)/4D0
        
        call res%init(mat_b)

    end function inv_3D4O2sym

end module mod_ten_3D4O2sym