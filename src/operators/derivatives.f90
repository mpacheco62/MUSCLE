module derivatives
    private
    public :: derivative
    interface derivative
        module procedure derivative_scalar_scalar
        module procedure derivative_scalar_3D2O
        module procedure derivative_scalar_3D2Osym
        module procedure derivative_3D2Osym_3D2Osym
    end interface

    public :: derivative2O
    interface derivative2O
        module procedure derivative2O_scalar_3D2Osym
    end interface

    interface
        pure function f_scalar_scalar(x)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64), intent(in) :: x
            real(real64) :: f_scalar_scalar
        end function f_scalar_scalar
        
        pure function f_scalar_3D2O(x)
            use tensors_types, only : ten_3D2O
            use, intrinsic :: iso_fortran_env
            implicit none
            type(ten_3D2O), intent(in) :: x
            real(real64) :: f_scalar_3D2O
        end function f_scalar_3D2O

        pure function f_scalar_3D2Osym(x)
            use tensors_types, only : ten_3D2Osym
            use, intrinsic :: iso_fortran_env
            implicit none
            type(ten_3D2Osym), intent(in) :: x
            real(real64) :: f_scalar_3D2Osym
        end function f_scalar_3D2Osym

        pure function f_3D2Osym_3D2Osym(x)
            use tensors_types, only : ten_3D2Osym
            use, intrinsic :: iso_fortran_env
            implicit none
            type(ten_3D2Osym), intent(in) :: x
            type(ten_3D2Osym) :: f_3D2Osym_3D2Osym
        end function f_3D2Osym_3D2Osym
    end interface

    contains

    pure function derivative_scalar_scalar(func, x, eps)
        use, intrinsic :: iso_fortran_env
        implicit none

        procedure(f_scalar_scalar) :: func
        real(real64), intent(in) :: x
        real(real64), intent(in), optional :: eps
        real(real64) :: derivative_scalar_scalar
        
        real(real64), parameter :: DIVEPS = 1D-7, MAX_EPS=1D-40  !! Relative and minimum absolute step size
        real(real64) :: x_var_forw, x_var_back
        real(real64) :: func_forw, func_back, func_val
        real(real64) :: epsr

        if (present(eps)) then
            epsr = eps
        else
            func_val = func(x)
            epsr = max(abs(x*DIVEPS), MAX_EPS)
        end if

        derivative_scalar_scalar = 0.0D0
        x_var_forw = x + epsr
        func_forw = func(x_var_forw)
        x_var_back = x - epsr
        func_back = func(x_var_back)
        derivative_scalar_scalar = (func_forw - func_back)/(2D0*epsr)
    end function derivative_scalar_scalar

    pure function derivative_scalar_3D2O(func, mat, eps)
        use tensors_types, only : ten_3D2O
        use, intrinsic :: iso_fortran_env
        implicit none

        procedure(f_scalar_3D2O) :: func
        type(ten_3D2O), intent(in) :: mat
        real(real64), intent(in), optional :: eps
        type(ten_3D2O) :: derivative_scalar_3D2O
        
        real(real64), parameter :: DIVEPS = 1D-7, MAX_EPS=1D-40  !! Relative and minimum absolute step size
        type(ten_3D2O) :: mat_var_forw, mat_var_back
        real(real64) :: val_func_forw, val_func_back, mat_norm
        real(real64) :: epsr
        integer :: i

        if (present(eps)) then
            epsr = eps
        else
            mat_norm = mat%norm()
            epsr = max(abs(mat_norm*DIVEPS), MAX_EPS)
        end if

        derivative_scalar_3D2O%vals = 0.0D0
        do i=1,9
            mat_var_forw = mat
            mat_var_forw%vals(i) = mat_var_forw%vals(i) + epsr
            val_func_forw = func(mat_var_forw)
            mat_var_back = mat
            mat_var_back%vals(i) = mat_var_back%vals(i) - epsr
            val_func_back = func(mat_var_back)
            derivative_scalar_3D2O%vals(i) = (val_func_forw - val_func_back)/(2D0*epsr)
        end do
    end function derivative_scalar_3D2O

    pure function derivative_scalar_3D2Osym(func, mat, eps)
        use tensors_types, only : ten_3D2Osym
        use, intrinsic :: iso_fortran_env
        implicit none

        procedure(f_scalar_3D2Osym) :: func
        type(ten_3D2Osym), intent(in) :: mat
        real(real64), intent(in), optional :: eps
        type(ten_3D2Osym) :: derivative_scalar_3D2OSym
        
        real(real64), parameter :: DIVEPS = 1D-6, MAX_EPS=1D-40  !! Relative and minimum absolute step size
        type(ten_3D2Osym) :: mat_var
        real(real64) :: val_func_forw, val_func_back, mat_norm
        real(real64) :: epsr

        if (present(eps)) then
            epsr = eps
        else
            mat_norm = mat%norm()
            epsr = max(abs(mat_norm*DIVEPS), MAX_EPS)
        end if

        derivative_scalar_3D2OSym%vals = 0.0D0

        mat_var = mat
        mat_var%vals(1) = mat%vals(1) + epsr; val_func_forw = func(mat_var)
        mat_var%vals(1) = mat%vals(1) - epsr; val_func_back = func(mat_var)
        derivative_scalar_3D2Osym%vals(1) = (val_func_forw - val_func_back)/(2D0*epsr)
        mat_var%vals(1) = mat%vals(1)

        mat_var%vals(2) = mat%vals(2) + epsr; val_func_forw = func(mat_var)
        mat_var%vals(2) = mat%vals(2) - epsr; val_func_back = func(mat_var)
        derivative_scalar_3D2OSym%vals(2) = (val_func_forw - val_func_back)/(2D0*epsr)
        mat_var%vals(2) = mat%vals(2)

        mat_var%vals(3) = mat%vals(3) + epsr; val_func_forw = func(mat_var)
        mat_var%vals(3) = mat%vals(3) - epsr; val_func_back = func(mat_var)
        derivative_scalar_3D2OSym%vals(3) = (val_func_forw - val_func_back)/(2D0*epsr)
        mat_var%vals(3) = mat%vals(3)

        mat_var%vals(4) = mat%vals(4) + epsr/2D0; val_func_forw = func(mat_var)
        mat_var%vals(4) = mat%vals(4) - epsr/2D0; val_func_back = func(mat_var)
        derivative_scalar_3D2OSym%vals(4) = (val_func_forw - val_func_back)/(2D0*epsr)
        mat_var%vals(4) = mat%vals(4)

        mat_var%vals(5) = mat%vals(5) + epsr/2D0; val_func_forw = func(mat_var)
        mat_var%vals(5) = mat%vals(5) - epsr/2D0; val_func_back = func(mat_var)
        derivative_scalar_3D2OSym%vals(5) = (val_func_forw - val_func_back)/(2D0*epsr)
        mat_var%vals(5) = mat%vals(5)

        mat_var%vals(6) = mat%vals(6) + epsr/2D0; val_func_forw = func(mat_var)
        mat_var%vals(6) = mat%vals(6) - epsr/2D0; val_func_back = func(mat_var)
        derivative_scalar_3D2OSym%vals(6) = (val_func_forw - val_func_back)/(2D0*epsr)
        mat_var%vals(6) = mat%vals(6)

    end function derivative_scalar_3D2OSym


    pure function derivative2O_scalar_3D2Osym(func, mat, eps) result(res)
        use tensors_types, only : ten_3D2Osym, ten_3D4O3sym,ten_3D4O2sym
        use, intrinsic :: iso_fortran_env
        implicit none

        procedure(f_scalar_3D2Osym) :: func
        type(ten_3D2Osym), intent(in) :: mat
        real(real64), intent(in), optional :: eps
        type(ten_3D4O2sym) :: a
        type(ten_3D4O3sym) :: res
        type(ten_3D2Osym) :: mat_var
        integer :: i,j
        real(real64), parameter :: DIVEPS = 1D-4, MAX_EPS=1D-40  !! Relative and minimum absolute step size
        real(real64) :: mat_norm, val_func_forw, val_func_back
        real(real64) :: epsr,f_pp, f_pm, f_mp, f_mm

        if (present(eps)) then
            epsr = eps
        else
            mat_norm = mat%norm()
            epsr = max(abs(mat_norm*DIVEPS), MAX_EPS)
        end if
        !  Voigt Matrix Layout (I, J):
        !  | (1,1) (1,2) (1,3) (1,4) (1,5) (1,6) |   <- xxxx, xxyy, xxzz, xxxy, xxyz, xxxz
        !  | (2,1) (2,2) (2,3) (2,4) (2,5) (2,6) |   <- yyxx, yyyy, yyzz, yyxy, yyyz, yyxz
        !  | (3,1) (3,2) (3,3) (3,4) (3,5) (3,6) |   <- zzxx, zzyy, zzzz, zzxy, zzyz, zzxz
        !  | (4,1) (4,2) (4,3) (4,4) (4,5) (4,6) |   <- xyxx, xyyy, xyzz, xyxy, xyyz, xyxz
        !  | (5,1) (5,2) (5,3) (5,4) (5,5) (5,6) |   <- yzxx, yzyy, yzzz, yzxy, yzyz, yzxz
        !  | (6,1) (6,2) (6,3) (6,4) (6,5) (6,6) |   <- xzxx, xzyy, xzzz, xzxy, xzyz, xzxz
        !--------------------------------------------------------------------------------------
        !! Voigt Matrix (Indices IJ):
        !! ```
        !! | 11 12 13 14 15 16 |
        !! |    22 23 24 25 26 |
        !! |       33 34 35 36 |
        !! |          44 45 46 |
        !! |             55 56 |
        !! |                66 |
        !! ```
        !! Storage order in `vals(1:21)`:
        !! (11, 22, 33, 44, 55, 66, 12, 23, 34, 45, 56, 13, 24, 35, 46, 14, 25, 36, 15, 26, 16)
        a%vals = 0.0D0
        res%vals = 0.0D0
        mat_var = mat
        !Para las variables xxxx,yyyy,zzzz (1,1),(2,2),(3,3) con el if
        !Para las otras variables, como xxyy, (1,2) (1,3),(2,3),(2,1),(3,1) (3,2) con el else
        do i=1,3
            do j=1,3
                if (i == j) then 
                     mat_var%vals(i) = mat%vals(i) + 2*epsr; val_func_forw = func(mat_var)
                     mat_var%vals(i) = mat%vals(i) - 2*epsr; val_func_back = func(mat_var)
                     a%vals(i,i) = (val_func_forw - 2D0*func(mat) + val_func_back)/(4*epsr*epsr)
                     mat_var%vals(i) = mat%vals(i)
                else
                    ! (+epsr, +epsr)
                    mat_var%vals(i) = mat%vals(i) + epsr
                    mat_var%vals(j) = mat%vals(j) + epsr   
                    f_pp = func(mat_var)
                    ! (+epsr, -epsr)
                    mat_var%vals(i) = mat%vals(i) + epsr
                    mat_var%vals(j) = mat%vals(j) - epsr
                    f_pm = func(mat_var)
                    ! (-epsr, +epsr)
                    mat_var%vals(i) = mat%vals(i) - epsr
                    mat_var%vals(j) = mat%vals(j) + epsr
                    f_mp = func(mat_var)
                    ! (-epsr, -epsr)
                    mat_var%vals(i) = mat%vals(i) - epsr
                    mat_var%vals(j) = mat%vals(j) - epsr
                    f_mm = func(mat_var)
                    ! derivada
                    a%vals(i,j) = (f_pp - f_pm - f_mp + f_mm) / (4D0*epsr*epsr)
                    mat_var%vals(i) = mat%vals(i)
                    mat_var%vals(j) = mat%vals(j)
                end if
            end do
        end do
        !------------------------------------------------------------------------------------
        !Para las otras variables con e/2 al inicio y e al final, como xyxx
        !(4,1),(4,2),(4,3),(5,1),(5,2),(5,3),(6,1),(6,2),(6,3)
        do i=4,6
            do j=1,3
                ! (+epsr, +epsr)
                mat_var%vals(i) = mat%vals(i) + epsr/2D0
                mat_var%vals(j) = mat%vals(j) + epsr  
                f_pp = func(mat_var)
                ! (+epsr, -epsr)
                mat_var%vals(i) = mat%vals(i) + epsr/2D0
                mat_var%vals(j) = mat%vals(j) - epsr
                f_pm = func(mat_var)
                ! (-epsr, +epsr)
                mat_var%vals(i) = mat%vals(i) - epsr/2D0
                mat_var%vals(j) = mat%vals(j) + epsr
                f_mp = func(mat_var)
                ! (-epsr, -epsr)
                mat_var%vals(i) = mat%vals(i) - epsr/2D0
                mat_var%vals(j) = mat%vals(j) - epsr
                f_mm = func(mat_var)
                ! derivada
                a%vals(i,j) = (f_pp - f_pm - f_mp + f_mm) / (4D0*epsr*epsr)
                mat_var%vals(i) = mat%vals(i)
                mat_var%vals(j) = mat%vals(j)
            end do
        end do
        !------------------------------------------------------------------------
        !Para las variables xyxy,yzyz,xzxz (4,4),(5,5),(6,6) con el if
        !Para variables con dos e/2 pero distintas ej yzxy, con el else
        !(4,5),(4,6),(5,6),(5,4),(6,4),(6,5)
        do i=4,6
            do j=4,6
                if (i == j) then
                    mat_var%vals(i) = mat%vals(i) + epsr; val_func_forw = func(mat_var)
                    mat_var%vals(i) = mat%vals(i) - epsr; val_func_back = func(mat_var)
                    a%vals(i,i) = (val_func_forw - 2D0*func(mat) + val_func_back)/(4*epsr*epsr)
                    mat_var%vals(i) = mat%vals(i)
                else
                    ! (+epsr, +epsr)
                    mat_var%vals(i) = mat%vals(i) + epsr/2D0
                    mat_var%vals(j) = mat%vals(j) + epsr/2D0   
                    f_pp = func(mat_var)
                    ! (+epsr, -epsr)
                    mat_var%vals(i) = mat%vals(i) + epsr/2D0
                    mat_var%vals(j) = mat%vals(j) - epsr/2D0
                    f_pm = func(mat_var)
                    ! (-epsr, +epsr)
                    mat_var%vals(i) = mat%vals(i) - epsr/2D0
                    mat_var%vals(j) = mat%vals(j) + epsr/2D0
                    f_mp = func(mat_var)
                    ! (-epsr, -epsr)
                    mat_var%vals(i) = mat%vals(i) - epsr/2D0
                    mat_var%vals(j) = mat%vals(j) - epsr/2D0
                    f_mm = func(mat_var)
                    ! derivada
                    a%vals(i,j) = (f_pp - f_pm - f_mp + f_mm) / (4D0*epsr*epsr)
                    mat_var%vals(i) = mat%vals(i)
                    mat_var%vals(j) = mat%vals(j)
                end if
            end do
        end do
        !! (11, 22, 33, 44, 55, 66, 12, 23, 34, 45, 56, 13, 24, 35, 46, 14, 25, 36, 15, 26, 16)
        res%vals(1)=a%vals(1,1)  ; res%vals(2)=a%vals(2,2)  ; res%vals(3)=a%vals(3,3)  ; res%vals(4)=a%vals(4,4)  
        res%vals(5)=a%vals(5,5)  ; res%vals(6)=a%vals(6,6)  ; res%vals(7)=a%vals(1,2)  ; res%vals(8)=a%vals(2,3)  
        res%vals(9)=a%vals(3,4)  ; res%vals(10)=a%vals(4,5) ; res%vals(11)=a%vals(5,6) ; res%vals(12)=a%vals(1,3)
        res%vals(13)=a%vals(2,4) ; res%vals(14)=a%vals(3,5) ; res%vals(15)=a%vals(4,6) ; res%vals(16)=a%vals(1,4) 
        res%vals(17)=a%vals(2,5) ; res%vals(18)=a%vals(3,6) ; res%vals(19)=a%vals(1,5) ; res%vals(20)=a%vals(2,6) 
        res%vals(21)=a%vals(1,6)
        return
    end function derivative2O_scalar_3D2Osym
    

    pure function derivative_3D2Osym_3D2Osym(func, mat, eps)
        use tensors_types, only : ten_3D2Osym, ten_3D4O2sym, &
                                  operator(+), operator(-), operator(*), operator(/), assignment(=)
        use, intrinsic :: iso_fortran_env
        implicit none

        procedure(f_3D2Osym_3D2Osym) :: func
        type(ten_3D2Osym), intent(in) :: mat
        real(real64), intent(in), optional :: eps
        type(ten_3D4O2sym) :: derivative_3D2Osym_3D2Osym
        
        type(ten_3D2Osym) :: mat_var, tmp
        type(ten_3D2Osym) :: val_func_forw, val_func_back
        real(real64) :: epsr, mat_norm

        real(real64), parameter :: DIVEPS = 1D-4, MAX_EPS=1D-35  !! Relative and minimum absolute step size
        

        if (present(eps)) then
            epsr = eps
        else
            mat_norm = mat%norm()
            epsr = max(abs(mat_norm*DIVEPS), MAX_EPS)
        end if

        derivative_3D2Osym_3D2Osym%vals = 0.0D0

        mat_var = mat
        mat_var%vals(1) = mat%vals(1) + epsr; val_func_forw = func(mat_var)
        mat_var%vals(1) = mat%vals(1) - epsr; val_func_back = func(mat_var)
        tmp = (val_func_forw - val_func_back)/(2D0*epsr)
        derivative_3D2Osym_3D2Osym%vals(:,1) = tmp%vals
        mat_var%vals(1) = mat%vals(1)

        mat_var%vals(2) = mat%vals(2) + epsr; val_func_forw = func(mat_var)
        mat_var%vals(2) = mat%vals(2) - epsr; val_func_back = func(mat_var)
        tmp = (val_func_forw - val_func_back)/(2D0*epsr)
        derivative_3D2Osym_3D2Osym%vals(:,2) = tmp%vals
        mat_var%vals(2) = mat%vals(2)

        mat_var%vals(3) = mat%vals(3) + epsr; val_func_forw = func(mat_var)
        mat_var%vals(3) = mat%vals(3) - epsr; val_func_back = func(mat_var)
        tmp = (val_func_forw - val_func_back)/(2D0*epsr)
        derivative_3D2Osym_3D2Osym%vals(:,3) = tmp%vals
        mat_var%vals(3) = mat%vals(3)

        mat_var%vals(4) = mat%vals(4) + epsr/2D0; val_func_forw = func(mat_var)
        mat_var%vals(4) = mat%vals(4) - epsr/2D0; val_func_back = func(mat_var)
        tmp = (val_func_forw - val_func_back)/(2D0*epsr)
        derivative_3D2Osym_3D2Osym%vals(:,4) = tmp%vals
        mat_var%vals(4) = mat%vals(4)

        mat_var%vals(5) = mat%vals(5) + epsr/2D0; val_func_forw = func(mat_var)
        mat_var%vals(5) = mat%vals(5) - epsr/2D0; val_func_back = func(mat_var)
        tmp = (val_func_forw - val_func_back)/(2D0*epsr)
        derivative_3D2Osym_3D2Osym%vals(:,5) = tmp%vals
        mat_var%vals(5) = mat%vals(5)

        mat_var%vals(6) = mat%vals(6) + epsr/2D0; val_func_forw = func(mat_var)
        mat_var%vals(6) = mat%vals(6) - epsr/2D0; val_func_back = func(mat_var)
        tmp = (val_func_forw - val_func_back)/(2D0*epsr)
        derivative_3D2Osym_3D2Osym%vals(:,6) = tmp%vals
        mat_var%vals(6) = mat%vals(6)

    end function derivative_3D2Osym_3D2Osym

end module