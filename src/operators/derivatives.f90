module derivatives
    private
    public :: derivative
    interface derivative
        module procedure derivative_scalar_3D2O
        module procedure derivative_scalar_3D2Osym
    end interface

    interface
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
    end interface

    contains

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
        real(real64) :: val_func_forw, val_func_back, val_func
        real(real64) :: epsr
        integer :: i

        if (present(eps)) then
            epsr = eps
        else
            val_func = func(mat)
            epsr = max(abs(val_func*DIVEPS), MAX_EPS)
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
        
        real(real64), parameter :: DIVEPS = 1D-7, MAX_EPS=1D-40  !! Relative and minimum absolute step size
        type(ten_3D2Osym) :: mat_var
        real(real64) :: val_func_forw, val_func_back, val_func
        real(real64) :: epsr

        if (present(eps)) then
            epsr = eps
        else
            val_func = func(mat)
            epsr = max(abs(val_func*DIVEPS), MAX_EPS)
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

end module