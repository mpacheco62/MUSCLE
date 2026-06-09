module test_muscle_math_derivatives_mod
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: fun_scalar_scalar_test1, fun_scalar_scalar_test2, fun_scalar_scalar_test3
    public :: fun_scalar_test1, fun_scalar_test2, fun_scalar_test3, energy_hooke
    public :: fun_tens_x, fun_tens_2x, fun_tens_tenx, fun_nonlinear_tens
    type, public :: mytype_test
        real(real64) :: a
        contains
        procedure, public :: fun => fun_object_scalar_test1
        procedure, public :: der => derivative_fun1
    end type

contains

    pure function fun_object_scalar_test1(self, x) result(res)
        use, intrinsic :: iso_fortran_env
        use muscle_tensors, only : ten_3D2Osym
        implicit none
        class(mytype_test), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res

        res =   self%a*(abs(x%vals(1)) + abs(x%vals(2)) + abs(x%vals(3)) &
            + 2*abs(x%vals(4)) + 2*abs(x%vals(5)) + 2*abs(x%vals(6)))
    end function

    pure function derivative_fun1(self, x) result(res)
        use, intrinsic :: iso_fortran_env
        use muscle_tensors, only : ten_3D2Osym
        use muscle_math_derivatives
        implicit none
        class(mytype_test), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res

        res = derivative(func=wrapper, mat=x)

        contains
            pure function wrapper(x1) result(res1)
                use, intrinsic :: iso_fortran_env
                use muscle_tensors, only : ten_3D2Osym
                implicit none
                type(ten_3D2Osym), intent(in) :: x1
                real(real64) :: res1
                res1 = self%fun(x1)
            end function wrapper
    end function


    pure function fun_scalar_scalar_test1(x) result(res)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64) :: res
        real(real64), intent(in) :: x

        res = abs(x)
    end function

    pure function fun_scalar_scalar_test2(x) result(res)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64) :: res
        real(real64), intent(in) :: x

        res = 3*abs(x)
    end function

    pure function fun_scalar_scalar_test3(x) result(res)
        use, intrinsic :: iso_fortran_env
        implicit none
        real(real64) :: res
        real(real64), intent(in) :: x

        res = sin(x)
    end function

    pure function fun_scalar_test1(x) result(res)
        use, intrinsic :: iso_fortran_env
        use muscle_tensors, only : ten_3D2Osym
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res

        res =   abs(x%vals(1)) + abs(x%vals(2)) + abs(x%vals(3)) &
            + 2*abs(x%vals(4)) + 2*abs(x%vals(5)) + 2*abs(x%vals(6))
    end function

    pure function fun_scalar_test2(x) result(res)
        use, intrinsic :: iso_fortran_env
        use muscle_tensors, only : ten_3D2Osym
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: y
        real(real64) :: res

        y%vals = 2*x%vals

        res =   abs(y%vals(1)) + abs(y%vals(2)) + abs(y%vals(3)) &
            + 2*abs(y%vals(4)) + 2*abs(y%vals(5)) + 2*abs(y%vals(6))
    end function

    pure function fun_scalar_test3(x) result(res)
        use, intrinsic :: iso_fortran_env
        use muscle_tensors, only : ten_3D2Osym
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res

        res =   abs(x%vals(1))*1 + abs(x%vals(2))*2 + abs(x%vals(3))*3 &
            + 2*abs(x%vals(4))*4 + 2*abs(x%vals(5))*5 + 2*abs(x%vals(6))*6
    end function

    pure function energy_hooke(strain) result(energy)
            use muscle_tensors
            implicit none
            type(ten_3D2Osym), intent(in) :: strain
            real(real64) :: energy
            type(ten_3D2Osym) :: strain2
            real(real64) :: E, nu, lam, mu

            E=100D0; nu=0.3D0
            lam = E*nu/((1D0+nu)*(1D0-2D0*nu))
            mu = E/(2D0*(1D0+nu))
            strain2 = strain%square()
            energy = mu*sum(strain2%vals(1:3)) + lam/2D0*sum(strain%vals(1:3))**2
    end function


    pure function fun_tens_x(x) result(res)
        use muscle_tensors, only : ten_3D2Osym
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res
        res = x
    end function

    pure function fun_tens_2x(x) result(res)
        use muscle_tensors, only : ten_3D2Osym, operator(*)
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res
        res = 2D0*x
    end function

    pure function fun_tens_tenx(x) result(res)
        use muscle_tensors, only : ten_3D2Osym, ten_3D4O3sym, operator(*), operator(.ddot.)
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res
        type(ten_3D4O3sym) :: C
        call C%init((/1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0,11D0,12D0,13D0,14D0,15D0,16D0,17D0,18D0,19D0,20D0,21D0/))
        res = C.ddot.x
    end function

    ! f(x) = tr(x) * x
    pure function fun_nonlinear_tens(x) result(res)
        use muscle_tensors, only : ten_3D2Osym
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res
        real(real64) :: trace        
        trace = x%vals(1) + x%vals(2) + x%vals(3)
        res%vals = trace * x%vals
    end function
end module test_muscle_math_derivatives_mod


! ********************** PROGRAM TEST ************************************
program test_muscle_math_derivatives
    use muscle_math_derivatives
    implicit none
    
    logical :: passed

    call test_derivative_scalar_scalar(passed)
    if (.not. passed) STOP 1

    call test_derivate_scalar_ten(passed)
    if (.not. passed) STOP 2

    call test_derivate_ten_ten(passed)
    if (.not. passed) STOP 3

    call test_object_derivate_scalar_ten(passed)
    if (.not. passed) STOP 4

    call test_derivate2O_scalar_ten(passed)
    if (.not. passed) STOP 5

    STOP 0
end program test_muscle_math_derivatives


subroutine test_derivative_scalar_scalar(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_math_derivatives
    use test_muscle_math_derivatives_mod
    implicit none
    
    logical, intent(out) :: passed

    real(real64) :: to_test
    real(real64) :: expected, result
    real(real64) :: pi
    pi = 4.0D0 * atan(1.0D0)
    

    to_test = 5D0
    expected = 1.0D0
    result = derivative(fun_scalar_scalar_test1, to_test)

    passed = abs(expected - result) < 1.0D-7

    if (.not. passed) print*, "Case 1 Derivate fun(scalar)=> scalar",  new_line('A'), &
                              "The value obtained is different from the expected one", new_line('A'), &
                              "The value obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    to_test = -5D0
    expected = -1.0D0
    result = derivative(fun_scalar_scalar_test1, to_test)

    passed = abs(expected - result) < 1.0D-7

    if (.not. passed) print*, "Case 2 Derivate fun(scalar)=> scalar",  new_line('A'), &
                              "The value obtained is different from the expected one", new_line('A'), &
                              "The value obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    to_test = 5D0
    expected = 3D0
    result = derivative(fun_scalar_scalar_test2, to_test)

    passed = abs(expected - result) < 1.0D-7

    if (.not. passed) print*, "Case 3 Derivate fun(scalar)=> scalar",  new_line('A'), &
                              "The value obtained is different from the expected one", new_line('A'), &
                              "The value obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    to_test = -5D0
    expected = -3D0
    result = derivative(fun_scalar_scalar_test2, to_test)

    passed = abs(expected - result) < 1.0D-7

    if (.not. passed) print*, "Case 4 Derivate fun(scalar)=> scalar",  new_line('A'), &
                              "The value obtained is different from the expected one", new_line('A'), &
                              "The value obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    to_test = pi
    expected = -1D0
    result = derivative(fun_scalar_scalar_test3, to_test)

    passed = abs(expected - result) < 1.0D-7

    if (.not. passed) print*, "Case 5 Derivate fun(scalar)=> scalar",  new_line('A'), &
                              "The value obtained is different from the expected one", new_line('A'), &
                              "The value obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return


    to_test = 0
    expected = 1D0
    result = derivative(fun_scalar_scalar_test3, to_test)

    passed = abs(expected - result) < 1.0D-7

    if (.not. passed) print*, "Case 6 Derivate fun(scalar)=> scalar",  new_line('A'), &
                              "The value obtained is different from the expected one", new_line('A'), &
                              "The value obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    ! call expected%init(vals=(/2D0, 2D0, 2D0, 2D0, 2D0, 2D0/))
    ! result = derivative(fun_scalar_test2, to_test)
    ! passed = expected .approx. result

    ! if (.not. passed) print*, "Case 2 Derivate",  new_line('A'), &
    !                           "The eigenvalues obtained is different from the expected one", new_line('A'), &
    !                           "The values obtained are:", result, new_line('A'), &
    !                           "The expected are:", expected
    ! if (.not. passed) return


    ! call expected%init(vals=(/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    ! result = derivative(fun_scalar_test3, to_test)
    ! passed = expected .approx. result

    ! if (.not. passed) print*, "Case 3 Derivate",  new_line('A'), &
    !                           "The eigenvalues obtained is different from the expected one", new_line('A'), &
    !                           "The values obtained are:", result, new_line('A'), &
    !                           "The expected are:", expected
    ! if (.not. passed) return
end subroutine


subroutine test_derivate_scalar_ten(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_math_derivatives
    use muscle_tensors
    use test_muscle_math_derivatives_mod
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym):: to_test
    type(ten_3D2Osym) :: expected, result

    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    result = derivative(fun_scalar_test1, to_test)
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 1 Derivate",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected
    if (.not. passed) return

    call expected%init(vals=(/2D0, 2D0, 2D0, 2D0, 2D0, 2D0/))
    result = derivative(fun_scalar_test2, to_test)
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 2 Derivate",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected
    if (.not. passed) return


    call expected%init(vals=(/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    result = derivative(fun_scalar_test3, to_test)
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 3 Derivate",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return
end subroutine


subroutine test_derivate_ten_ten(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_math_derivatives
    use muscle_tensors
    use test_muscle_math_derivatives_mod
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test
    type(iden_4O4T) :: I4O4T
    type(iden_2O) :: I2O
    type(ten_3D4O2sym) :: expected, result

    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    expected = I4O4T
    result = derivative(fun_tens_x, to_test)
    passed = expected .approx. result

    if (.not. passed) print*, "Case 1 Derivate tensor=>tensor",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    to_test%vals = (/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/)
    expected = I4O4T
    result = derivative(fun_tens_x, to_test)
    passed = expected .approx. result

    if (.not. passed) print*, "Case 2 Derivate tensor=>tensor",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return


    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    expected = 2D0*I4O4T
    result = derivative(fun_tens_2x, to_test)
    passed = expected .approx. result

    if (.not. passed) print*, "Case 3 Derivate tensor=>tensor",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    to_test%vals = (/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/)
    expected = 2D0*I4O4T
    result = derivative(fun_tens_2x, to_test) 
    passed = expected .approx. result

    if (.not. passed) print*, "Case 4 Derivate tensor=>tensor",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return


    call expected%init(xxxx=1D0,   yyyy=2D0,   zzzz=3D0,  &
                       xyxy=4D0,   yzyz=5D0,   xzxz=6D0,  &
                       xxyy=7D0,   yyzz=8D0,              &
                       yyxx=7D0,   zzyy=8D0,              &
                       zzxy=9D0,   xyyz=10D0,  yzxz=11D0, &
                       xyzz=9D0,   yzxy=10D0,  xzyz=11D0, &
                       xxzz=12D0,                         &
                       zzxx=12D0,                         &
                       yyxy=13D0,  zzyz=14D0,  xyxz=15D0, &
                       xyyy=13D0,  yzzz=14D0,  xzxy=15D0, &
                       xxxy=16D0,  yyyz=17D0,  zzxz=18D0, &
                       xyxx=16D0,  yzyy=17D0,  xzzz=18D0, &
                       xxyz=19D0,  yyxz=20D0,  xxxz=21D0, &
                       yzxx=19D0,  xzyy=20D0,  xzxx=21D0  &
                       )
    result = derivative(fun_tens_tenx, to_test)
    passed = expected .approx. result

    if (.not. passed) print*, "Case 5 Derivate tensor=>tensor",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return


    expected = sum(to_test%vals(1:3))*I4O4T 
    expected = expected + (to_test .tdot. I2O)  ! TODO suma de I4O3TS y I4O4TS
    result = derivative(fun_nonlinear_tens, to_test)
    passed = expected .approx. result

    if (.not. passed) print*, "Case 6 Derivate tensor=>tensor",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

end subroutine



subroutine test_object_derivate_scalar_ten(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use test_muscle_math_derivatives_mod
    implicit none
    
    logical, intent(out) :: passed

    type(mytype_test) :: obj
    type(ten_3D2Osym) :: to_test
    type(ten_3D2Osym) :: expected, result

    obj%a = 5D0
    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected%init(vals=(/5D0, 5D0, 5D0, 5D0, 5D0, 5D0/))
    result = obj%der(to_test)

    passed = result%is_approx(expected, tol=1.0D-10)


    if (.not. passed) print*, "Case 1 Derivate object",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return


end subroutine



subroutine test_derivate2O_scalar_ten(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_math_derivatives
    use muscle_tensors
    use test_muscle_math_derivatives_mod
    use mod_elasticity_linear
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym):: to_test
    type(ten_3D4O3sym) :: expected, result
    real(real64) :: E, nu, lam, mu

    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    expected%vals = 0.0D0
    result = derivative2O(fun_scalar_test1, to_test)

    passed = expected%is_approx(result, tol=1.0D-8)

    if (.not. passed) print*, "Case 1 second derivative",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    E=100D0; nu=0.3D0
    lam = E*nu/((1D0+nu)*(1D0-2D0*nu))
    mu = E/(2D0*(1D0+nu))
    call to_test%init(vals=(/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init(xxxx=2D0*mu+lam, yyyy=2D0*mu+lam, zzzz=2D0*mu+lam, &
                       xyxy=mu,         yzyz=mu,         xzxz=mu,         &
                       xxyy=lam,        yyzz=lam,                         &
                       zzxy=0D0,        xyyz=0D0,        yzxz=0D0,        &
                       xxzz=lam,                                          &
                       yyxy=0D0,        zzyz=0D0,        xyxz=0D0,        &
                       xxxy=0D0,        yyyz=0D0,        zzxz=0D0,        &
                       xxyz=0D0,        yyxz=0D0,        xxxz=0D0         &
                       )
    result = derivative2O(energy_hooke, to_test)

    passed = expected%is_approx(result, tol=1.0D-8)
    if (.not. passed) print*, "Case 2 second derivative",  new_line('A'), &
                              "The hook law by derivative is different to analitical", new_line('A'), &
                              "The values obtained are:", new_line('A'),  &
                              result%vals(:), new_line('A'), &
                              "The expected are:", new_line('A'), &
                              expected%vals(:), new_line('A')
                             
    if (.not. passed) return
end subroutine