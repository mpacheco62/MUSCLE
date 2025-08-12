module test_derivatives_mod
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: fun_scalar_test1, fun_scalar_test2, fun_scalar_test3, energy_hooke
    type, public :: mytype_test
        real(real64) :: a
        contains
        procedure, public :: fun => fun_object_scalar_test1
        procedure, public :: der => derivative_fun1
    end type

contains
    pure function fun_object_scalar_test1(self, x) result(res)
        use, intrinsic :: iso_fortran_env
        use tensors_types, only : ten_3D2Osym
        implicit none
        class(mytype_test), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res

        res =   self%a*(abs(x%vals(1)) + abs(x%vals(2)) + abs(x%vals(3)) &
            + 2*abs(x%vals(4)) + 2*abs(x%vals(5)) + 2*abs(x%vals(6)))
    end function

    pure function derivative_fun1(self, x) result(res)
        use, intrinsic :: iso_fortran_env
        use tensors_types, only : ten_3D2Osym
        use derivatives
        implicit none
        class(mytype_test), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res

        res = derivative(func=wrapper, mat=x)

        contains
            pure function wrapper(x1) result(res1)
                use, intrinsic :: iso_fortran_env
                use tensors_types, only : ten_3D2Osym
                implicit none
                type(ten_3D2Osym), intent(in) :: x1
                real(real64) :: res1
                res1 = self%fun(x1)
            end function wrapper
    end function



    pure function fun_scalar_test1(x) result(res)
        use, intrinsic :: iso_fortran_env
        use tensors_types, only : ten_3D2Osym
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res

        res =   abs(x%vals(1)) + abs(x%vals(2)) + abs(x%vals(3)) &
            + 2*abs(x%vals(4)) + 2*abs(x%vals(5)) + 2*abs(x%vals(6))
    end function

    pure function fun_scalar_test2(x) result(res)
        use, intrinsic :: iso_fortran_env
        use tensors_types, only : ten_3D2Osym
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
        use tensors_types, only : ten_3D2Osym
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res

        res =   abs(x%vals(1))*1 + abs(x%vals(2))*2 + abs(x%vals(3))*3 &
            + 2*abs(x%vals(4))*4 + 2*abs(x%vals(5))*5 + 2*abs(x%vals(6))*6
    end function

    pure function energy_hooke(strain) result(energy)
            use tensors_types
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
end module test_derivatives_mod


! ********************** PROGRAM TEST ************************************
program test_derivatives
    use derivatives
    implicit none
    
    logical :: passed

    call test_derivate(passed)
    if (.not. passed) STOP 1

    call test_object_derivate(passed)
    if (.not. passed) STOP 2

    call test_derivate2O(passed)
    if (.not. passed) STOP 1

    STOP 0
end program test_derivatives

subroutine test_derivate(passed)
    use, intrinsic :: iso_fortran_env
    use derivatives
    use tensors_types
    use test_derivatives_mod
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym):: to_test
    type(ten_3D2Osym) :: expected, result

    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    result = derivative(fun_scalar_test1, to_test)

    passed = expected .isequal. result

    if (.not. passed) print*, "Case 1",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    call expected%init(vals=(/2D0, 2D0, 2D0, 2D0, 2D0, 2D0/))
    result = derivative(fun_scalar_test2, to_test)
    passed = expected .isequal. result

    if (.not. passed) print*, "Case 2",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return


    call expected%init(vals=(/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    result = derivative(fun_scalar_test3, to_test)
    passed = expected .isequal. result

    if (.not. passed) print*, "Case 3",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return
end subroutine




subroutine test_object_derivate(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    use test_derivatives_mod
    implicit none
    
    logical, intent(out) :: passed

    type(mytype_test) :: obj
    type(ten_3D2Osym) :: to_test
    type(ten_3D2Osym) :: expected, result

    obj%a = 5D0
    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected%init(vals=(/5D0, 5D0, 5D0, 5D0, 5D0, 5D0/))
    result = obj%der(to_test)

    passed = expected .isequal. result

    if (.not. passed) print*, "Case 1",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return


end subroutine



subroutine test_derivate2O(passed)
    use, intrinsic :: iso_fortran_env
    use derivatives
    use tensors_types
    use test_derivatives_mod
    use mod_elasticity_linear
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym):: to_test
    type(ten_3D4O2sym) :: expected, result
    real(real64) :: E, nu, lam, mu

    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    expected%vals = 0.0D0
    result = derivative2O(fun_scalar_test1, to_test)

    passed = expected .isequal. result

    if (.not. passed) print*, "Case 1 second derivative",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    E=100D0; nu=0.3D0
    lam = E*nu/((1D0+nu)*(1D0-2D0*nu))
    mu = E/(2D0*(1D0+nu))
    call to_test%init(vals=(/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init(xxxx=2D0*mu+lam, xxyy=lam,        xxzz=lam,        xxxy=0D0, xxyz=0D0, xxxz=0D0,  &
                       yyxx=lam,        yyyy=2D0*mu+lam, yyzz=lam,        yyxy=0D0, yyyz=0D0, yyxz=0D0,  &
                       zzxx=lam,        zzyy=lam,        zzzz=2D0*mu+lam, zzxy=0D0, zzyz=0D0, zzxz=0D0,  &
                       xyxx=0D0,        xyyy=0D0,        xyzz=0D0,        xyxy=mu,  xyyz=0D0, xyxz=0D0,  &
                       yzxx=0D0,        yzyy=0D0,        yzzz=0D0,        yzxy=0D0, yzyz=mu,  yzxz=0D0,  &
                       xzxx=0D0,        xzyy=0D0,        xzzz=0D0,        xzxy=0D0, xzyz=0D0, xzxz=mu    &
                       )
    result = derivative2O(energy_hooke, to_test)

    passed = expected .isequal. result
    if (.not. passed) print*, "Case 2 second derivative",  new_line('A'), &
                              "The hook law by derivative is different to analitical", new_line('A'), &
                              "The values obtained are:", new_line('A'),  &
                              result%vals(1,:), new_line('A'), &
                              result%vals(2,:), new_line('A'), &
                              result%vals(3,:), new_line('A'), &
                              result%vals(4,:), new_line('A'), &
                              result%vals(5,:), new_line('A'), &
                              result%vals(6,:), new_line('A'), &
                              "The expected are:", new_line('A'), &
                              expected%vals(1,:), new_line('A'), &
                              expected%vals(2,:), new_line('A'), &
                              expected%vals(3,:), new_line('A'), &
                              expected%vals(4,:), new_line('A'), &
                              expected%vals(5,:), new_line('A'), &
                              expected%vals(6,:), new_line('A'), new_line('A'), &

                              "Differences:", new_line('A'), &
                              expected%vals(1,:)-result%vals(1,:), new_line('A'), &
                              expected%vals(2,:)-result%vals(2,:), new_line('A'), &
                              expected%vals(3,:)-result%vals(3,:), new_line('A'), &
                              expected%vals(4,:)-result%vals(4,:), new_line('A'), &
                              expected%vals(5,:)-result%vals(5,:), new_line('A'), &
                              expected%vals(6,:)-result%vals(6,:), new_line('A')
    if (.not. passed) return


    passed = .false.
end subroutine