module test_derivatives_mod
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: fun_scalar_test1, fun_scalar_test2, fun_scalar_test3
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

    ! call test_derivate2O(passed)
    ! if (.not. passed) STOP 1

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
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym):: to_test
    type(ten_3D4O2sym) :: expected, result

    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    expected%vals = 0.0D0
    result = derivative2O(fun_scalar_test1, to_test)

    passed = expected .isequal. result

    if (.not. passed) print*, "Case 1 second derivative",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

end subroutine