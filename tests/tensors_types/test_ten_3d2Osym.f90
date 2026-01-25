program test_3D2Osym
    use tensors_types
    implicit none
    
    logical :: passed

    call test_ten_3D2Osym_approx(passed)
    if (.not. passed) STOP 1

    call test_ten_3D2Osym_sum(passed)
    if (.not. passed) STOP 2

    call test_ten_3D2Osym_sub(passed)
    if (.not. passed) STOP 3

    call test_ten_3D2Osym_mul(passed)
    if (.not. passed) STOP 4

    call test_ten_3D2Osym_div(passed)
    if (.not. passed) STOP 5

    call test_ten_3D2Osym_dev(passed)
    if (.not. passed) STOP 6

    call test_ten_3D2Osym_ddot(passed)
    if (.not. passed) STOP 7

    call test_ten_3D2Osym_components(passed)
    if (.not. passed) STOP 8

    call test_ten_3D2Osym_dot(passed)
    if (.not. passed) STOP 9

    call test_ten_3D2Osym_square(passed)
    if (.not. passed) STOP 10

    print*, "Hola!", passed
    STOP 0
end program test_3D2Osym

subroutine test_ten_3D2Osym_approx(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2

    call to_test1%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0)
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0)

    passed = to_test1 .approx. to_test2
    if (.not. passed) return

    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6.0001D0)
    passed = .not. (to_test1 .approx. to_test2)
    if (.not. passed) return

    call to_test1%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0)
    call to_test2%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0)
    passed = to_test1 .approx. to_test2
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_sum(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    type(ten_3D2Osym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (to_test1 + to_test1) .approx. to_test1
    if (.not. passed) return
    
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test1 + to_test2) .approx. to_test2
    if (.not. passed) return

    passed = (to_test2 + to_test1) .approx. to_test2
    if (.not. passed) return

    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    call expected_result%init((/12D0, 14D0, 16D0, 18D0, 20D0, 22D0/))
    passed = (to_test2 + to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1 + to_test2) .approx. expected_result
    if (.not. passed) return

end subroutine

subroutine test_ten_3D2Osym_sub(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    type(ten_3D2Osym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (to_test1 - to_test1) .approx. to_test1
    if (.not. passed) return
    
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test2 - to_test1) .approx. to_test2
    if (.not. passed) return

    passed = (to_test1 - to_test2) .approx. (-to_test2)
    if (.not. passed) return

    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    call expected_result%init((/10D0, 10D0, 10D0, 10D0, 10D0, 10D0/))
    passed = (to_test1 - to_test2) .approx. expected_result
    if (.not. passed) return

    passed = (to_test2 - to_test1) .approx. (-expected_result)
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_mul(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (10D0*to_test1) .approx. to_test1
    if (.not. passed) return

    passed = (to_test1*10D0) .approx. to_test1
    if (.not. passed) return

    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected_result%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (0D0*to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1*0D0) .approx. expected_result
    if (.not. passed) return

    call expected_result%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0/))
    passed = (2D0*to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1*2D0) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_div(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (to_test1/2D0) .approx. to_test1
    if (.not. passed) return

    call to_test1%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0/))
    call expected_result%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test1/2D0) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_dev(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (.dev. to_test1) .approx. to_test1
    if (.not. passed) return

    call to_test1%init((/5D0, 5D0, 5D0, 0D0, 0D0, 0D0/))
    call expected_result%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (.dev. to_test1) .approx. expected_result
    if (.not. passed) return

    call to_test1%init((/5D0, 5D0, 5D0, 1D0, 1D0, 1D0/))
    call expected_result%init((/0D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    passed = (.dev. to_test1) .approx. expected_result
    if (.not. passed) return

    call to_test1%init((/3D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    call expected_result%init((/2D0, -1D0, -1D0, 1D0, 1D0, 1D0/))
    passed = (.dev. to_test1) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_ddot(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    
    passed = .false.

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    if ((to_test1 .ddot. to_test1) < 1D-10) passed = .true.
    if (.not. passed) return

    passed = .false.
    if ((to_test1 .ddot. to_test2) < 1D-10) passed = .true.
    if (.not. passed) return

    passed = .false.
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    if (abs((to_test1 .ddot. to_test2) - 528D0) < 1D-10) passed = .true.
    if (.not. passed) return

    passed = .false.
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    if (abs((to_test2 .ddot. to_test1) - 528D0) < 1D-10) passed = .true.
    print*, to_test1 .ddot. to_test2
    if (.not. passed) return
end subroutine



subroutine test_ten_3D2Osym_components(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    real(real64) :: expected
    
    passed = .false.

    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))

    expected = 1D0
    if (abs(to_test1%xx() - expected) < 1D-10) passed = .true.
    if (.not. passed) return

    expected = 2D0
    if (abs(to_test1%yy() - expected) < 1D-10) passed = .true.
    if (.not. passed) return

    expected = 3D0
    if (abs(to_test1%zz() - expected) < 1D-10) passed = .true.
    if (.not. passed) return

    expected = 4D0
    if (abs(to_test1%xy() - expected) < 1D-10) passed = .true.
    if (.not. passed) return

    expected = 5D0
    if (abs(to_test1%yz() - expected) < 1D-10) passed = .true.
    if (.not. passed) return

    expected = 6D0
    if (abs(to_test1%xz() - expected) < 1D-10) passed = .true.
    if (.not. passed) return

end subroutine


subroutine test_ten_3D2Osym_dot(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    type(ten_3D2O) :: expected, result
    
    passed = .false.

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) return


    call to_test1%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    call expected%init((/1D0, 0D0, 0D0, &
                         0D0, 1D0, 0D0, &
                         0D0, 0D0, 1D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) return

    call to_test1%init((/2D0, 2D0, 2D0, 0D0, 0D0, 0D0/))
    call expected%init((/4D0, 0D0, 0D0, &
                         0D0, 4D0, 0D0, &
                         0D0, 0D0, 4D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) return

    call to_test1%init((/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected%init((/3D0, 3D0, 3D0, &
                         3D0, 3D0, 3D0, &
                         3D0, 3D0, 3D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) return

    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected%init((/53D0, 42D0, 44D0, &
                         42D0, 45D0, 49D0, &
                         44D0, 49D0, 70D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) return

    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call to_test2%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init((/0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0/))
    result = to_test1 * to_test2
    passed = result .approx. expected
    if (.not. passed) return

    result = to_test2 * to_test1
    passed = result .approx. expected
    if (.not. passed) return


    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call to_test2%init((/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected%init((/11D0, 11D0, 14D0, &
                         11D0, 11D0, 14D0, &
                         11D0, 11D0, 14D0/))
    result = to_test1 * to_test2
    passed = result .approx. expected
    if (.not. passed) return

    call expected%init((/11D0, 11D0, 11D0, &
                         11D0, 11D0, 11D0, &
                         14D0, 14D0, 14D0/))
    result = to_test2 * to_test1
    passed = result .approx. expected
    if (.not. passed) return

end subroutine


subroutine test_ten_3D2Osym_square(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, expected, result
    
    passed = .false.

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result = to_test1%square()
    passed = result .approx. expected
    if (.not. passed) return

    call to_test1%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    call expected%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    result = to_test1%square()
    passed = result .approx. expected
    if (.not. passed) return

    call to_test1%init((/2D0, 2D0, 2D0, 0D0, 0D0, 0D0/))
    call expected%init((/4D0, 4D0, 4D0, 0D0, 0D0, 0D0/))
    result = to_test1%square()
    passed = result .approx. expected
    if (.not. passed) return

    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected%init((/53D0, 45D0, 70D0, 42D0, 49D0, 44D0/))
    result = to_test1%square()
    passed = result .approx. expected
    if (.not. passed) return

end subroutine