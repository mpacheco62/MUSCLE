program test_3D2Osym
    use tensors_types
    implicit none
    
    logical :: passed

    call test_ten_3D2Osym_isequal(passed)
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

    print*, "Hola!", passed
    STOP 0
end program test_3D2Osym

subroutine test_ten_3D2Osym_isequal(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2

    call to_test1%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0)
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0)

    passed = to_test1 .isequal. to_test2
    if (.not. passed) return

    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6.0001D0)
    passed = .not. (to_test1 .isequal. to_test2)
    if (.not. passed) return

    call to_test1%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0)
    call to_test2%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0)
    passed = to_test1 .isequal. to_test2
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
    passed = (to_test1 + to_test1) .isequal. to_test1
    if (.not. passed) return
    
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test1 + to_test2) .isequal. to_test2
    if (.not. passed) return

    passed = (to_test2 + to_test1) .isequal. to_test2
    if (.not. passed) return

    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    call expected_result%init((/12D0, 14D0, 16D0, 18D0, 20D0, 22D0/))
    passed = (to_test2 + to_test1) .isequal. expected_result
    if (.not. passed) return

    passed = (to_test1 + to_test2) .isequal. expected_result
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
    passed = (to_test1 - to_test1) .isequal. to_test1
    if (.not. passed) return
    
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test2 - to_test1) .isequal. to_test2
    if (.not. passed) return

    passed = (to_test1 - to_test2) .isequal. (-to_test2)
    if (.not. passed) return

    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    call expected_result%init((/10D0, 10D0, 10D0, 10D0, 10D0, 10D0/))
    passed = (to_test1 - to_test2) .isequal. expected_result
    if (.not. passed) return

    passed = (to_test2 - to_test1) .isequal. (-expected_result)
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
    passed = (10D0*to_test1) .isequal. to_test1
    if (.not. passed) return

    passed = (to_test1*10D0) .isequal. to_test1
    if (.not. passed) return

    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected_result%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (0D0*to_test1) .isequal. expected_result
    if (.not. passed) return

    passed = (to_test1*0D0) .isequal. expected_result
    if (.not. passed) return

    call expected_result%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0/))
    passed = (2D0*to_test1) .isequal. expected_result
    if (.not. passed) return

    passed = (to_test1*2D0) .isequal. expected_result
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
    passed = (to_test1/2D0) .isequal. to_test1
    if (.not. passed) return

    call to_test1%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0/))
    call expected_result%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test1/2D0) .isequal. expected_result
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
    passed = (.dev. to_test1) .isequal. to_test1
    if (.not. passed) return

    call to_test1%init((/5D0, 5D0, 5D0, 0D0, 0D0, 0D0/))
    call expected_result%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (.dev. to_test1) .isequal. expected_result
    if (.not. passed) return

    call to_test1%init((/5D0, 5D0, 5D0, 1D0, 1D0, 1D0/))
    call expected_result%init((/0D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    passed = (.dev. to_test1) .isequal. expected_result
    if (.not. passed) return

    call to_test1%init((/3D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    call expected_result%init((/2D0, -1D0, -1D0, 1D0, 1D0, 1D0/))
    passed = (.dev. to_test1) .isequal. expected_result
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