program test_basic_operations
    use basic_operations
    implicit none
    
    logical :: passed

    call test_eigenvalues(passed)
    if (.not. passed) STOP 1

    ! call test_deviatoric_3x3_1(passed)
    ! if (.not. passed) STOP 1

    ! call test_deviatoric_6_2(passed)
    ! if (.not. passed) STOP 2

    ! call test_dp_sym_3x3_3(passed)
    ! if (.not. passed) STOP 3

    ! call test_dp_sym_6_4(passed)
    ! if (.not. passed) STOP 4

    ! print*, "Hola!", passed
    STOP 0
end program test_basic_operations

subroutine test_eigenvalues(passed)
    use, intrinsic :: iso_fortran_env
    use basic_operations
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym):: to_test
    real(real64), dimension(3) :: expected, result

    call to_test%init(vals=(/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    expected = (/ 1D0, 1D0, 1D0/)
    result = eigenvals(to_test)

    passed = abs(norm2(expected-result)) .le. 1D-7

    if (.not. passed) print*, "Case 1",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return
    


    call to_test%init(vals=(/1D0, 1D0, 2D0, 0D0, 0D0, 0D0/))
    expected = (/ 2D0, 1D0, 1D0/)
    result = eigenvals(to_test)

    passed = abs(norm2(expected-result)) .le. 1D-7

    if (.not. passed) print*, "Case 2",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 0D0, 0D0/))
    expected = (/ 2D0, 1D0, 0D0/)
    result = eigenvals(to_test)

    passed = abs(norm2(expected-result)) .le. 1D-7

    if (.not. passed) print*, "Case 3",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    call to_test%init(vals=(/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    expected = (/ 3D0, 0D0, 0D0/)
    result = eigenvals(to_test)

    passed = abs(norm2(expected-result)) .le. 1D-7

    if (.not. passed) print*, "Case 4",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    call to_test%init(vals=(/1D0, 1D0, 1D0, 2D0, 2D0, 2D0/))
    expected = (/ 5D0, -1D0, -1D0/)
    result = eigenvals(to_test)

    passed = abs(norm2(expected-result)) .le. 1D-7

    if (.not. passed) print*, "Case 5",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

    call to_test%init(vals=(/1D0, 2D0, 3D0, 1D0, 1D0, 1D0/))
    expected = (/ 4.21431974337753D0, 1.46081112718911D0, 0.32486912943D0/)
    result = eigenvals(to_test)

    passed = abs(norm2(expected-result)) .le. 1D-7

    if (.not. passed) print*, "Case 6",  new_line('A'), &
                              "The eigenvalues obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected
    if (.not. passed) return

end subroutine