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

    print*, "Hola!", passed
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

    call to_test%init(vals=(/1D0, 1D0, 2D0, 0D0, 0D0, 0D0/))
    expected = (/ 1D0, 1D0, 1D0/)
    result = eigenvals(to_test)

    print*, "Autovalores"
    print*, result

    passed = .false.
end subroutine
