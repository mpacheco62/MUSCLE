program test_I2O
    use muscle_tensors
    implicit none
    logical :: passed

    print*, "Running tests for iden_2O..."

    call test_iden_2O_sum(passed)
    if (.not. passed) STOP 1

    call test_iden_2O_sub(passed)
    if (.not. passed) STOP 2

    call test_iden_2O_mul(passed)
    if (.not. passed) STOP 3

    call test_iden_2O_div(passed)
    if (.not. passed) STOP 4

    print*, "All tests for iden_2O passed successfully!"

end program test_I2O


subroutine test_iden_2O_sum(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed
    type(iden_2O)  :: to_test1, to_test2
    type(iden_2OS) :: to_test3, result, expected

    ! 1.- Sum of two standard identities (I + I = 2.0*I)
    call expected%init(2D0)
    result = to_test1 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Sum of two iden_2O failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 2.- Sum of two different standard identity instances
    result = to_test1 + to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Sum of two iden_2O (different instances) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 3.- Sum of standard (1.0) and scaled (4.0) identities
    call expected%init(5D0)
    call to_test3%init(4D0)
    result = to_test1 + to_test3
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Sum of iden_2O and iden_2OS failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 4.- Commutativity: Sum of scaled (4.0) and standard (1.0) identities
    result = to_test3 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: Sum of iden_2OS and iden_2O failed (Commutativity)", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

end subroutine


subroutine test_iden_2O_sub(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed
    type(iden_2O)  :: to_test1, to_test2
    type(iden_2OS) :: to_test3, result, expected

    ! 1.- Subtraction of same identity (I - I = 0.0)
    call expected%init(0D0)
    result = to_test1 - to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Subtraction of same iden_2O failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 2.- Subtraction of different identity instances
    result = to_test1 - to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Subtraction of two iden_2O instances failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 3.- Unary subtraction (-I = -1.0*I)
    call expected%init(-1D0)
    result = -to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Unary minus of iden_2O failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 4.- Mixed subtraction: iden_2O (1.0) - iden_2OS (5.0) = -4.0
    call expected%init(-4D0)
    call to_test3%init(5D0)
    result = to_test1 - to_test3
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: Subtraction (iden_2O - iden_2OS) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 5.- Mixed subtraction: iden_2OS (5.0) - iden_2O (1.0) = 4.0
    call expected%init(4D0)
    result = to_test3 - to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: Subtraction (iden_2OS - iden_2O) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

end subroutine

subroutine test_iden_2O_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed
    type(iden_2O)  :: to_test1, result_I2O
    type(iden_2OS) :: result, expected

    ! 1.- Left scalar multiplication (10.0 * I)
    call expected%init(10D0)
    result = 10D0 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Multiplication (real * iden_2O) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 2.- Right scalar multiplication (I * 5.0)
    call expected%init(5D0)
    result = to_test1 * 5D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Multiplication (iden_2O * real) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 3.- Identity multiplication (I * I = I)
    ! As iden_2O has no data, we check if the operation returns the correct type 
    ! and is compatible with another iden_2O instance.
    result_I2O = to_test1 * to_test1
    ! if (.not. passed) then
    !     print*, "3.- Error: Multiplication of two iden_2O (I * I) failed"
    !     return
    ! end if

end subroutine

subroutine test_iden_2O_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed
    type(iden_2O)  :: to_test1
    type(iden_2OS) :: result, expected

    ! 1.- Scalar division (I / 10.0)
    call expected%init(0.1D0)
    result = to_test1 / 10D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Division of iden_2O by real failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

end subroutine