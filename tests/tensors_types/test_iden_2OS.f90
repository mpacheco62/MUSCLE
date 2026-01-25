program test_I2OS
    use tensors_types
    implicit none
    logical :: passed

    ! real*8 :: a
    ! call test123(a, 5D0)

    call test_iden_2OS_sum(passed)
    if (.not. passed) STOP 1

    call test_iden_2OS_sub(passed)
    if (.not. passed) STOP 2

    call test_iden_2OS_mul(passed)
    if (.not. passed) STOP 3

    call test_iden_2OS_div(passed)
    if (.not. passed) STOP 4

    STOP 0
end program test_I2OS


subroutine test_iden_2OS_sum(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    type(iden_2OS) :: to_test1, to_test2
    type(iden_2OS) :: expected, result

    call to_test1%init(0D0)
    call expected%init(0D0)
    result = to_test1 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Sum of two zero-valued iden_2OS (0 + 0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if
    
    call to_test2%init(1D0)
    call expected%init(1D0)
    result = (to_test1 + to_test2)
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Sum of iden_2OS (0.0 + 1.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

    result = (to_test2 + to_test1)
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Sum of iden_2OS (1.0 + 0.0) failed (Commutativity)", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

    call to_test1%init(2D0)
    call expected%init(3D0)
    result = (to_test1 + to_test2) 
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: Sum of iden_2OS (2.0 + 1.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

    call expected%init(3D0)
    result = (to_test2 + to_test1) 
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: Sum of iden_2OS (1.0 + 2.0) failed (Commutativity)", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

end subroutine

subroutine test_iden_2OS_sub(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    type(iden_2OS) :: to_test1, to_test2
    type(iden_2OS) :: result, expected

    ! 1.- Subtraction of two zero-valued scaled identities (0.0 - 0.0)
    call to_test1%init(0D0)
    call expected%init(0D0)
    result = to_test1 - to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Subtraction of zero-valued iden_2OS (0.0 - 0.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if
    
    ! 2.- Subtraction: 1.0 - 0.0
    call to_test2%init(1D0)
    call expected%init(1D0)
    result = to_test2 - to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Subtraction of iden_2OS (1.0 - 0.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

    ! 3.- Subtraction: 0.0 - 1.0 (Negative result check)
    call expected%init(-1D0)
    result = to_test1 - to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Subtraction of iden_2OS (0.0 - 1.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

    ! 4.- Subtraction: 3.0 - 1.0
    call to_test1%init(3D0)
    call expected%init(2D0)
    result = to_test1 - to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: Subtraction of iden_2OS (3.0 - 1.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

    ! 5.- Unary Minus: -(3.0)
    call expected%init(-3D0)
    result = -to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: Unary minus of iden_2OS (-3.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

end subroutine

subroutine test_iden_2OS_mul(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    type(iden_2OS) :: to_test1, to_test2
    type(iden_2OS) :: expected, result

    ! 1.- Scalar multiplication: real * iden_2OS (10.0 * 0.0)
    call to_test1%init(0D0)
    call expected%init(0D0)
    result = 10D0 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Multiplication (real * iden_2OS) with zero failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 2.- Scalar multiplication: iden_2OS * real (0.0 * 10.0)
    result = to_test1 * 10D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Multiplication (iden_2OS * real) with zero failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 3.- Multiplication by zero: iden_2OS * 0.0
    call to_test1%init(1D0)
    call expected%init(0D0)
    result = to_test1 * 0D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Multiplication by zero (iden_2OS * 0.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 4.- Multiplication by zero: 0.0 * iden_2OS
    result = 0D0 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: Multiplication by zero (0.0 * iden_2OS) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 5.- Scalar multiplication: 2.0 * iden_2OS (2.0 * 1.0)
    call expected%init(2D0)
    result = 2D0 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: Scalar multiplication (2.0 * iden_2OS) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 6.- Scalar multiplication: iden_2OS * 2.0 (1.0 * 2.0)
    result = to_test1 * 2D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "6.- Error: Scalar multiplication (iden_2OS * 2.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 7.- Multiplication of two scaled identities (5.0 * 6.0)
    call to_test1%init(5D0)
    call to_test2%init(6D0)
    call expected%init(30D0)
    result = to_test1 * to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "7.- Error: Multiplication of two iden_2OS (5.0 * 6.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

end subroutine

subroutine test_iden_2OS_div(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    type(iden_2OS) :: to_test1
    type(iden_2OS) :: result, expected

    ! 1.- Division of zero-valued identity: 0.0 / 10.0
    call to_test1%init(0D0)
    call expected%init(0D0)
    result = to_test1 / 10D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Division of zero-valued iden_2OS (0.0 / 10.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 2.- Division by identity: 2.0 / 1.0
    call to_test1%init(2D0)
    call expected%init(2D0)
    result = to_test1 / 1D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Division of iden_2OS by 1.0 (2.0 / 1.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val
        return
    end if

    ! 3.- Standard division: 2.0 / 2.0
    call expected%init(1D0)
    result = to_test1 / 2D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Division of iden_2OS (2.0 / 2.0) failed", new_line('A'), &
                "Expected:", expected%val, " Actual:", result%val, " Diff:", result%val - expected%val
        return
    end if

end subroutine
