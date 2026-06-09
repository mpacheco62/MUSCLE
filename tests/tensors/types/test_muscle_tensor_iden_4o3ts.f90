program test_muscle_tensor_iden_4o3ts
    use muscle_tensors
    implicit none
    logical :: passed

    print*, "Starting iden_4O3TS tests..."

    call test_iden_4O3TS_mul(passed)
    if (.not. passed) stop 1

    call test_iden_4O3TS_div(passed)
    if (.not. passed) stop 2

    print*, "All iden_4O3TS tests passed successfully!"
    stop 0
end program test_muscle_tensor_iden_4o3ts

! ==========================================
! Subroutine: test_iden_4O3TS_mul
! ==========================================
subroutine test_iden_4O3TS_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed
    
    type(iden_4O3TS) :: iden, res
    real(real64) :: factor, expected_val
    real(real64), parameter :: EPS=1D-10

    ! 1.- Scalar * Tensor (Left multiplication)
    iden%val = 5.0D0
    factor = 2.0D0
    expected_val = 10.0D0
    res = factor * iden
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "1.- Error: Left multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 2.- Tensor * Scalar (Right multiplication)
    factor = -3.0D0
    expected_val = -15.0D0
    res = iden * factor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "2.- Error: Right multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 3.- Multiplication by zero
    factor = 0.0D0
    expected_val = 0.0D0
    res = iden * factor
    passed = abs(res%val - expected_val) < 1D-12
    if (.not. passed) then
        print*, "3.- Error: Multiplication by zero failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if
end subroutine test_iden_4O3TS_mul

! ==========================================
! Subroutine: test_iden_4O3TS_div
! ==========================================
subroutine test_iden_4O3TS_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed
    
    type(iden_4O3TS) :: iden, res
    real(real64) :: divisor, expected_val
    real(real64), parameter :: EPS=1D-10

    ! 1.- Standard scalar division
    iden%val = 10.0D0
    divisor = 4.0D0
    expected_val = 2.5D0
    res = iden / divisor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "1.- Error: Scalar division failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 2.- Division of negative value
    iden%val = -100.0D0
    divisor = 10.0D0
    expected_val = -10.0D0
    res = iden / divisor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "2.- Error: Division of negative value failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if
end subroutine test_iden_4O3TS_div