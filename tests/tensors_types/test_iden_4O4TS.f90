program test_iden_4O4TS
    use mod_iden_4O4TS
    implicit none
    logical :: passed

    print*, "Starting iden_4O4TS (Scaled Sym I4) tests..."

    call test_iden_4O4TS_mul(passed)
    if (.not. passed) stop 1

    call test_iden_4O4TS_div(passed)
    if (.not. passed) stop 2

    print*, "All iden_4O4TS tests passed successfully!"
    stop 0
end program test_iden_4O4TS

! ==========================================
! Subroutine: test_iden_4O4TS_mul
! ==========================================
subroutine test_iden_4O4TS_mul(passed)
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4TS
    implicit none
    logical, intent(out) :: passed
    
    type(iden_4O4TS) :: iden, res
    real(real64) :: factor, expected_val
    real(real64), parameter :: EPS=1D-10

    ! 1.- Tensor * Scalar
    iden%val = 4.0D0
    factor = 0.5D0
    expected_val = 2.0D0
    res = iden * factor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "1.- Error: Right multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 2.- Scalar * Tensor
    factor = 10.0D0
    expected_val = 40.0D0
    res = factor * iden
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "2.- Error: Left multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 3.- Multiplication with negative scaling
    iden%val = -1.0D0
    factor = 5.0D0
    expected_val = -5.0D0
    res = iden * factor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "3.- Error: Negative scaling multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if
end subroutine test_iden_4O4TS_mul

! ==========================================
! Subroutine: test_iden_4O4TS_div
! ==========================================
subroutine test_iden_4O4TS_div(passed)
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4TS
    implicit none
    logical, intent(out) :: passed
    
    type(iden_4O4TS) :: iden, res
    real(real64) :: divisor, expected_val
    real(real64), parameter :: EPS=1D-10

    ! 1.- Standard division
    iden%val = 5.0D0
    divisor = 10.0D0
    expected_val = 0.5D0
    res = iden / divisor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "1.- Error: Scalar division failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 2.- Division by negative scalar
    iden%val = 1.0D0
    divisor = -2.0D0
    expected_val = -0.5D0
    res = iden / divisor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "2.- Error: Division by negative failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if
end subroutine test_iden_4O4TS_div