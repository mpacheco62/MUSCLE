program test_iden_4O3T
    use muscle_tensor_iden_4o3t
    implicit none
    logical :: passed

    print*, "Starting iden_4O3T (Symbolic) tests..."

    call test_iden_4O3T_mul(passed)
    if (.not. passed) stop 1

    call test_iden_4O3T_div(passed)
    if (.not. passed) stop 2

    print*, "All iden_4O3T tests passed successfully!"
    stop 0
end program test_iden_4O3T

! ==========================================
! Subroutine: test_iden_4O3T_mul
! ==========================================
subroutine test_iden_4O3T_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensor_iden_4o3t
    use muscle_tensor_iden_4o3ts
    implicit none
    logical, intent(out) :: passed
    
    type(iden_4O3T)  :: id_sym
    type(iden_4O3TS) :: res
    real(real64)     :: factor, expected_val
    real(real64), parameter :: EPS=1D-10


    ! 1.- Left multiplication: real64 * iden_4O3T -> iden_4O3TS
    factor = 15.5D0
    expected_val = 15.5D0
    res = factor * id_sym
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "1.- Error: Left symbolic multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 2.- Right multiplication: iden_4O3T * real64 -> iden_4O3TS
    factor = -2.0D0
    expected_val = -2.0D0
    res = id_sym * factor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "2.- Error: Right symbolic multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 3.- Multiplication by zero
    factor = 0.0D0
    expected_val = 0.0D0
    res = id_sym * factor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "3.- Error: Symbolic multiplication by zero failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if
end subroutine test_iden_4O3T_mul

! ==========================================
! Subroutine: test_iden_4O3T_div
! ==========================================
subroutine test_iden_4O3T_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensor_iden_4o3t
    use muscle_tensor_iden_4o3ts
    implicit none
    logical, intent(out) :: passed
    
    type(iden_4O3T)  :: id_sym
    type(iden_4O3TS) :: res
    real(real64)     :: divisor, expected_val
    real(real64), parameter :: EPS=1D-10

    ! 1.- Symbolic division: iden_4O3T / real64 -> iden_4O3TS (val = 1/a)
    divisor = 4.0D0
    expected_val = 0.25D0
    res = id_sym / divisor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "1.- Error: Symbolic division failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 2.- Division by one (Identity check)
    divisor = 1.0D0
    expected_val = 1.0D0
    res = id_sym / divisor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "2.- Error: Symbolic division by 1 failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if
end subroutine test_iden_4O3T_div