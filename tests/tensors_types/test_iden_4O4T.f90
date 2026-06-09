program test_iden_4O4T
    use muscle_tensor_iden_4o4t
    implicit none
    logical :: passed

    print*, "Starting iden_4O4T (Symbolic Sym I4) tests..."

    call test_iden_4O4T_mul(passed)
    if (.not. passed) stop 1

    call test_iden_4O4T_div(passed)
    if (.not. passed) stop 2

    print*, "All iden_4O4T tests passed successfully!"
    stop 0
end program test_iden_4O4T

! ==========================================
! Subroutine: test_iden_4O4T_mul
! ==========================================
subroutine test_iden_4O4T_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensor_iden_4o4t
    use muscle_tensor_iden_4o4ts
    implicit none
    logical, intent(out) :: passed
    
    type(iden_4O4T)  :: id_sym
    type(iden_4O4TS) :: res
    real(real64)     :: factor, expected_val
    real(real64), parameter :: EPS=1D-10

    ! 1.- Left multiplication: real64 * iden_4O4T
    factor = 3.14159D0
    expected_val = 3.14159D0
    res = factor * id_sym
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "1.- Error: Left symbolic multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 2.- Right multiplication: iden_4O4T * real64
    factor = 0.001D0
    expected_val = 0.001D0
    res = id_sym * factor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "2.- Error: Right symbolic multiplication failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if
end subroutine test_iden_4O4T_mul

! ==========================================
! Subroutine: test_iden_4O4T_div
! ==========================================
subroutine test_iden_4O4T_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensor_iden_4o4t
    use muscle_tensor_iden_4o4ts
    implicit none
    logical, intent(out) :: passed
    real(real64), parameter :: EPS=1D-10
    
    type(iden_4O4T)  :: id_sym
    type(iden_4O4TS) :: res
    real(real64)     :: divisor, expected_val

    ! 1.- Symbolic division: I^sym / a -> res%val = 1/a
    divisor = 2.0D0
    expected_val = 0.5D0
    res = id_sym / divisor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "1.- Error: Symbolic division failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if

    ! 2.- Division by small number
    divisor = 0.1D0
    expected_val = 10.0D0
    res = id_sym / divisor
    passed = abs(res%val - expected_val) < EPS
    if (.not. passed) then
        print*, "2.- Error: Division by small number failed", new_line('A'), &
                "Expected val:", expected_val, " Actual:", res%val
        return
    end if
end subroutine test_iden_4O4T_div