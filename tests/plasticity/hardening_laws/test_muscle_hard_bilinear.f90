program test_muscle_hard_bilinear
    implicit none
    
    logical :: passed

    call test_bilinear(passed)
    if (.not. passed) STOP 1

    call test_dbilinear(passed)
    if (.not. passed) STOP 2

    call test_ddbilinear(passed)
    if (.not. passed) STOP 3

    print*, "Passed!", passed
    STOP 0
end program test_muscle_hard_bilinear

subroutine test_bilinear(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_bilinear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Bilinear_hardening) :: hard

    real(real64) :: result
    real(real64) :: expected_result1

    hard = Bilinear_hardening(K=1D0, y0=0D0)

    expected_result1 = 0D0
    result = hard%stress(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = hard%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    hard = Bilinear_hardening(K=1D0, y0=1D0)

    expected_result1 = 1D0
    result = hard%stress(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = hard%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine



subroutine test_dbilinear(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_bilinear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Bilinear_hardening) :: hard

    real(real64) :: result
    real(real64) :: expected_result1

    hard = Bilinear_hardening(K=1D0, y0=0D0)

    expected_result1 = 1D0
    result = hard%dstress_dep(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = hard%dstress_dep(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    hard = Bilinear_hardening(K=2D0, y0=10D0)

    expected_result1 = 2.0D0
    result = hard%dstress_dep(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = hard%dstress_dep(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine


subroutine test_ddbilinear(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_bilinear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Bilinear_hardening) :: hard

    real(real64) :: result
    real(real64) :: expected_result1

    hard = Bilinear_hardening(K=10D0, y0=10D0)

    expected_result1 = 0D0
    result = hard%ddstress_ddep(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 0D0
    result = hard%ddstress_ddep(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 0D0
    result = hard%ddstress_ddep(10D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************
end subroutine