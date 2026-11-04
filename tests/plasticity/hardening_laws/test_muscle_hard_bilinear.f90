! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
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
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_bilinear_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_bilinear_dp) :: hard

    real(real64) :: result
    real(real64) :: expected_result1

    hard = hard_bilinear_dp(K=1D0, y0=0D0)

    expected_result1 = 0D0
    result = hard_stress(hard, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = hard_stress(hard, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    hard = hard_bilinear_dp(K=1D0, y0=1D0)

    expected_result1 = 1D0
    result = hard_stress(hard, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = hard_stress(hard, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine



subroutine test_dbilinear(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_bilinear_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_bilinear_dp) :: hard

    real(real64) :: result
    real(real64) :: expected_result1

    hard = hard_bilinear_dp(K=1D0, y0=0D0)

    expected_result1 = 1D0
    result = hard_dstress_dep(hard, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = hard_dstress_dep(hard, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    hard = hard_bilinear_dp(K=2D0, y0=10D0)

    expected_result1 = 2.0D0
    result = hard_dstress_dep(hard, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = hard_dstress_dep(hard, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine


subroutine test_ddbilinear(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_bilinear_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_bilinear_dp) :: hard

    real(real64) :: result
    real(real64) :: expected_result1

    hard = hard_bilinear_dp(K=10D0, y0=10D0)

    expected_result1 = 0D0
    result = hard_ddstress_ddep(hard, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 0D0
    result = hard_ddstress_ddep(hard, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 0D0
    result = hard_ddstress_ddep(hard, 10D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************
end subroutine