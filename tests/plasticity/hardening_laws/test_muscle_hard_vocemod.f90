! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_hard_vocemod
    implicit none

    logical :: passed

    call test_vocemod(passed)
    if (.not. passed) STOP 1

    call test_dvocemod(passed)
    if (.not. passed) STOP 2

    call test_ddvocemod(passed)
    if (.not. passed) STOP 3

    print*, "Passed!", passed
    STOP 0
end program test_muscle_hard_vocemod

! Law of muscle_hard_vocemod: sigma(ep) = Sy + k*ep + q*(1 - exp(-n*ep)). With Sy = 200, k = 500,
! q = 150 and n = 20: at ep = 0 only Sy, k + q*n and -q*n**2 remain; at ep = ln(2)/n,
! exp(-n*ep) = 1/2 checks the rate n. The expected values follow by hand.

subroutine test_vocemod(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_vocemod_dp
    implicit none

    real(real64), parameter :: EPS = 1e-10
    logical, intent(out) :: passed

    type(hard_vocemod_dp) :: voce
    real(real64) :: result, expected_result1

    voce = hard_vocemod_dp(sy=200D0, k=500D0, q=150D0, n=20D0)

    ! Initial yield stress: the Voce term is still zero
    expected_result1 = 200D0
    result = hard_stress(voce, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    ! Half of the saturation range q: 200 + 25*ln(2) + 75
    expected_result1 = 275D0 + 25D0*log(2D0)
    result = hard_stress(voce, log(2D0)/20D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return
end subroutine


subroutine test_dvocemod(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_vocemod_dp
    implicit none

    real(real64), parameter :: EPS = 1e-10
    logical, intent(out) :: passed

    type(hard_vocemod_dp) :: voce
    real(real64) :: result, expected_result1

    voce = hard_vocemod_dp(sy=200D0, k=500D0, q=150D0, n=20D0)

    ! Initial slope k + q*n
    expected_result1 = 3500D0
    result = hard_dstress_dep(voce, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    ! k + q*n/2
    expected_result1 = 2000D0
    result = hard_dstress_dep(voce, log(2D0)/20D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return
end subroutine


subroutine test_ddvocemod(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_vocemod_dp
    implicit none

    real(real64), parameter :: EPS = 1e-10
    logical, intent(out) :: passed

    type(hard_vocemod_dp) :: voce
    real(real64) :: result, expected_result1

    voce = hard_vocemod_dp(sy=200D0, k=500D0, q=150D0, n=20D0)

    ! -q*n**2: negative, the slope decreases towards k
    expected_result1 = -60000D0
    result = hard_ddstress_ddep(voce, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    ! -q*n**2/2
    expected_result1 = -30000D0
    result = hard_ddstress_ddep(voce, log(2D0)/20D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return
end subroutine
