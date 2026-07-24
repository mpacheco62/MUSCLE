! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_dp
    implicit none
    logical :: passed

    call test_dp_tension_uniaxial(passed)
    if (.not. passed) stop 1

    call test_dp_compression_uniaxial(passed)
    if (.not. passed) stop 2

    call test_dp_shear_calibration(passed)
    if (.not. passed) stop 3

    call test_dp_pure_hydrostatic(passed)
    if (.not. passed) stop 4

    stop 0
end program test_muscle_yield_dp


subroutine test_dp_tension_uniaxial(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, expected, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-7

    beta = 20.0D0
    K_ratio = 0.8D0   ! (0.778 <= K <= 1.0)
    expected = 150.0D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)

    call stress%init(xx=expected, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tension test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", res - expected
    else
        print *, "PASS: Uniaxial Tension Test"
    end if
end subroutine test_dp_tension_uniaxial


subroutine test_dp_compression_uniaxial(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, applied_stress, expected, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-7

    beta = 20.0D0
    K_ratio = 0.8D0
    applied_stress = -150.0D0
    expected = 150.0D0 

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)

    call stress%init(xx=applied_stress, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", res - expected
    else
        print *, "PASS: Uniaxial Compression Test"
    end if
end subroutine test_dp_compression_uniaxial


subroutine test_dp_shear_calibration(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, expected, beta, K_ratio, tau, q_val
    real(real64), parameter :: EPS = 1.0D-7

    beta = 10.0D0
    K_ratio = 0.8D0
    tau = 80.0D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)

    call stress%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=tau, yz=0.0D0, xz=0.0D0)

    q_val = sqrt(3.0D0) * tau
    expected = 0.5D0 * q_val * (1.0D0 + 1.0D0 / K_ratio)

    res = dp%stress_eq(stress)

    passed = abs(res - expected) < EPS
    if (.not. passed) then
        print *, "FAIL: Pure Shear test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", res - expected
    else
        print *, "PASS: Pure Shear Test"
    end if
end subroutine test_dp_shear_calibration


subroutine test_dp_pure_hydrostatic(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress
    real(real64) :: res, expected, beta, K_ratio, p_hydro, tanbeta, PI
    real(real64), parameter :: EPS = 1.0D-7

    PI = acos(-1.0D0)
    beta = 30.0D0
    K_ratio = 1.0D0
    p_hydro = 100.0D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)

    call stress%init(xx=p_hydro, yy=p_hydro, zz=p_hydro, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    tanbeta = tan(beta * PI / 180.0D0)
    expected = (tanbeta * (1.0D0 + tanbeta / 3.0D0)**(-1)) * p_hydro

    res = dp%stress_eq(stress)

    passed = abs(res - expected) < EPS
    if (.not. passed) then
        print *, "FAIL: Pure Hydrostatic test"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", res - expected
    else
        print *, "PASS: Pure Hydrostatic Test"
    end if
end subroutine test_dp_pure_hydrostatic