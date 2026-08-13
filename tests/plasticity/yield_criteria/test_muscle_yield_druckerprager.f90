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

    call test_dp_tension_uniaxial_comp_abaqus(passed)
    if (.not. passed) stop 5

    call test_dp_compression_uniaxial_comp_abaqus(passed)
    if (.not. passed) stop 6

    call test_dp_biaxial_comp_abaqus(passed)
    if (.not. passed) stop 7

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

subroutine test_dp_tension_uniaxial_comp_abaqus(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_hard_bilinear
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-3

    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------
    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    expected = bilinear%stress(0.09344D0)
    expected_tan = -2.581929D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=393.4D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tension test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (TENSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tangent test (TENSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (TENSIONDEF)"
    end if

    ! ---------------------------------------------------
    ! --------  Compression type  -----------------------
    ! ---------------------------------------------------
    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    expected = bilinear%stress(0.06688D0)
    expected_tan = -2.58221D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=260.8D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tension test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (COMPRESSIONDEF)"
    end if

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tangent test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (COMPRESSIONDEF)"
    end if

    ! ---------------------------------------------------
    ! --------  Shear type  -----------------------------
    ! ---------------------------------------------------
    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    expected = bilinear%stress(0.08032D0)
    expected_tan = -2.58221D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=299.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tension test (SHEARDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Tension Test (SHEARDEF)"
    end if

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Tangent test (SHEARDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Tangent Test (SHEARDEF)"
    end if
end subroutine test_dp_tension_uniaxial_comp_abaqus

subroutine test_dp_compression_uniaxial_comp_abaqus(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_hard_bilinear
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-3

    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)
    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.0112D0)
    expected_tan = -1.5180952

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=-437.7D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (TENSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression Tangent test (TENSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (TENSIONDEF)"
    end if

    ! ---------------------------------------------------
    ! --------  Compression type  -----------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.00858D0)
    expected_tan = -1.518584D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=-308.58D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (COMPRESSIONDEF)"
    end if

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression Tangent test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (COMPRESSIONDEF)"
    end if

    ! ---------------------------------------------------
    ! --------  Shear type  -----------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.01013D0)
    expected_tan = -1.5198555

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=-342.9D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression test (SHEARDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Uniaxial Compression Test (SHEARDEF)"
    end if

    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan) < EPS*1.5 ! Falta de decimales
    if (.not. passed) then
        print *, "FAIL: Uniaxial Compression Tangent test (SHEARDEF)"
        print *, "  Expected:", expected_tan
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan)
    else
        print *, "PASS: Uniaxial Compression Tangent Test (SHEARDEF)"
    end if
end subroutine test_dp_compression_uniaxial_comp_abaqus


subroutine test_dp_biaxial_comp_abaqus(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_druckerprager
    use muscle_hard_bilinear
    implicit none
    logical, intent(out) :: passed

    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: stress, tangent
    type(Bilinear_hardening) :: bilinear
    real(real64) :: res, expected, expected_tan1, expected_tan2, beta, K_ratio
    real(real64), parameter :: EPS = 1.0D-4

    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0   ! (0.778 <= K <= 1.0)

    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.03497D0)
    expected_tan1 = -0.471306D0
    expected_tan2 = 0.55486D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)
    call stress%init(xx=-994.63D0, yy=-1711.35D0, zz=-901.47D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial test (TENSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (TENSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan1)/expected_tan1 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent1 test (TENSIONDEF)"
        print *, "  Expected:", expected_tan1
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan1)
    else
        print *, "PASS: Biaxial Tangent1 Test (TENSIONDEF)"
    end if

    res = tangent%xx()/tangent%zz()
    passed = abs(res - expected_tan2)/expected_tan2 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent2 test (TENSIONDEF)"
        print *, "  Expected:", expected_tan2
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan2)
    else
        print *, "PASS: Biaxial Tangent2 Test (TENSIONDEF)"
    end if

    ! ---------------------------------------------------
    ! --------  Compression type  -----------------------
    ! ---------------------------------------------------
    expected = bilinear%stress(0.02531D0)
    expected_tan1 = -0.472822D0
    expected_tan2 = 0.55775577D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_COMPRESSION)
    call stress%init(xx=-1048.47D0, yy=-1647.97D0, zz=-970.69D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (COMPRESSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan1)/expected_tan1 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent1 test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan1
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan1)
    else
        print *, "PASS: Biaxial Tangent1 Test (COMPRESSIONDEF)"
    end if

    res = tangent%xx()/tangent%zz()
    passed = abs(res - expected_tan2)/expected_tan2 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent2 test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan2
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan2)
    else
        print *, "PASS: Biaxial Tangent2 Test (COMPRESSIONDEF)"
    end if

    ! ! ---------------------------------------------------
    ! ! --------  Shear type  -----------------------------
    ! ! ---------------------------------------------------
    expected = bilinear%stress(0.03031D0)
    expected_tan1 = -0.4724789D0
    expected_tan2 = 0.55682D0

    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_SHEAR)
    call stress%init(xx=-1033.74D0, yy=-1665.32D0, zz=-951.75D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    res = dp%stress_eq(stress)

    passed = abs(res - expected)/expected < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial test (COMPRESSIONDEF)"
        print *, "  Expected:", expected
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected)
    else
        print *, "PASS: Biaxial Test (COMPRESSIONDEF)"
    end if
    
    tangent = dp%dstressEq_dstress(stress)
    res = tangent%xx()/tangent%yy()
    passed = abs(res - expected_tan1)/expected_tan1 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent1 test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan1
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan1)
    else
        print *, "PASS: Biaxial Tangent1 Test (COMPRESSIONDEF)"
    end if

    res = tangent%xx()/tangent%zz()
    passed = abs(res - expected_tan2)/expected_tan2 < EPS
    if (.not. passed) then
        print *, "FAIL: Biaxial Tangent2 test (COMPRESSIONDEF)"
        print *, "  Expected:", expected_tan2
        print *, "  Got     :", res
        print *, "  Diff    :", (res - expected_tan2)
    else
        print *, "PASS: Biaxial Tangent2 Test (COMPRESSIONDEF)"
    end if

end subroutine test_dp_biaxial_comp_abaqus