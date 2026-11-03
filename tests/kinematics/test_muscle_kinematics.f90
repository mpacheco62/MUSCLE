! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_kinematics
    use, intrinsic :: iso_fortran_env
    implicit none

    logical :: passed

    print*, "========================================================"
    print*, "  RUNNING FULL KINEMATICS SUBSYSTEM TEST SUITE          "
    print*, "========================================================"

    print*, "1. Small Strain Kinematics..."
    call test_small_strain_kinematics(passed)
    if (.not. passed) stop 1

    print*, "2. Total Lagrangian Kinematics (Green-Lagrange)..."
    call test_total_lagrangian_kinematics(passed)
    if (.not. passed) stop 2

    print*, "3. Updated Lagrangian Kinematics (Euler-Almansi)..."
    call test_updated_lagrangian_kinematics(passed)
    if (.not. passed) stop 3

    print*, "4. Material Hencky Logarithmic Kinematics..."
    call test_material_logarithmic_kinematics(passed)
    if (.not. passed) stop 4

    print*, "5. Spatial Hencky Logarithmic Kinematics..."
    call test_spatial_logarithmic_kinematics(passed)
    if (.not. passed) stop 5

    print*, "6. Corotational Frame Kinematics (ANSYS Option B)..."
    call test_corotational_kinematics(passed)
    if (.not. passed) stop 6

    print*, "7. Stress Conversion Matrix Round-Trip..."
    call test_stress_conversions_roundtrip(passed)
    if (.not. passed) stop 7

    print*, "8. 4th-Order Tangent Modulus Transformations..."
    call test_tangent_modulus_transformations(passed)
    if (.not. passed) stop 8

    print*, "9. Inverse quantities of finite-strain kinematics (F^-1, b^-1)..."
    call test_finite_inverse_quantities(passed)
    if (.not. passed) stop 9

    print*, "========================================================"
    print*, "  ALL 9 KINEMATICS TESTS PASSED SUCCESSFULLY!          "
    print*, "========================================================"
    stop 0
end program test_muscle_kinematics


! =========================================================================
! 1. Small Strain Kinematics
! =========================================================================
subroutine test_small_strain_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    implicit none
    logical, intent(out) :: passed

    type(kin_small_strain_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2O)                :: F
    type(ten_3D2Osym)             :: S_material, sigma_spatial, strain_expected, strain_actual
    type(ten_3D4O2sym)            :: C_material, c_spatial

    passed = .FALSE.

    call F%init(xx=1.002D0, xy=0.004D0, xz=0.0D0, yx=0.0D0, yy=0.997D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.001D0)
    call kin_update(kin, F, st, status)

    call strain_expected%init(xx=0.002D0, yy=-0.003D0, zz=0.001D0, xy=0.002D0, yz=0.0D0, xz=0.0D0)
    strain_actual = st%strain
    if (.not. (strain_actual .approx. strain_expected)) return

    if (abs(st%J - 1.0D0) > 1.0D-10 .or. kin_is_finite_strain(kin)) return

    call S_material%init(xx=100.0D0, yy=-50.0D0, zz=20.0D0, xy=15.0D0, yz=-5.0D0, xz=8.0D0)
    sigma_spatial = kin_to_cauchy(kin, st, S_material)
    if (.not. (sigma_spatial .approx. S_material)) return

    passed = .TRUE.
end subroutine test_small_strain_kinematics


! =========================================================================
! 2. Total Lagrangian Kinematics
! =========================================================================
subroutine test_total_lagrangian_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    implicit none
    logical, intent(out) :: passed

    type(kin_total_lagrangian_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2O)                     :: F
    type(ten_3D2Osym)                  :: C_expected, C_actual, E_expected, E_actual
    type(ten_3D2Osym)                  :: S_pk2, sigma_expected, sigma_actual

    passed = .FALSE.

    call F%init(xx=1.2D0, xy=0.2D0, xz=0.0D0, yx=0.0D0, yy=0.9D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin_update(kin, F, st, status)

    ! C = F^T * F
    call C_expected%init(xx=1.44D0, yy=0.85D0, zz=1.0D0, xy=0.24D0, yz=0.0D0, xz=0.0D0)
    C_actual = st%C
    if (.not. (C_actual .approx. C_expected)) return

    ! E_GL = 0.5 * (C - I)
    call E_expected%init(xx=0.22D0, yy=-0.075D0, zz=0.0D0, xy=0.12D0, yz=0.0D0, xz=0.0D0)
    E_actual = st%strain
    if (.not. (E_actual .approx. E_expected)) return

    ! Push-Forward PK2 -> Cauchy
    call S_pk2%init(xx=100.0D0, yy=50.0D0, zz=-20.0D0, xy=10.0D0, yz=0.0D0, xz=0.0D0)
    call sigma_expected%init(xx=139.6296296296296D0, yy=37.5D0, zz=-18.5185185185185D0, xy=18.3333333333333D0, yz=0.0D0, xz=0.0D0)

    sigma_actual = kin_to_cauchy(kin, st, S_pk2)
    if (.not. (sigma_actual .approx. sigma_expected)) return

    passed = .TRUE.
end subroutine test_total_lagrangian_kinematics


! =========================================================================
! 3. Updated Lagrangian Kinematics
! =========================================================================
subroutine test_updated_lagrangian_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    implicit none
    logical, intent(out) :: passed

    type(kin_updated_lagrangian_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2O)                      :: F
    type(ten_3D2Osym)                   :: e_expected, e_actual
    type(ten_3D2Osym)                   :: sigma_spatial, S_pk2, sigma_expected

    passed = .FALSE.

    call F%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, yx=0.0D0, yy=1.0D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin_update(kin, F, st, status)

    ! e_EA = 0.5 * (I - b^-1) = diag(0.375, 0, 0)
    call e_expected%init(xx=0.375D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    e_actual = st%strain
    if (.not. (e_actual .approx. e_expected)) return

    ! to_Cauchy in UL is identity
    call sigma_spatial%init(xx=100.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    if (.not. (kin_to_cauchy(kin, st, sigma_spatial) .approx. sigma_spatial)) return

    ! Geometric push_forward_stress PK2 -> Cauchy
    call S_pk2%init(xx=50.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call sigma_expected%init(xx=100.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    if (.not. (kin_pk2_to_cauchy(st, S_pk2) .approx. sigma_expected)) return

    passed = .TRUE.
end subroutine test_updated_lagrangian_kinematics


! =========================================================================
! 4. Material Logarithmic Kinematics
! =========================================================================
subroutine test_material_logarithmic_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    implicit none
    logical, intent(out) :: passed

    type(kin_material_logarithmic_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2O)                        :: F
    type(ten_3D2Osym)                     :: H_expected, H_actual, tau_R, sigma_expected, sigma_actual

    passed = .FALSE.

    ! F = diag(2.0, 0.5, 1.0)
    call F%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, yx=0.0D0, yy=0.5D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin_update(kin, F, st, status)

    ! H = ln(U) = diag(ln(2.0), ln(0.5), 0.0)
    call H_expected%init(xx=0.6931471805599453D0, yy=-0.6931471805599453D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    H_actual = st%strain
    if (.not. (H_actual .approx. H_expected)) return

    ! to_Cauchy converts Rotated Kirchhoff stress tau_R to Cauchy: sigma = (1/J) * R * tau_R * R^T
    ! For diagonal F, R = I, J = 1.0 => sigma = tau_R
    call tau_R%init(xx=100.0D0, yy=50.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call sigma_expected%init(xx=100.0D0, yy=50.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    sigma_actual = kin_to_cauchy(kin, st, tau_R)
    if (.not. (sigma_actual .approx. sigma_expected)) return

    passed = .TRUE.
end subroutine test_material_logarithmic_kinematics


! =========================================================================
! 5. Spatial Logarithmic Kinematics
! =========================================================================
subroutine test_spatial_logarithmic_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    implicit none
    logical, intent(out) :: passed

    type(kin_spatial_logarithmic_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2O)                       :: F
    type(ten_3D2Osym)                    :: h_expected, h_actual, tau_spatial, sigma_expected, sigma_actual

    passed = .FALSE.

    ! F = diag(2.0, 0.5, 1.0) => J = 1.0
    call F%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, yx=0.0D0, yy=0.5D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin_update(kin, F, st, status)

    ! h = ln(V) = diag(ln(2.0), ln(0.5), 0.0)
    call h_expected%init(xx=0.6931471805599453D0, yy=-0.6931471805599453D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    h_actual = st%strain
    if (.not. (h_actual .approx. h_expected)) return

    ! to_Cauchy converts Kirchhoff stress tau -> Cauchy: sigma = tau / J
    call tau_spatial%init(xx=100.0D0, yy=50.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call sigma_expected%init(xx=100.0D0, yy=50.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    sigma_actual = kin_to_cauchy(kin, st, tau_spatial)
    if (.not. (sigma_actual .approx. sigma_expected)) return

    passed = .TRUE.
end subroutine test_spatial_logarithmic_kinematics


! =========================================================================
! 6. Corotational Kinematics (ANSYS Option B / LS-DYNA)
! =========================================================================
subroutine test_corotational_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    implicit none
    logical, intent(out) :: passed

    type(kin_corotational_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2Osym)             :: strain_unrotated, sigma_unrotated, sigma_global, sigma_expected
    type(ten_3D2O)                :: R_90deg

    passed = .FALSE.

    ! 90 degree rotation around Z axis:
    ! R = [ 0, -1, 0 ; 1, 0, 0 ; 0, 0, 1 ]
    call R_90deg%init(xx=0.0D0, xy=-1.0D0, xz=0.0D0, &
                      yx=1.0D0, yy=0.0D0,  yz=0.0D0, &
                      zx=0.0D0, zy=0.0D0,  zz=1.0D0)

    call strain_unrotated%init(xx=0.1D0, yy=-0.05D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    call kin_corot_update(kin, strain_unrotated, st, status, R=R_90deg)

    ! Unrotated stress: sigma_hat = diag(100, 20, 0)
    call sigma_unrotated%init(xx=100.0D0, yy=20.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    ! Rotated global Cauchy stress: sigma = R * sigma_hat * R^T
    ! xx_global = yy_unrotated = 20, yy_global = xx_unrotated = 100
    call sigma_expected%init(xx=20.0D0, yy=100.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    sigma_global = kin_to_cauchy(kin, st, sigma_unrotated)

    if (.not. (sigma_global .approx. sigma_expected)) then
        print*, "FAIL: Corotational stress rotation R * S * R^T failed", new_line('A'), &
                "Expected:", sigma_expected%vals, new_line('A'), &
                "Actual  :", sigma_global%vals
        return
    end if

    passed = .TRUE.
end subroutine test_corotational_kinematics


! =========================================================================
! 7. Stress Conversions Matrix Round-Trip
! =========================================================================
subroutine test_stress_conversions_roundtrip(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    implicit none
    logical, intent(out) :: passed

    type(kin_total_lagrangian_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2O)                     :: F
    type(ten_3D2Osym)                  :: S_orig, S_recovered, sigma, tau
    type(ten_3D2O)                     :: P_pk1

    passed = .FALSE.

    call F%init(xx=1.1D0, xy=0.2D0, xz=0.0D0, yx=0.1D0, yy=0.95D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin_update(kin, F, st, status)

    call S_orig%init(xx=120.0D0, yy=40.0D0, zz=-10.0D0, xy=15.0D0, yz=0.0D0, xz=0.0D0)

    ! 1. PK2 -> Cauchy -> PK2
    sigma = kin_pk2_to_cauchy(st, S_orig)
    S_recovered = kin_cauchy_to_pk2(st, sigma)
    if (.not. (S_recovered .approx. S_orig)) return

    ! 2. PK2 -> Kirchhoff -> PK2
    tau = kin_pk2_to_kirchhoff(st, S_orig)
    S_recovered = kin_kirchhoff_to_pk2(st, tau)
    if (.not. (S_recovered .approx. S_orig)) return

    ! 3. PK2 -> 1st PK (P) -> PK2
    P_pk1 = kin_pk2_to_pk1(st, S_orig)
    S_recovered = kin_pk1_to_pk2(st, P_pk1)
    if (.not. (S_recovered .approx. S_orig)) return

    passed = .TRUE.
end subroutine test_stress_conversions_roundtrip


! =========================================================================
! 8. 4th-Order Tangent Modulus Transformations
! =========================================================================
subroutine test_tangent_modulus_transformations(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    implicit none
    logical, intent(out) :: passed

    type(kin_total_lagrangian_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2O)                     :: F
    type(ten_3D4O2sym)                 :: C_mat, c_spat, C_recovered

    passed = .FALSE.

    ! F = diag(2.0, 1.0, 1.0) => J = 2.0
    call F%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, yx=0.0D0, yy=1.0D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin_update(kin, F, st, status)

    C_mat%vals = 0.0D0
    C_mat%vals(1,1) = 200.0D0
    C_mat%vals(2,2) = 200.0D0

    ! Push-Forward: c_1111 = (1/J) * F_11^4 * C_1111 = (1/2) * 16 * 200 = 1600.0
    c_spat = kin_push_forward_tangent(st, C_mat)

    if (abs(c_spat%vals(1,1) - 1600.0D0) > 1.0D-6) then
        print*, "FAIL: 4th-order push-forward tangent incorrect", new_line('A'), &
                "Expected c_1111 = 1600.0, got:", c_spat%vals(1,1)
        return
    end if

    ! Pull-Back: C_recovered = J * F_11^-4 * c_1111 = 2.0 * (0.5)^4 * 1600.0 = 200.0
    C_recovered = kin_pull_back_tangent(st, c_spat)

    if (.not. (C_recovered .approx. C_mat)) then
        print*, "FAIL: 4th-order tangent round-trip pull-back failed"
        return
    end if

    passed = .TRUE.
end subroutine test_tangent_modulus_transformations


! =========================================================================
! 9. Inverse quantities of finite-strain kinematics
! =========================================================================
subroutine test_finite_inverse_quantities(passed)
    !! F^-1 and b^-1 = F^-T F^-1 computed by `kin_update` (through the tensor inverse, M33INV):
    !!
    !! 1. General F with an exactly known inverse: F = L U with unit lower/upper triangular factors
    !!    (L21 = 1/2, L31 = -1/4, L32 = 3/4, U12 = 1/4, U13 = 1/2, U23 = -1/2), so det F = 1 and
    !!    F^-1 = U^-1 L^-1 (closed form). All entries are short dyadic fractions, so F^-1 and b^-1
    !!    are exact in binary: compared exactly (tol = 0).
    !! 2. Ill-conditioned F = Q diag(1e5, 1, 1e-5) Q^T (Q: rotation of 0.7 rad about (1,2,2)/3),
    !!    cond ~ 1e10. Reference: inverse of the same double-precision F in quad precision. A
    !!    backward-stable inverse has a relative error O(kappa eps), so tol = 10 kappa_inf eps.
    !!    The former hand-written adjugate in the OO update_F had no LU fallback and missed this bound
    !!    (error ~1e-4 relative at cond 1e10).
    use, intrinsic :: iso_fortran_env, only : real64, real128
    use muscle_tensors
    use muscle_kin_laws
    use muscle_kin_finite
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(kin_total_lagrangian_dp) :: kin
    type(kin_state_dp) :: st
    integer :: status
    type(ten_3D2O) :: F, F_inv_exact
    type(ten_3D2Osym) :: b_inv_exact
    real(real64) :: L(3,3), U(3,3), L_inv(3,3), U_inv(3,3), g(3,3), q(3,3), n(3), d(3, 3), f3(3,3)
    real(real64) :: ct, sth, kappa, x(3,3)
    real(real128) :: xr(3,3)
    real(real64), parameter :: THETA = 0.7D0

    passed = .true.

    ! --- 1. General F with an exact inverse ---
    L     = reshape((/ 1D0, 0.5D0, -0.25D0,   0D0, 1D0, 0.75D0,   0D0, 0D0, 1D0 /), (/3,3/))
    U     = reshape((/ 1D0, 0D0, 0D0,   0.25D0, 1D0, 0D0,   0.5D0, -0.5D0, 1D0 /), (/3,3/))
    L_inv = reshape((/ 1D0, -0.5D0, 0.5D0*0.75D0 + 0.25D0,   0D0, 1D0, -0.75D0,   0D0, 0D0, 1D0 /), (/3,3/))
    U_inv = reshape((/ 1D0, 0D0, 0D0,   -0.25D0, 1D0, 0D0,   0.25D0*(-0.5D0) - 0.5D0, 0.5D0, 1D0 /), (/3,3/))
    g = matmul(U_inv, L_inv)
    call F%init(reshape(matmul(L, U), (/9/)))
    call F_inv_exact%init(reshape(g, (/9/)))
    x = matmul(transpose(g), g)                                   ! b^-1 = F^-T F^-1
    call b_inv_exact%init(xx=x(1,1), yy=x(2,2), zz=x(3,3), xy=x(1,2), yz=x(2,3), xz=x(1,3))

    call kin_update(kin, F, st, status)
    call check_approx(passed, "kin_update: F^-1 (F = L U, exact)", st%F_inv, F_inv_exact, tol=0D0)
    call check_approx(passed, "kin_update: b^-1 = F^-T F^-1 (exact)", st%b_inv, b_inv_exact, tol=0D0)
    call check_approx(passed, "kin_update: J = det F = 1", st%J, 1D0, tol=0D0)

    ! --- 2. Ill-conditioned F ---
    n = (/1D0, 2D0, 2D0/) / 3D0
    ct = cos(THETA); sth = sin(THETA)
    q(1,:) = (/ ct + n(1)**2*(1-ct),        n(1)*n(2)*(1-ct) - n(3)*sth, n(1)*n(3)*(1-ct) + n(2)*sth /)
    q(2,:) = (/ n(2)*n(1)*(1-ct) + n(3)*sth, ct + n(2)**2*(1-ct),        n(2)*n(3)*(1-ct) - n(1)*sth /)
    q(3,:) = (/ n(3)*n(1)*(1-ct) - n(2)*sth, n(3)*n(2)*(1-ct) + n(1)*sth, ct + n(3)**2*(1-ct)        /)
    d = 0D0
    d(1,1) = 1D5; d(2,2) = 1D0; d(3,3) = 1D-5
    f3 = matmul(q, matmul(d, transpose(q)))
    call F%init(reshape(f3, (/9/)))
    call quad_inverse_3(f3, xr)
    x = real(xr, real64)
    call F_inv_exact%init(reshape(x, (/9/)))
    kappa = maxval(sum(abs(f3), dim=2)) * maxval(sum(abs(x), dim=2))

    call kin_update(kin, F, st, status)
    call check_approx(passed, "kin_update: F^-1 of an ill-conditioned F (cond ~1e10)", st%F_inv, &
                      F_inv_exact, tol=10D0*kappa*epsilon(1D0))

contains

    subroutine quad_inverse_3(a, ainv)
        !! Gauss-Jordan with partial pivoting in quad precision (reference inverse).
        real(real64), intent(in) :: a(3,3)
            !! Matrix (double precision).
        real(real128), intent(out) :: ainv(3,3)
            !! Its inverse, computed in quad precision.
        real(real128) :: w(3,6), tmp(6), fac
        integer :: ii, kk, p
        w = 0; w(:,1:3) = real(a, real128)
        do ii = 1, 3
            w(ii,3+ii) = 1
        end do
        do kk = 1, 3
            p = kk - 1 + maxloc(abs(w(kk:3,kk)), 1)
            tmp = w(kk,:); w(kk,:) = w(p,:); w(p,:) = tmp
            w(kk,:) = w(kk,:) / w(kk,kk)
            do ii = 1, 3
                if (ii /= kk) then
                    fac = w(ii,kk); w(ii,:) = w(ii,:) - fac*w(kk,:)
                end if
            end do
        end do
        ainv = w(:,4:6)
    end subroutine quad_inverse_3
end subroutine test_finite_inverse_quantities
