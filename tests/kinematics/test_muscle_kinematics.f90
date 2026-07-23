! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_kinematics
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kinematics
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

    print*, "========================================================"
    print*, "  ALL 8 KINEMATICS TESTS PASSED SUCCESSFULLY!          "
    print*, "========================================================"
    stop 0
end program test_muscle_kinematics


! =========================================================================
! 1. Small Strain Kinematics
! =========================================================================
subroutine test_small_strain_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kinematics
    implicit none
    logical, intent(out) :: passed

    type(Small_strain_kinematics) :: kin
    type(ten_3D2O)                :: F
    type(ten_3D2Osym)             :: S_material, sigma_spatial, strain_expected, strain_actual
    type(ten_3D4O2sym)            :: C_material, c_spatial

    passed = .FALSE.

    call F%init(xx=1.002D0, xy=0.004D0, xz=0.0D0, yx=0.0D0, yy=0.997D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.001D0)
    call kin%update(F)

    call strain_expected%init(xx=0.002D0, yy=-0.003D0, zz=0.001D0, xy=0.002D0, yz=0.0D0, xz=0.0D0)
    strain_actual = kin%get_strain()
    if (.not. (strain_actual .approx. strain_expected)) return

    if (abs(kin%jacobian() - 1.0D0) > 1.0D-10 .or. kin%is_finite_strain()) return

    call S_material%init(xx=100.0D0, yy=-50.0D0, zz=20.0D0, xy=15.0D0, yz=-5.0D0, xz=8.0D0)
    sigma_spatial = kin%to_Cauchy(S_material)
    if (.not. (sigma_spatial .approx. S_material)) return

    passed = .TRUE.
end subroutine test_small_strain_kinematics


! =========================================================================
! 2. Total Lagrangian Kinematics
! =========================================================================
subroutine test_total_lagrangian_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kinematics
    implicit none
    logical, intent(out) :: passed

    type(Total_lagrangian_kinematics) :: kin
    type(ten_3D2O)                     :: F
    type(ten_3D2Osym)                  :: C_expected, C_actual, E_expected, E_actual
    type(ten_3D2Osym)                  :: S_pk2, sigma_expected, sigma_actual

    passed = .FALSE.

    call F%init(xx=1.2D0, xy=0.2D0, xz=0.0D0, yx=0.0D0, yy=0.9D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin%update(F)

    ! C = F^T * F
    call C_expected%init(xx=1.44D0, yy=0.85D0, zz=1.0D0, xy=0.24D0, yz=0.0D0, xz=0.0D0)
    C_actual = kin%get_RightCauchyGreen()
    if (.not. (C_actual .approx. C_expected)) return

    ! E_GL = 0.5 * (C - I)
    call E_expected%init(xx=0.22D0, yy=-0.075D0, zz=0.0D0, xy=0.12D0, yz=0.0D0, xz=0.0D0)
    E_actual = kin%get_strain()
    if (.not. (E_actual .approx. E_expected)) return

    ! Push-Forward PK2 -> Cauchy
    call S_pk2%init(xx=100.0D0, yy=50.0D0, zz=-20.0D0, xy=10.0D0, yz=0.0D0, xz=0.0D0)
    call sigma_expected%init(xx=139.6296296296296D0, yy=37.5D0, zz=-18.5185185185185D0, xy=18.3333333333333D0, yz=0.0D0, xz=0.0D0)

    sigma_actual = kin%to_Cauchy(S_pk2)
    if (.not. (sigma_actual .approx. sigma_expected)) return

    passed = .TRUE.
end subroutine test_total_lagrangian_kinematics


! =========================================================================
! 3. Updated Lagrangian Kinematics
! =========================================================================
subroutine test_updated_lagrangian_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kinematics
    implicit none
    logical, intent(out) :: passed

    type(Updated_lagrangian_kinematics) :: kin
    type(ten_3D2O)                      :: F
    type(ten_3D2Osym)                   :: e_expected, e_actual
    type(ten_3D2Osym)                   :: sigma_spatial, S_pk2, sigma_expected

    passed = .FALSE.

    call F%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, yx=0.0D0, yy=1.0D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin%update(F)

    ! e_EA = 0.5 * (I - b^-1) = diag(0.375, 0, 0)
    call e_expected%init(xx=0.375D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    e_actual = kin%get_strain()
    if (.not. (e_actual .approx. e_expected)) return

    ! to_Cauchy in UL is identity
    call sigma_spatial%init(xx=100.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    if (.not. (kin%to_Cauchy(sigma_spatial) .approx. sigma_spatial)) return

    ! Geometric push_forward_stress PK2 -> Cauchy
    call S_pk2%init(xx=50.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call sigma_expected%init(xx=100.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    if (.not. (kin%push_forward_stress(S_pk2) .approx. sigma_expected)) return

    passed = .TRUE.
end subroutine test_updated_lagrangian_kinematics


! =========================================================================
! 4. Material Logarithmic Kinematics
! =========================================================================
subroutine test_material_logarithmic_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kinematics
    implicit none
    logical, intent(out) :: passed

    type(Material_logarithmic_kinematics) :: kin
    type(ten_3D2O)                        :: F
    type(ten_3D2Osym)                     :: H_expected, H_actual, tau_R, sigma_expected, sigma_actual

    passed = .FALSE.

    ! F = diag(2.0, 0.5, 1.0)
    call F%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, yx=0.0D0, yy=0.5D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin%update(F)

    ! H = ln(U) = diag(ln(2.0), ln(0.5), 0.0)
    call H_expected%init(xx=0.6931471805599453D0, yy=-0.6931471805599453D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    H_actual = kin%get_strain()
    if (.not. (H_actual .approx. H_expected)) return

    ! to_Cauchy converts Rotated Kirchhoff stress tau_R to Cauchy: sigma = (1/J) * R * tau_R * R^T
    ! For diagonal F, R = I, J = 1.0 => sigma = tau_R
    call tau_R%init(xx=100.0D0, yy=50.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call sigma_expected%init(xx=100.0D0, yy=50.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    sigma_actual = kin%to_Cauchy(tau_R)
    if (.not. (sigma_actual .approx. sigma_expected)) return

    passed = .TRUE.
end subroutine test_material_logarithmic_kinematics


! =========================================================================
! 5. Spatial Logarithmic Kinematics
! =========================================================================
subroutine test_spatial_logarithmic_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kinematics
    implicit none
    logical, intent(out) :: passed

    type(Spatial_logarithmic_kinematics) :: kin
    type(ten_3D2O)                       :: F
    type(ten_3D2Osym)                    :: h_expected, h_actual, tau_spatial, sigma_expected, sigma_actual

    passed = .FALSE.

    ! F = diag(2.0, 0.5, 1.0) => J = 1.0
    call F%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, yx=0.0D0, yy=0.5D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin%update(F)

    ! h = ln(V) = diag(ln(2.0), ln(0.5), 0.0)
    call h_expected%init(xx=0.6931471805599453D0, yy=-0.6931471805599453D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    h_actual = kin%get_strain()
    if (.not. (h_actual .approx. h_expected)) return

    ! to_Cauchy converts Kirchhoff stress tau -> Cauchy: sigma = tau / J
    call tau_spatial%init(xx=100.0D0, yy=50.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    call sigma_expected%init(xx=100.0D0, yy=50.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    sigma_actual = kin%to_Cauchy(tau_spatial)
    if (.not. (sigma_actual .approx. sigma_expected)) return

    passed = .TRUE.
end subroutine test_spatial_logarithmic_kinematics


! =========================================================================
! 6. Corotational Kinematics (ANSYS Option B / LS-DYNA)
! =========================================================================
subroutine test_corotational_kinematics(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kinematics
    implicit none
    logical, intent(out) :: passed

    type(Corotational_kinematics) :: kin
    type(ten_3D2Osym)             :: strain_unrotated, sigma_unrotated, sigma_global, sigma_expected
    type(ten_3D2O)                :: R_90deg

    passed = .FALSE.

    ! 90 degree rotation around Z axis:
    ! R = [ 0, -1, 0 ; 1, 0, 0 ; 0, 0, 1 ]
    call R_90deg%init(xx=0.0D0, xy=-1.0D0, xz=0.0D0, &
                      yx=1.0D0, yy=0.0D0,  yz=0.0D0, &
                      zx=0.0D0, zy=0.0D0,  zz=1.0D0)

    call strain_unrotated%init(xx=0.1D0, yy=-0.05D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    call kin%update(strain = strain_unrotated, R = R_90deg)

    ! Unrotated stress: sigma_hat = diag(100, 20, 0)
    call sigma_unrotated%init(xx=100.0D0, yy=20.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    ! Rotated global Cauchy stress: sigma = R * sigma_hat * R^T
    ! xx_global = yy_unrotated = 20, yy_global = xx_unrotated = 100
    call sigma_expected%init(xx=20.0D0, yy=100.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    sigma_global = kin%to_Cauchy(sigma_unrotated)

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
    use muscle_kinematics
    implicit none
    logical, intent(out) :: passed

    type(Total_lagrangian_kinematics) :: kin
    type(ten_3D2O)                     :: F
    type(ten_3D2Osym)                  :: S_orig, S_recovered, sigma, tau
    type(ten_3D2O)                     :: P_pk1

    passed = .FALSE.

    call F%init(xx=1.1D0, xy=0.2D0, xz=0.0D0, yx=0.1D0, yy=0.95D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin%update(F)

    call S_orig%init(xx=120.0D0, yy=40.0D0, zz=-10.0D0, xy=15.0D0, yz=0.0D0, xz=0.0D0)

    ! 1. PK2 -> Cauchy -> PK2
    sigma = kin%convert_PK2_to_Cauchy(S_orig)
    S_recovered = kin%convert_Cauchy_to_PK2(sigma)
    if (.not. (S_recovered .approx. S_orig)) return

    ! 2. PK2 -> Kirchhoff -> PK2
    tau = kin%convert_PK2_to_Kirchhoff(S_orig)
    S_recovered = kin%convert_Kirchhoff_to_PK2(tau)
    if (.not. (S_recovered .approx. S_orig)) return

    ! 3. PK2 -> 1st PK (P) -> PK2
    P_pk1 = kin%convert_PK2_to_PK1(S_orig)
    S_recovered = kin%convert_PK1_to_PK2(P_pk1)
    if (.not. (S_recovered .approx. S_orig)) return

    passed = .TRUE.
end subroutine test_stress_conversions_roundtrip


! =========================================================================
! 8. 4th-Order Tangent Modulus Transformations
! =========================================================================
subroutine test_tangent_modulus_transformations(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_kinematics
    implicit none
    logical, intent(out) :: passed

    type(Total_lagrangian_kinematics) :: kin
    type(ten_3D2O)                     :: F
    type(ten_3D4O2sym)                 :: C_mat, c_spat, C_recovered

    passed = .FALSE.

    ! F = diag(2.0, 1.0, 1.0) => J = 2.0
    call F%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, yx=0.0D0, yy=1.0D0, yz=0.0D0, zx=0.0D0, zy=0.0D0, zz=1.0D0)
    call kin%update(F)

    C_mat%vals = 0.0D0
    C_mat%vals(1,1) = 200.0D0
    C_mat%vals(2,2) = 200.0D0

    ! Push-Forward: c_1111 = (1/J) * F_11^4 * C_1111 = (1/2) * 16 * 200 = 1600.0
    c_spat = kin%push_forward_tangent(C_mat)

    if (abs(c_spat%vals(1,1) - 1600.0D0) > 1.0D-6) then
        print*, "FAIL: 4th-order push-forward tangent incorrect", new_line('A'), &
                "Expected c_1111 = 1600.0, got:", c_spat%vals(1,1)
        return
    end if

    ! Pull-Back: C_recovered = J * F_11^-4 * c_1111 = 2.0 * (0.5)^4 * 1600.0 = 200.0
    C_recovered = kin%pull_back_tangent(c_spat)

    if (.not. (C_recovered .approx. C_mat)) then
        print*, "FAIL: 4th-order tangent round-trip pull-back failed"
        return
    end if

    passed = .TRUE.
end subroutine test_tangent_modulus_transformations