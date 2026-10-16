! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_plastic_history
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_plastic_history
    implicit none

    logical :: passed

    print*, "--------------------------------------------------------"
    print*, "Running Material History Tests..."
    print*, "--------------------------------------------------------"

    call test_history_transactions(passed)
    if (.not. passed) stop 1

    call test_fea_array_packing(passed)
    if (.not. passed) stop 2

    call test_fea_slot_mapping(passed)
    if (.not. passed) stop 3

    call test_history_partial_init(passed)
    if (.not. passed) stop 4

    print*, "Material History Tests PASSED successfully!"
    stop 0
end program test_muscle_plastic_history


subroutine test_history_transactions(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_plastic_history
    implicit none
    logical, intent(out) :: passed

    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain_p0, strain_p_new, stress_new
    real(real64)      :: strain_pf0, strain_pf_new
    real(real64), parameter :: EPS = 1.0D-10

    passed = .FALSE.

    ! 1. Initialize History (t_n)
    call strain_p0%init(xx=0.01D0, yy=-0.005D0, zz=-0.005D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    strain_pf0 = 0.01D0
    call history%init(strain_p = strain_p0, strain_pf = strain_pf0)

    ! 2. Simulate Return Mapping updating state_np1 (candidate)
    call strain_p_new%init(xx=0.05D0, yy=-0.025D0, zz=-0.025D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
    strain_pf_new = 0.05D0
    call stress_new%init(xx=250.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)

    history%state_np1%strain_p  = strain_p_new
    history%state_np1%strain_pf = strain_pf_new
    history%state_np1%stress    = stress_new

    ! 3. Test Rollback (FEA global step failed => state_np1 must reset to state_n)
    call history%rollback()

    if (.not. (history%state_np1%strain_p .approx. strain_p0)) then
        print*, "FAIL: Rollback failed to restore strain_p"
        return
    end if

    if (abs(history%state_np1%strain_pf - strain_pf0) > EPS) then
        print*, "FAIL: Rollback failed to restore strain_pf"
        return
    end if

    ! 4. Re-apply candidate update and Commit (FEA global step converged)
    history%state_np1%strain_p  = strain_p_new
    history%state_np1%strain_pf = strain_pf_new
    history%state_np1%stress    = stress_new

    call history%commit() ! state_n becomes state_np1

    if (.not. (history%state_n%strain_p .approx. strain_p_new)) then
        print*, "FAIL: Commit failed to promote state_np1 to state_n"
        return
    end if

    if (abs(history%state_n%strain_pf - strain_pf_new) > EPS) then
        print*, "FAIL: Commit failed to promote strain_pf"
        return
    end if

    passed = .TRUE.
end subroutine test_history_transactions


subroutine test_fea_array_packing(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_plastic_history
    implicit none
    logical, intent(out) :: passed

    type(Plastic_material_history) :: history
    real(real64) :: hsv_input(13), hsv_output(13)
    real(real64), parameter :: EPS = 1.0D-10

    passed = .FALSE.

    ! 1. Mock ANSYS/LS-DYNA input array hsv_input
    ! Stress: (100, 20, -10, 5, 0, 0)
    hsv_input(1:6)  = (/ 100.0D0, 20.0D0, -10.0D0, 5.0D0, 0.0D0, 0.0D0 /)
    ! Plastic Strain: (0.02, -0.01, -0.01, 0.005, 0, 0)
    hsv_input(7:12) = (/ 0.02D0, -0.01D0, -0.01D0, 0.005D0, 0.0D0, 0.0D0 /)
    ! Equivalent Plastic Strain: 0.02
    hsv_input(13)   = 0.02D0

    ! 2. Unpack into history
    call history%unpack_from_fea(hsv_input)

    if (abs(history%state_n%strain_pf - 0.02D0) > EPS) then
        print*, "FAIL: Unpack equivalent plastic strain failed"
        return
    end if

    ! 3. Pack back to output array
    call history%pack_to_fea(hsv_output)

    if (maxval(abs(hsv_output - hsv_input)) > EPS) then
        print*, "FAIL: Pack/Unpack roundtrip mismatch"
        return
    end if

    passed = .TRUE.
end subroutine test_fea_array_packing

subroutine test_fea_slot_mapping(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_plastic_history
    implicit none
    logical, intent(out) :: passed

    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: stress_n, strain_p_n, stress_new, strain_p_new, packed_stress, packed_strain_p
    real(real64) :: hsv_input(15), hsv_output(15)
    real(real64), parameter :: EPS = 1.0D-10

    passed = .FALSE.

    ! 1. All 13 slots distinct and nonzero, so unpack cannot swap two of them (yz and xz
    !    included) unnoticed. Slots 14 and 15 belong to the FE program.
    hsv_input(1:6)   = (/ 100.0D0, 20.0D0, -10.0D0, 5.0D0, 7.0D0, -3.0D0 /)
    hsv_input(7:12)  = (/ 0.02D0, -0.011D0, -0.009D0, 0.004D0, -0.006D0, 0.003D0 /)
    hsv_input(13)    = 0.025D0
    hsv_input(14:15) = (/ 99.0D0, -99.0D0 /)
    call stress_n%init(xx=100.0D0, yy=20.0D0, zz=-10.0D0, xy=5.0D0, yz=7.0D0, xz=-3.0D0)
    call strain_p_n%init(xx=0.02D0, yy=-0.011D0, zz=-0.009D0, xy=0.004D0, yz=-0.006D0, xz=0.003D0)

    call history%unpack_from_fea(hsv_input)

    if (.not. history%state_n%stress%is_approx(stress_n, tol=EPS)) then
        print*, "FAIL: Unpack stress slots (xx, yy, zz, xy, yz, xz)"
        return
    end if

    if (.not. history%state_n%strain_p%is_approx(strain_p_n, tol=EPS)) then
        print*, "FAIL: Unpack plastic strain slots (xx, yy, zz, xy, yz, xz)"
        return
    end if

    ! 2. pack_to_fea writes the candidate state_np1, not state_n, and leaves slots 14-15 as they are.
    !    Every candidate field differs from state_n and has distinct shears.
    call stress_new%init(xx=250.0D0, yy=-40.0D0, zz=15.0D0, xy=12.0D0, yz=-8.0D0, xz=6.0D0)
    call strain_p_new%init(xx=0.05D0, yy=-0.03D0, zz=-0.02D0, xy=0.008D0, yz=-0.007D0, xz=0.009D0)
    history%state_np1%stress = stress_new
    history%state_np1%strain_p = strain_p_new
    history%state_np1%strain_pf = 0.125D0
    hsv_output = hsv_input
    call history%pack_to_fea(hsv_output)
    call packed_stress%init(hsv_output(1:6))

    if (.not. packed_stress%is_approx(stress_new, tol=EPS)) then
        print*, "FAIL: Pack does not write the candidate stress"
        return
    end if

    call packed_strain_p%init(hsv_output(7:12))
    if (.not. packed_strain_p%is_approx(strain_p_new, tol=EPS) .or. abs(hsv_output(13) - 0.125D0) > EPS) then
        print*, "FAIL: Pack does not write the candidate plastic strains"
        return
    end if

    if (abs(hsv_output(14) - 99.0D0) > EPS .or. abs(hsv_output(15) + 99.0D0) > EPS) then
        print*, "FAIL: Pack changed slots beyond 13"
        return
    end if

    ! 3. Commit and rollback also carry the stress (the existing test checks strains only)
    call history%commit()

    if (.not. history%state_n%stress%is_approx(stress_new, tol=EPS)) then
        print*, "FAIL: Commit failed to promote the stress"
        return
    end if

    history%state_np1%stress = stress_n
    call history%rollback()

    if (.not. history%state_np1%stress%is_approx(stress_new, tol=EPS)) then
        print*, "FAIL: Rollback failed to restore the stress"
        return
    end if

    passed = .TRUE.
end subroutine test_fea_slot_mapping

subroutine test_history_partial_init(passed)
    ! Called by the test program: omitted strain_p must preserve two distinct prior tensors.
    use muscle_tensors
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_test_utils, only : check_approx
    implicit none
    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain_p_n, strain_p_np1, stress

    passed = .true.
    ! Partial init: different committed and candidate plastic strains expose an accidental copy or reset.
    call strain_p_n%init(xx=0.02D0, yy=-0.011D0, zz=-0.009D0, &
                         xy=0.004D0, yz=-0.006D0, xz=0.003D0)
    strain_p_np1 = 2D0*strain_p_n
    call history%init(strain_p=strain_p_n, strain_pf=0.025D0)
    history%state_np1%strain_p = strain_p_np1
    call stress%init(xx=120D0, yy=-30D0, zz=15D0, xy=8D0, yz=-11D0, xz=6D0)
    call history%init(stress=stress, strain_pf=0.125D0)
    ! init only copies values (no arithmetic), so every comparison is exact: tol = 0
    call check_approx(passed, 'partial init: committed plastic strain', history%state_n%strain_p, strain_p_n, tol=0D0)
    call check_approx(passed, 'partial init: candidate plastic strain', history%state_np1%strain_p, strain_p_np1, tol=0D0)
    call check_approx(passed, 'partial init: committed stress', history%state_n%stress, stress, tol=0D0)
    call check_approx(passed, 'partial init: candidate stress', history%state_np1%stress, stress, tol=0D0)
    call check_approx(passed, 'partial init: committed equivalent strain', history%state_n%strain_pf, 0.125D0, tol=0D0)
    call check_approx(passed, 'partial init: candidate equivalent strain', history%state_np1%strain_pf, 0.125D0, tol=0D0)
end subroutine test_history_partial_init
