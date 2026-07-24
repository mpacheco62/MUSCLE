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