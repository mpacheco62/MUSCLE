! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_operator_transform
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none

    logical :: passed

    print*, "Running tests for operator .transform. (A .transform. S)..."

    call test_transform_identity(passed)
    if (.not. passed) stop 1

    call test_transform_scaling(passed)
    if (.not. passed) stop 2

    call test_transform_general(passed)
    if (.not. passed) stop 3

    print*, "All .transform. operator tests passed successfully!"
    stop 0
end program test_operator_transform


subroutine test_transform_identity(passed)
    !! Test 1: Transformation with Identity Tensor (F = I => res = S)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F_identity
    type(ten_3D2Osym) :: S, res

    ! Initialize F as Identity (Column-major storage)
    call F_identity%init(xx=1.0D0, xy=0.0D0, xz=0.0D0, &
                         yx=0.0D0, yy=1.0D0, yz=0.0D0, &
                         zx=0.0D0, zy=0.0D0, zz=1.0D0)

    ! Initialize S with arbitrary values
    call S%init(xx=10.0D0, yy=-5.0D0, zz=2.0D0, xy=3.0D0, yz=-1.0D0, xz=4.0D0)

    ! Perform transformation
    res = F_identity .transform. S

    passed = res .approx. S
    if (.not. passed) then
        print*, "FAIL: Identity Transformation (F = I)", new_line('A'), &
                "Expected:", S%vals, new_line('A'), &
                "Actual  :", res%vals
    end if
end subroutine test_transform_identity


subroutine test_transform_scaling(passed)
    !! Test 2: Diagonal Stretch Transformation (F = diag(2, 3, 4))
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F_stretch
    type(ten_3D2Osym) :: S, res, expected

    ! F = diag(2, 3, 4)
    call F_stretch%init(xx=2.0D0, xy=0.0D0, xz=0.0D0, &
                        yx=0.0D0, yy=3.0D0, yz=0.0D0, &
                        zx=0.0D0, zy=0.0D0, zz=4.0D0)

    ! S = (1, 2, 3, 4, 5, 6) in Voigt (xx, yy, zz, xy, yz, xz)
    call S%init(xx=1.0D0, yy=2.0D0, zz=3.0D0, xy=4.0D0, yz=5.0D0, xz=6.0D0)

    ! Analytically:
    ! R_11 = 2^2 * 1 = 4
    ! R_22 = 3^2 * 2 = 18
    ! R_33 = 4^2 * 3 = 48
    ! R_12 = 2 * 3 * 4 = 24
    ! R_23 = 3 * 4 * 5 = 60
    ! R_13 = 2 * 4 * 6 = 48
    call expected%init(xx=4.0D0, yy=18.0D0, zz=48.0D0, xy=24.0D0, yz=60.0D0, xz=48.0D0)

    res = F_stretch .transform. S

    passed = res .approx. expected
    if (.not. passed) then
        print*, "FAIL: Diagonal Scaling Transformation", new_line('A'), &
                "Expected:", expected%vals, new_line('A'), &
                "Actual  :", res%vals
    end if
end subroutine test_transform_scaling


subroutine test_transform_general(passed)
    !! Test 3: General Shear + Extension Deformation Gradient Transformation
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2O) :: F
    type(ten_3D2Osym) :: S, res, expected

    ! F = [1.0, 0.5, 0.0 ; 0.0, 1.2, 0.0 ; 0.0, 0.0, 0.8]
    call F%init(xx=1.0D0, xy=0.5D0, xz=0.0D0, &
                yx=0.0D0, yy=1.2D0, yz=0.0D0, &
                zx=0.0D0, zy=0.0D0, zz=0.8D0)

    ! S = [100.0, 10.0, 0.0 ; 10.0, 50.0, 0.0 ; 0.0, 0.0, -20.0]
    call S%init(xx=100.0D0, yy=50.0D0, zz=-20.0D0, xy=10.0D0, yz=0.0D0, xz=0.0D0)

    ! Analytical calculation of F * S * F^T:
    ! T = F * S = [105.0, 35.0, 0.0 ; 12.0, 60.0, 0.0 ; 0.0, 0.0, -16.0]
    ! R = T * F^T:
    ! R_11 = 105.0 * 1.0 + 35.0 * 0.5 = 122.5
    ! R_22 = 60.0 * 1.2 = 72.0
    ! R_33 = -16.0 * 0.8 = -12.8
    ! R_12 = 35.0 * 1.2 = 42.0
    ! R_23 = 0.0
    ! R_13 = 0.0
    call expected%init(xx=122.5D0, yy=72.0D0, zz=-12.8D0, xy=42.0D0, yz=0.0D0, xz=0.0D0)

    res = F .transform. S

    passed = res .approx. expected
    if (.not. passed) then
        print*, "FAIL: General Transformation (F * S * F^T)", new_line('A'), &
                "Expected:", expected%vals, new_line('A'), &
                "Actual  :", res%vals, new_line('A'), &
                "Diff    :", res%vals - expected%vals
    end if
end subroutine test_transform_general