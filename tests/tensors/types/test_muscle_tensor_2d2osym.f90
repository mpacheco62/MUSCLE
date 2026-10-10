! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_tensor_2d2osym
    !! Unit tests for `ten_2D2Osym` (symmetric 2D second-order tensor, Voigt: xx, yy, zz, xy).
    !!
    !! All expected values are exact closed-form results obtained by hand; the inputs are
    !! small integers/dyadic fractions, so every operation is exact (or accurate to a few ulp)
    !! in double precision. The default relative tolerance of `.approx.` (1e-12) is therefore
    !! several orders of magnitude above the round-off level (~1e-16).
    use muscle_tensors
    implicit none
    logical :: passed

    call test_init_access(passed)
    if (.not. passed) STOP 1

    call test_arithmetic(passed)
    if (.not. passed) STOP 2

    call test_dev_ddot_trace(passed)
    if (.not. passed) STOP 3

    call test_square(passed)
    if (.not. passed) STOP 4

    call test_norm_approx(passed)
    if (.not. passed) STOP 5

    call test_det_inv(passed)
    if (.not. passed) STOP 6

    call test_assign(passed)
    if (.not. passed) STOP 7

    print *, "Passed!"
    STOP 0
end program test_muscle_tensor_2d2osym


subroutine test_init_access(passed)
    !! Both `init` forms store (xx, yy, zz, xy) and the accessors return them.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, b

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call b%init((/1D0, 2D0, 3D0, 4D0/))
    call check_approx(passed, "init(xx=..) vs init(array)", a, b)
    call check_approx(passed, "accessor xx()", a%xx(), 1D0, tol=0D0)
    call check_approx(passed, "accessor yy()", a%yy(), 2D0, tol=0D0)
    call check_approx(passed, "accessor zz()", a%zz(), 3D0, tol=0D0)
    call check_approx(passed, "accessor xy()", a%xy(), 4D0, tol=0D0)
end subroutine test_init_access


subroutine test_arithmetic(passed)
    !! +, -, unary -, scalar *, / on exact values.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, b, expected

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call b%init(xx=5D0, yy=-1D0, zz=0.5D0, xy=-2D0)

    call expected%init(xx=6D0, yy=1D0, zz=3.5D0, xy=2D0)
    call check_approx(passed, "a + b", a + b, expected)

    call expected%init(xx=-4D0, yy=3D0, zz=2.5D0, xy=6D0)
    call check_approx(passed, "a - b", a - b, expected)

    call expected%init(xx=-1D0, yy=-2D0, zz=-3D0, xy=-4D0)
    call check_approx(passed, "-a", -a, expected)

    call expected%init(xx=2D0, yy=4D0, zz=6D0, xy=8D0)
    call check_approx(passed, "2*a", 2D0*a, expected)
    call check_approx(passed, "a*2", a*2D0, expected)

    call expected%init(xx=0.5D0, yy=1D0, zz=1.5D0, xy=2D0)
    call check_approx(passed, "a/2", a/2D0, expected)
end subroutine test_arithmetic


subroutine test_dev_ddot_trace(passed)
    !! dev(A) = A - tr(A)/3 I; A:B = sum A_ij B_ij (xy counted twice); I:A = tr(A).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, b, expected
    type(iden_2O) :: I2
    real(real64), parameter :: TOL = 1D-14  ! exact integer arithmetic; round-off ~1e-16

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=6D0, xy=4D0)        ! tr = 9
    call expected%init(xx=-2D0, yy=-1D0, zz=3D0, xy=4D0)
    call check_approx(passed, ".dev. a", .dev. a, expected)

    call b%init(xx=2D0, yy=3D0, zz=-1D0, xy=0.5D0)
    ! 1*2 + 2*3 + 6*(-1) + 2*4*0.5 = 2 + 6 - 6 + 4 = 6
    call check_approx(passed, "a .ddot. b", a .ddot. b, 6D0, tol=TOL)

    call check_approx(passed, "I .ddot. a (trace)", I2 .ddot. a, 9D0, tol=TOL)
    call check_approx(passed, "a .ddot. I (trace)", a .ddot. I2, 9D0, tol=TOL)
end subroutine test_dev_ddot_trace


subroutine test_square(passed)
    !! A.A for A = [[1,4,0],[4,2,0],[0,0,3]] = [[17,12,0],[12,20,0],[0,0,9]].
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, expected

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call expected%init(xx=17D0, yy=20D0, zz=9D0, xy=12D0)
    call check_approx(passed, "a%square()", a%square(), expected)
end subroutine test_square


subroutine test_norm_approx(passed)
    !! Modified L1 norm and the relative tolerance of is_approx.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, b

    passed = .true.
    call a%init(xx=1D0, yy=-2D0, zz=3D0, xy=-4D0)
    call check_approx(passed, "a%norm() = 1 + 2 + 3 + 2*4", a%norm(), 14D0, tol=1D-14)

    ! Perturbation of 1e-9 relative: rejected by default (1e-12), accepted with tol=1e-8
    b = a
    b%vals(4) = b%vals(4) * (1D0 + 1D-9)
    call check_true(passed, ".approx. rejects a 1e-9 relative perturbation", .not. (a .approx. b))
    call check_true(passed, "is_approx(tol=1e-8) accepts a 1e-9 relative perturbation", &
                    a%is_approx(b, tol=1D-8))
end subroutine test_norm_approx


subroutine test_det_inv(passed)
    !! det A = (xx*yy - xy^2)*zz and A.A^-1 = I, with the analytical inverse
    !! A = [[2,1,0],[1,3,0],[0,0,4]] -> det = (6-1)*4 = 20,
    !! A^-1 = [[3/5,-1/5,0],[-1/5,2/5,0],[0,0,1/4]].
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, ainv, expected, zero
    type(ten_2D2O) :: prod, ident

    passed = .true.
    call a%init(xx=2D0, yy=3D0, zz=4D0, xy=1D0)
    call check_approx(passed, "a%det()", a%det(), 20D0, tol=1D-13)

    call expected%init(xx=0.6D0, yy=0.4D0, zz=0.25D0, xy=-0.2D0)
    ainv = .inv. a
    call check_approx(passed, ".inv. a", ainv, expected)
    call check_approx(passed, "a%inv()", a%inv(), expected)

    call ident%init(xx=1D0, xy=0D0, yx=0D0, yy=1D0, zz=1D0)
    prod = a * ainv
    call check_approx(passed, "a * inv(a) = I", prod, ident, tol=1D-12)

    ! Singular tensor (zz = 0): the inverse is defined as the zero tensor
    call a%init(xx=2D0, yy=3D0, zz=0D0, xy=1D0)
    zero = 0D0
    call check_approx(passed, ".inv. (singular) = 0", .inv. a, zero)
end subroutine test_det_inv


subroutine test_assign(passed)
    !! Scalar and identity assignments.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, expected
    type(iden_2O) :: I2
    type(iden_2OS) :: I2S

    passed = .true.
    a = 3D0
    call expected%init(xx=3D0, yy=3D0, zz=3D0, xy=3D0)
    call check_approx(passed, "a = 3 (scalar assignment)", a, expected)

    a = I2
    call expected%init(xx=1D0, yy=1D0, zz=1D0, xy=0D0)
    call check_approx(passed, "a = I", a, expected)

    call I2S%init(2.5D0)
    a = I2S
    call expected%init(xx=2.5D0, yy=2.5D0, zz=2.5D0, xy=0D0)
    call check_approx(passed, "a = 2.5 I", a, expected)
end subroutine test_assign
