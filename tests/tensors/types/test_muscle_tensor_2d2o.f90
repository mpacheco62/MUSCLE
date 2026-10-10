! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_tensor_2d2o
    !! Unit tests for `ten_2D2O` (general 2D second-order tensor, storage: xx, yx, xy, yy, zz).
    !!
    !! Expected values are exact closed-form results computed by hand from small integer
    !! inputs, so the default relative tolerance of `.approx.` (1e-12) is far above round-off.
    use muscle_tensors
    implicit none
    logical :: passed

    call test_init(passed)
    if (.not. passed) STOP 1

    call test_arithmetic(passed)
    if (.not. passed) STOP 2

    call test_dev_ddot(passed)
    if (.not. passed) STOP 3

    call test_det_transpose(passed)
    if (.not. passed) STOP 4

    call test_products(passed)
    if (.not. passed) STOP 5

    call test_identity_ops(passed)
    if (.not. passed) STOP 6

    call test_mixed_sym(passed)
    if (.not. passed) STOP 7

    print *, "Passed!"
    STOP 0
end program test_muscle_tensor_2d2o


subroutine test_init(passed)
    !! init2 takes (xx, xy, yx, yy, zz) row-wise and stores (xx, yx, xy, yy, zz).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: a, b

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call b%init((/1D0, 3D0, 2D0, 4D0, 5D0/))
    call check_approx(passed, "init(xx=..) vs init(array)", a, b)
end subroutine test_init


subroutine test_arithmetic(passed)
    !! +, -, unary -, scalar *, / and scalar assignment.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: a, b, expected

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call b%init(xx=1D0, xy=1D0, yx=1D0, yy=1D0, zz=1D0)

    call expected%init(xx=2D0, xy=3D0, yx=4D0, yy=5D0, zz=6D0)
    call check_approx(passed, "a + b", a + b, expected)
    call expected%init(xx=0D0, xy=1D0, yx=2D0, yy=3D0, zz=4D0)
    call check_approx(passed, "a - b", a - b, expected)
    call expected%init(xx=-1D0, xy=-2D0, yx=-3D0, yy=-4D0, zz=-5D0)
    call check_approx(passed, "-a", -a, expected)
    call expected%init(xx=2D0, xy=4D0, yx=6D0, yy=8D0, zz=10D0)
    call check_approx(passed, "2*a", 2D0*a, expected)
    call check_approx(passed, "a*2", a*2D0, expected)
    call check_approx(passed, "(2a)/2", expected/2D0, a)
    b = 7D0
    call check_approx(passed, "b = 7 (scalar assignment)", b%vals, spread(7D0, 1, 5), tol=0D0)
end subroutine test_arithmetic


subroutine test_dev_ddot(passed)
    !! dev(A) with tr(A) = 1 + 4 + 4 = 9; A:B = sum A_ij B_ij; norm = sum |A_ij|.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: a, b, expected
    real(real64), parameter :: TOL = 1D-14  ! exact integer arithmetic

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=4D0)
    call expected%init(xx=-2D0, xy=2D0, yx=3D0, yy=1D0, zz=1D0)
    call check_approx(passed, ".dev. a", .dev. a, expected)
    call b%init(xx=1D0, xy=-1D0, yx=2D0, yy=0D0, zz=0.5D0)
    ! 1*1 + 2*(-1) + 3*2 + 4*0 + 4*0.5 = 1 - 2 + 6 + 0 + 2 = 7
    call check_approx(passed, "a .ddot. b", a .ddot. b, 7D0, tol=TOL)
    call check_approx(passed, "a%norm()", a%norm(), 14D0, tol=TOL)
end subroutine test_dev_ddot


subroutine test_det_transpose(passed)
    !! det = (xx*yy - xy*yx)*zz; transpose swaps xy and yx (method and generic).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: a, expected

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call check_approx(passed, "a%det() = (4 - 6)*5", a%det(), -10D0, tol=1D-13)
    call expected%init(xx=1D0, xy=3D0, yx=2D0, yy=4D0, zz=5D0)
    call check_approx(passed, "transpose(a)", transpose(a), expected)
    call check_approx(passed, "a%transpose()", a%transpose(), expected)
end subroutine test_det_transpose


subroutine test_products(passed)
    !! Single contraction: [[1,2],[3,4]].[[0,1],[1,0]] = [[2,1],[4,3]], zz: 5*2 = 10;
    !! det(A.B) = det(A) det(B).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: a, b, ab, expected

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call b%init(xx=0D0, xy=1D0, yx=1D0, yy=0D0, zz=2D0)
    call expected%init(xx=2D0, xy=1D0, yx=4D0, yy=3D0, zz=10D0)
    ab = a * b
    call check_approx(passed, "a * b", ab, expected)
    call check_approx(passed, "det(a*b) = det(a) det(b)", ab%det(), a%det()*b%det(), tol=1D-12)
end subroutine test_products


subroutine test_identity_ops(passed)
    !! I and cI with a general 2D tensor: +, -, *, :.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: a, expected
    type(iden_2O) :: I2
    type(iden_2OS) :: I2S
    real(real64), parameter :: TOL = 1D-14

    passed = .true.
    call I2S%init(2D0)
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)

    call expected%init(xx=2D0, xy=2D0, yx=3D0, yy=5D0, zz=6D0)
    call check_approx(passed, "a + I", a + I2, expected)
    call check_approx(passed, "I + a", I2 + a, expected)
    call expected%init(xx=0D0, xy=2D0, yx=3D0, yy=3D0, zz=4D0)
    call check_approx(passed, "a - I", a - I2, expected)
    call check_approx(passed, "I - a", I2 - a, -expected)
    call expected%init(xx=3D0, xy=2D0, yx=3D0, yy=6D0, zz=7D0)
    call check_approx(passed, "a + 2I", a + I2S, expected)
    call check_approx(passed, "2I + a", I2S + a, expected)
    call expected%init(xx=-1D0, xy=2D0, yx=3D0, yy=2D0, zz=3D0)
    call check_approx(passed, "a - 2I", a - I2S, expected)
    call check_approx(passed, "2I - a", I2S - a, -expected)
    call check_approx(passed, "I * a", I2 * a, a)
    call check_approx(passed, "a * I", a * I2, a)
    call check_approx(passed, "2I * a", I2S * a, 2D0*a)
    call check_approx(passed, "a * 2I", a * I2S, 2D0*a)
    call check_approx(passed, "I .ddot. a", I2 .ddot. a, 10D0, tol=TOL)
    call check_approx(passed, "a .ddot. I", a .ddot. I2, 10D0, tol=TOL)
    call check_approx(passed, "2I .ddot. a", I2S .ddot. a, 20D0, tol=TOL)
    call check_approx(passed, "a .ddot. 2I", a .ddot. I2S, 20D0, tol=TOL)
end subroutine test_identity_ops


subroutine test_mixed_sym(passed)
    !! General/symmetric mixing: assignment, +, -, :, and I^S : A = sym(A).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2O) :: a, s_as_gen, expected
    type(ten_2D2Osym) :: s, sym_a
    type(iden_4O4T) :: I4
    type(iden_4O4TS) :: I4S

    passed = .true.
    call a%init(xx=1D0, xy=2D0, yx=3D0, yy=4D0, zz=5D0)
    call s%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)

    s_as_gen = s
    call expected%init(xx=1D0, xy=4D0, yx=4D0, yy=2D0, zz=3D0)
    call check_approx(passed, "2D2O = 2D2Osym", s_as_gen, expected)

    call check_approx(passed, "a + s", a + s, a + s_as_gen)
    call check_approx(passed, "s + a", s + a, a + s_as_gen)
    call check_approx(passed, "a - s", a - s, a - s_as_gen)
    call check_approx(passed, "s - a", s - a, s_as_gen - a)
    call check_approx(passed, "a .ddot. s", a .ddot. s, a .ddot. s_as_gen, tol=1D-13)
    call check_approx(passed, "s .ddot. a", s .ddot. a, a .ddot. s_as_gen, tol=1D-13)

    call sym_a%init(xx=1D0, yy=4D0, zz=5D0, xy=2.5D0)   ! (xy + yx)/2 = 2.5
    call check_approx(passed, "I^S .ddot. a = sym(a)", I4 .ddot. a, sym_a)
    call check_approx(passed, "a .ddot. I^S = sym(a)", a .ddot. I4, sym_a)
    I4S%val = 2D0
    call check_approx(passed, "2I^S .ddot. a", I4S .ddot. a, 2D0*sym_a)
    call check_approx(passed, "a .ddot. 2I^S", a .ddot. I4S, 2D0*sym_a)
    call check_approx(passed, "2I^S .ddot. s", I4S .ddot. s, 2D0*s)
    call check_approx(passed, "s .ddot. 2I^S", s .ddot. I4S, 2D0*s)
end subroutine test_mixed_sym
