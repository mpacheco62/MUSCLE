! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_tensor_2d4o
    !! Unit tests for the 2D fourth-order tensors without major symmetry:
    !! `ten_2D4O2sym` (4x4 Voigt, minor symmetries) and `ten_2D4O` (5x5, general).
    use muscle_tensors
    implicit none
    logical :: passed

    call test_2D4O2sym_basic(passed)
    if (.not. passed) STOP 1

    call test_2D4O2sym_inverse(passed)
    if (.not. passed) STOP 2

    call test_2D4O_get_set_init(passed)
    if (.not. passed) STOP 3

    call test_2D4O_arith_inverse(passed)
    if (.not. passed) STOP 4

    call test_2D4O_convert(passed)
    if (.not. passed) STOP 5

    print *, "Passed!"
    STOP 0
end program test_muscle_tensor_2d4o


subroutine test_2D4O2sym_basic(passed)
    !! init forms, arithmetic, norm and conversion to the fully symmetric type
    !! (exact small-integer values).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O2sym) :: a, b, r
    type(ten_2D4O3sym) :: expected3
    real(real64) :: m(4,4)

    passed = .true.
    m = reshape((/ 1D0,  5D0,  9D0, 13D0, &
                   2D0,  6D0, 10D0, 14D0, &
                   3D0,  7D0, 11D0, 15D0, &
                   4D0,  8D0, 12D0, 16D0 /), (/4,4/))
    call a%init(m)
    call b%init(xxxx=1D0,  xxyy=2D0,  xxzz=3D0,  xxxy=4D0,  &
                yyxx=5D0,  yyyy=6D0,  yyzz=7D0,  yyxy=8D0,  &
                zzxx=9D0,  zzyy=10D0, zzzz=11D0, zzxy=12D0, &
                xyxx=13D0, xyyy=14D0, xyzz=15D0, xyxy=16D0)
    call check_approx(passed, "init(matrix) vs init(xxxx=..)", a, b)

    r = a + b;   call check_approx(passed, "a + b", r%vals, 2D0*m, tol=0D0)
    r = a - b;   call check_approx(passed, "a - b", r%vals, 0D0*m, tol=0D0)
    r = -a;      call check_approx(passed, "-a", r%vals, -m, tol=0D0)
    r = 2D0*a;   call check_approx(passed, "2*a", r%vals, 2D0*m, tol=0D0)
    r = a*2D0;   call check_approx(passed, "a*2", r%vals, 2D0*m, tol=0D0)
    r = a/2D0;   call check_approx(passed, "a/2", r%vals, 0.5D0*m, tol=0D0)
    call check_approx(passed, "a%norm() = sum(1..16)", a%norm(), 136D0, tol=1D-15)

    ! Upper triangle in storage order (11,22,33,44,12,23,34,13,24,14)
    call expected3%init((/1D0, 6D0, 11D0, 16D0, 2D0, 7D0, 12D0, 3D0, 8D0, 4D0/))
    call check_approx(passed, "a%convert_3sym()", a%convert_3sym(), expected3)
end subroutine test_2D4O2sym_basic


subroutine test_2D4O2sym_inverse(passed)
    !! Property: for an invertible minor-symmetric A, inv(A) : (A : x) = x and
    !! A : (inv(A) : x) = x for any symmetric x (the shear scaling of the inverse is
    !! what makes this hold with the tensorial-shear .ddot.).
    !! A is diagonally dominant (cond ~ 3), so LU round-off is ~1e-15; tol 1e-12.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O2sym) :: a, ainv
    type(ten_2D2Osym) :: x

    passed = .true.
    call a%init(xxxx=10D0, xxyy=1D0,  xxzz=2D0,  xxxy=0.5D0, &
                yyxx=-1D0, yyyy=8D0,  yyzz=1D0,  yyxy=1D0,   &
                zzxx=0.5D0, zzyy=2D0, zzzz=9D0,  zzxy=-1D0,  &
                xyxx=1D0,  xyyy=-0.5D0, xyzz=0.25D0, xyxy=4D0)
    ainv = .inv. a
    call x%init(xx=1D0, yy=-2D0, zz=3D0, xy=0.7D0)

    call check_approx(passed, "inv(a) : (a : x) = x", ainv .ddot. (a .ddot. x), x, tol=1D-12)
    call check_approx(passed, "a : (inv(a) : x) = x", a .ddot. (ainv .ddot. x), x, tol=1D-12)
end subroutine test_2D4O2sym_inverse


subroutine test_2D4O_get_set_init(passed)
    !! Index mapping (1,1)->1, (2,1)->2, (1,2)->3, (2,2)->4, (3,3)->5; out-of-plane
    !! coupling indices read as zero and are ignored on set; init from a 3x3x3x3 array.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O) :: a, b
    real(real64) :: C4(3,3,3,3), expected(5,5), got_get(5,5)
    integer :: i, j, k, l, row, col
    integer, parameter :: ip(5) = (/1, 2, 1, 2, 3/)  ! (i,j) of each 5D index
    integer, parameter :: jp(5) = (/1, 1, 2, 2, 3/)

    passed = .true.
    a = 0D0
    call a%set(1, 2, 2, 1, 7D0)
    call a%set(3, 3, 1, 1, -2D0)
    call a%set(1, 3, 1, 1, 99D0)          ! out-of-plane coupling: ignored
    call check_approx(passed, "set(1,2,2,1) -> vals(3,2)", a%vals(3, 2), 7D0, tol=0D0)
    call check_approx(passed, "set(3,3,1,1) -> vals(5,1)", a%vals(5, 1), -2D0, tol=0D0)
    call check_approx(passed, "get(1,2,2,1)", a%get(1, 2, 2, 1), 7D0, tol=0D0)
    call check_approx(passed, "get(1,3,1,1) (out of plane) = 0", a%get(1, 3, 1, 1), 0D0, tol=0D0)
    call check_approx(passed, "a%norm() after set (99 ignored)", a%norm(), 9D0, tol=1D-14)

    ! C4_ijkl = 1000 i + 100 j + 10 k + l (only the 2D entries survive)
    do l = 1, 3; do k = 1, 3; do j = 1, 3; do i = 1, 3
        C4(i,j,k,l) = 1000D0*i + 100D0*j + 10D0*k + l
    end do; end do; end do; end do
    call b%init(C4)
    call check_approx(passed, "init(C4) then get(2,3,1,1) (out of plane) = 0", &
                      b%get(2, 3, 1, 1), 0D0, tol=0D0)
    ! All 25 stored components must be exact copies of C4 (integer values, so tol = 0)
    do col = 1, 5
        do row = 1, 5
            expected(row, col) = C4(ip(row), jp(row), ip(col), jp(col))
            got_get(row, col) = b%get(ip(row), jp(row), ip(col), jp(col))
        end do
    end do
    call check_approx(passed, "init(C4): vals", b%vals, expected, tol=0D0)
    call check_approx(passed, "init(C4): get(i,j,k,l)", got_get, expected, tol=0D0)
end subroutine test_2D4O_get_set_init


subroutine test_2D4O_arith_inverse(passed)
    !! Arithmetic and .inv.: A * inv(A) = I_5 (matrix product of the 5x5 representations).
    !! A is diagonally dominant (cond ~ 3), LU round-off ~1e-15; tol 1e-12.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O) :: a, ainv, zero, r
    real(real64) :: id5(5,5)
    integer :: i

    passed = .true.
    a = 0.5D0
    do i = 1, 5
        a%vals(i, i) = 5D0 + i
    end do
    r = a + a;   call check_approx(passed, "a + a", r%vals, 2D0*a%vals, tol=0D0)
    r = a - a;   call check_approx(passed, "a - a", r%vals, 0D0*a%vals, tol=0D0)
    r = -a;      call check_approx(passed, "-a", r%vals, -a%vals, tol=0D0)
    r = 3D0*a;   call check_approx(passed, "3*a", r%vals, 3D0*a%vals, tol=0D0)
    r = a*3D0;   call check_approx(passed, "a*3", r%vals, 3D0*a%vals, tol=0D0)
    r = a/2D0;   call check_approx(passed, "a/2", r%vals, 0.5D0*a%vals, tol=0D0)
    call check_true(passed, "a .approx. a", a .approx. a)

    id5 = 0D0
    do i = 1, 5
        id5(i, i) = 1D0
    end do
    ainv = .inv. a
    call check_approx(passed, "a * inv(a) = I_5", matmul(a%vals, ainv%vals), id5, tol=1D-12)

    ! Singular matrix -> zero tensor
    zero = 0D0
    a = 1D0
    call check_approx(passed, ".inv. (singular) = 0", .inv. a, zero)
end subroutine test_2D4O_arith_inverse


subroutine test_2D4O_convert(passed)
    !! convert_2sym samples rows/cols (xx, yy, zz, xy) = 5D indices (1, 4, 5, 3).
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O) :: a
    type(ten_2D4O2sym) :: c
    real(real64) :: expected(4,4)
    integer :: i, j
    integer, parameter :: map4(4) = (/1, 4, 5, 3/)

    passed = .true.
    do j = 1, 5
        do i = 1, 5
            a%vals(i, j) = 10D0*i + j
        end do
    end do
    do j = 1, 4
        do i = 1, 4
            expected(i, j) = 10D0*map4(i) + map4(j)
        end do
    end do
    c = a%convert_2sym()
    call check_approx(passed, "a%convert_2sym()", c%vals, expected, tol=0D0)
end subroutine test_2D4O_convert
