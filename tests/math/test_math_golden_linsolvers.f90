! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_math_golden_linsolvers
    !! Golden test of the public routines of `muscle_math_linsolvers` (double precision):
    !! `gauss_solve`, `gauss_seidel` and `mklsolve` (dgesv of the LAPACK the build links), on
    !! well-conditioned systems of order 1..20 (diagonally dominant, and general ones that need
    !! row pivoting), plus the error paths (size mismatch, singular matrix, zero diagonal, no
    !! convergence within `maxit`). Solutions, iteration counts and `info` (integers recorded as
    !! `real(dp)`) are compared bit for bit with `tests/reference/math/math_linsolvers_dp.txt`,
    !! captured from the original code (`tests/reference/math/regenerate.sh`).
    !!
    !! Note: the `mklsolve` records depend on the LAPACK implementation linked (reference LAPACK
    !! or an optimized one may round differently); regenerate them only on the reference code.
    use muscle_kinds, only : dp
    use muscle_math_linsolvers, only : gauss_solve, gauss_seidel, mklsolve
    use muscle_test_helpers, only : park_miller, rnd
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    implicit none

    integer, parameter :: NSIZES = 9
    integer, parameter :: SIZES(NSIZES) = [1, 2, 3, 4, 6, 7, 9, 12, 20]
    real(dp), parameter :: GS_TOL = 1.0e-12_dp
        !! Residual tolerance of Gauss-Seidel.
    integer, parameter :: GS_MAXIT = 500
        !! Iteration limit of Gauss-Seidel.
    type(park_miller) :: rng
    real(dp), allocatable :: a(:,:), b(:), x(:), a0(:,:), b0(:)
    character(len=32) :: c
    integer :: k, n, info, its, i

    call golden_init('math/math_linsolvers_dp.txt')
    rng = park_miller(9)

    do k = 1, NSIZES
        n = SIZES(k)
        write(c, '(" | n = ", i0)') n
        call new_system(n, .true.)
        call golden_check('input dominant A, b' // trim(c), [reshape(a, [n*n]), b])
        ! gauss_solve (its arguments are intent(inout): record that they are not modified)
        x = 0.0_dp
        call gauss_solve(a, b, x, info)
        call golden_check('gauss_solve dominant x, info, A&b unchanged' // trim(c), &
                          [x, real(info, dp), merge(1.0_dp, 0.0_dp, all(a == a0) .and. all(b == b0))])
        ! gauss_seidel from x = 0
        x = 0.0_dp
        call gauss_seidel(a, b, x, GS_TOL, GS_MAXIT, its, info)
        call golden_check('gauss_seidel dominant x0 = 0: x, its, info' // trim(c), [x, real(its, dp), real(info, dp)])
        ! gauss_seidel from a random x
        do i = 1, n
            x(i) = rnd(rng)
        end do
        call golden_check('input gauss_seidel x0' // trim(c), x)
        call gauss_seidel(a, b, x, GS_TOL, GS_MAXIT, its, info)
        call golden_check('gauss_seidel dominant x0 random: x, its, info' // trim(c), &
                          [x, real(its, dp), real(info, dp)])
        ! mklsolve
        x = 0.0_dp
        call mklsolve(a, b, x, GS_TOL, GS_MAXIT, its, info)
        call golden_check('mklsolve dominant x, its, info' // trim(c), [x, real(its, dp), real(info, dp)])

        ! General (pivoting) system, direct solvers only
        call new_system(n, .false.)
        call golden_check('input general A, b' // trim(c), [reshape(a, [n*n]), b])
        x = 0.0_dp
        call gauss_solve(a, b, x, info)
        call golden_check('gauss_solve general x, info' // trim(c), [x, real(info, dp)])
        x = 0.0_dp
        call mklsolve(a, b, x, GS_TOL, GS_MAXIT, its, info)
        call golden_check('mklsolve general x, its, info' // trim(c), [x, real(its, dp), real(info, dp)])
    end do

    ! Gauss-Seidel stopped by maxit (info = 1): x after 3 sweeps
    call new_system(6, .true.)
    x = 0.0_dp
    call gauss_seidel(a, b, x, GS_TOL, 3, its, info)
    call golden_check('gauss_seidel maxit = 3: x, its, info', [x, real(its, dp), real(info, dp)])

    ! Zero diagonal in row 2 (gauss_seidel info = 2, x(1) already updated); pivoting solves it
    call new_system(4, .true.)
    a(2, 2) = 0.0_dp
    x = 0.0_dp
    call gauss_seidel(a, b, x, GS_TOL, GS_MAXIT, its, info)
    call golden_check('gauss_seidel zero diagonal: x, its, info', [x, real(its, dp), real(info, dp)])
    x = 0.0_dp
    call gauss_solve(a, b, x, info)
    call golden_check('gauss_solve zero diagonal: x, info', [x, real(info, dp)])

    ! Singular matrices: zero first column, and a repeated row
    call new_system(5, .false.)
    a(:, 1) = 0.0_dp
    call gauss_solve(a, b, x, info)
    call golden_check('gauss_solve zero column: info', [real(info, dp)])
    call mklsolve(a, b, x, GS_TOL, GS_MAXIT, its, info)
    call golden_check('mklsolve zero column: its, info', [real(its, dp), real(info, dp)])
    call new_system(3, .false.)
    a(3, :) = 0.0_dp
    call gauss_solve(a, b, x, info)
    call golden_check('gauss_solve zero last row: info', [real(info, dp)])

    ! Size mismatch (info = -1)
    call new_system(3, .true.)
    deallocate(x)
    allocate(x(4))
    call gauss_solve(a, b, x, info)
    call golden_check('gauss_solve size mismatch: info', [real(info, dp)])
    call gauss_seidel(a, b, x, GS_TOL, GS_MAXIT, its, info)
    call golden_check('gauss_seidel size mismatch: its, info', [real(its, dp), real(info, dp)])
    call mklsolve(a, b, x, GS_TOL, GS_MAXIT, its, info)
    call golden_check('mklsolve size mismatch: its, info', [real(its, dp), real(info, dp)])

    call golden_finish()

contains

    subroutine new_system(n, dominant)
        !! Pseudo-random system of order n: entries in (-1, 1); diagonally dominant (diagonal
        !! + (n + 1), sign kept) or general (diagonal scaled by 0.1, so that pivoting is needed).
        integer, intent(in) :: n
            !! Order.
        logical, intent(in) :: dominant
            !! Kind of matrix.
        integer :: i, j
        if (allocated(a)) deallocate(a, b, x, a0, b0)
        allocate(a(n,n), b(n), x(n), a0(n,n), b0(n))
        do j = 1, n
            do i = 1, n
                a(i, j) = rnd(rng)
            end do
        end do
        do i = 1, n
            if (dominant) then
                a(i, i) = a(i, i) + sign(real(n + 1, dp), a(i, i))
            else
                a(i, i) = 0.1_dp*a(i, i)
            end if
            b(i) = 10.0_dp*rnd(rng)
        end do
        a0 = a
        b0 = b
    end subroutine new_system

end program test_math_golden_linsolvers
