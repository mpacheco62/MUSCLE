! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_test_utils
    !! Assertions shared by the unit tests.
    !!
    !! Every check receives the accumulated `passed` flag, a label naming the operation under
    !! test, the obtained value and the expected one. The flag is AND-ed with the result, so a
    !! test can run all of its checks and report every failure instead of stopping at the first.
    !!
    !! On failure the check prints the label, the acceptance criterion with its numeric
    !! threshold, and one line per component outside the tolerance with the obtained value,
    !! the expected value and the absolute and relative differences:
    !!
    !! ```text
    !!  FAIL: a + b [ten_2D2Osym]
    !!    criterion: |got - expected| <= tol * scale = 1.00E-12 * 6.00E+00 = 6.00E-12
    !!    scale = max(max|got|, max|expected|)
    !!    1 of 4 component(s) outside tolerance:
    !!    comp        got                       expected                  |diff|     |diff|/scale
    !!    xy           2.000000000000000E+00     2.500000000000000E+00   5.00E-01   8.33E-02
    !! ```
    !!
    !! Acceptance criteria:
    !!
    !! - Tensors: the verdict is the type's own `is_approx` (same default tolerance as
    !!   `.approx.`), i.e. \( |a_I - b_I| \le tol \cdot \max(\max|a|, \max|b|, 10^{-30}) \).
    !! - Real arrays: the same component-wise criterion; `tol = 0` means exact equality.
    !! - Real scalars: \( |a - b| \le tol \cdot \max(1, |a|, |b|) \) (absolute below 1,
    !!   relative above); `tol = 0` means exact equality.
    use, intrinsic :: iso_fortran_env, only : real64, output_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_tensors
    implicit none
    private

    public :: check_approx, check_true

    real(real64), parameter :: DEFAULT_TOL = 1.0D-12
        !! Default relative tolerance (same as `.approx.` for most tensor types).
    real(real64), parameter :: DEFAULT_TOL_4O2SYM = 1.0D-8
        !! Default tolerance of `.approx.` for the minor-symmetric 4th-order types.
    real(real64), parameter :: EPS_ABS = 1.0D-30
        !! Lower bound of the scale (avoids a zero threshold when both operands vanish).

    character(len=2), parameter :: N2SYM(4) = [character(len=2) :: 'xx', 'yy', 'zz', 'xy']
        !! Voigt labels of `ten_2D2Osym` (also rows/columns of `ten_2D4O2sym`).
    character(len=2), parameter :: N2GEN(5) = [character(len=2) :: 'xx', 'yx', 'xy', 'yy', 'zz']
        !! Storage labels of `ten_2D2O` (also rows/columns of `ten_2D4O`).
    character(len=4), parameter :: N2D4O3(10) = [character(len=4) :: 'xxxx', 'yyyy', 'zzzz', &
        'xyxy', 'xxyy', 'yyzz', 'zzxy', 'xxzz', 'yyxy', 'xxxy']
        !! Storage labels of `ten_2D4O3sym`.
    character(len=2), parameter :: N3SYM(6) = [character(len=2) :: 'xx', 'yy', 'zz', 'xy', 'yz', 'xz']
        !! Voigt labels of `ten_3D2Osym`.
    character(len=2), parameter :: N3GEN(9) = [character(len=2) :: 'xx', 'yx', 'zx', 'xy', 'yy', &
        'zy', 'xz', 'yz', 'zz']
        !! Column-major labels of `ten_3D2O`.

    interface check_approx
        !! Approximate (or, with `tol = 0`, exact) comparison with a detailed failure report.
        module procedure check_real
        module procedure check_real_1d
        module procedure check_real_2d
        module procedure check_2D2Osym
        module procedure check_2D2O
        module procedure check_2D4O3sym
        module procedure check_2D4O2sym
        module procedure check_2D4O
        module procedure check_3D2Osym
        module procedure check_3D2O
    end interface check_approx

contains

    subroutine check_true(passed, label, cond)
        !! Checks a logical property (e.g. that two tensors are *not* approximately equal).
        logical, intent(inout) :: passed
            !! Accumulated test result; set to `.false.` if `cond` is false.
        character(len=*), intent(in) :: label
            !! Description of the property, printed on failure.
        logical, intent(in) :: cond
            !! Property to check.
        if (.not. cond) then
            print '(a)', ' FAIL: ' // label
            print '(a)', '   expected condition evaluated to .false.'
            flush(output_unit)
        end if
        passed = passed .and. cond
    end subroutine check_true

    subroutine check_real(passed, label, got, expected, tol)
        !! Scalar comparison: \( |a - b| \le tol \cdot \max(1, |a|, |b|) \).
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        real(real64), intent(in) :: got
            !! Obtained value.
        real(real64), intent(in) :: expected
            !! Expected value.
        real(real64), optional, intent(in) :: tol
            !! Tolerance (default `1E-12`; `0` = exact).
        real(real64) :: t, scale, diff
        logical :: ok

        t = DEFAULT_TOL
        if (present(tol)) t = tol
        scale = max(1.0D0, abs(got), abs(expected))
        diff = abs(got - expected)
        ok = diff <= t*scale
        if (.not. ok) then
            print '(a)', ' FAIL: ' // label
            call print_criterion(t, scale, 'max(1, |got|, |expected|)')
            print '(a)', '   got                       expected                  |diff|     |diff|/scale'
            print '(3x, es24.15, 2x, es24.15, 2x, es9.2, 2x, es9.2)', got, expected, diff, diff/scale
            flush(output_unit)
        end if
        passed = passed .and. ok
    end subroutine check_real

    subroutine check_real_1d(passed, label, got, expected, tol)
        !! Component-wise comparison of two vectors relative to their largest component.
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        real(real64), intent(in) :: got(:)
            !! Obtained values.
        real(real64), intent(in) :: expected(:)
            !! Expected values (same size as `got`).
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default `1E-12`; `0` = exact).
        real(real64) :: t
        character(len=16) :: names(size(got))
        integer :: i

        t = DEFAULT_TOL
        if (present(tol)) t = tol
        do i = 1, size(got)
            write(names(i), '("(", i0, ")")') i
        end do
        call compare(passed, label, got, expected, names, t)
    end subroutine check_real_1d

    subroutine check_real_2d(passed, label, got, expected, tol)
        !! Component-wise comparison of two matrices relative to their largest entry.
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        real(real64), intent(in) :: got(:,:)
            !! Obtained values.
        real(real64), intent(in) :: expected(:,:)
            !! Expected values (same shape as `got`).
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default `1E-12`; `0` = exact).
        real(real64) :: t
        character(len=16) :: names(size(got))
        integer :: i, j

        t = DEFAULT_TOL
        if (present(tol)) t = tol
        do j = 1, size(got, 2)
            do i = 1, size(got, 1)
                write(names(i + (j-1)*size(got, 1)), '("(", i0, ",", i0, ")")') i, j
            end do
        end do
        call compare(passed, label, reshape(got, [size(got)]), reshape(expected, [size(expected)]), &
                     names, t)
    end subroutine check_real_2d

    subroutine check_2D2Osym(passed, label, got, expected, tol)
        !! Compares two `ten_2D2Osym` (verdict: `is_approx`).
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        type(ten_2D2Osym), intent(in) :: got
            !! Obtained tensor.
        type(ten_2D2Osym), intent(in) :: expected
            !! Expected tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default: that of `.approx.`, `1E-12`).
        real(real64) :: t
        t = DEFAULT_TOL
        if (present(tol)) t = tol
        call report(passed, label // ' [ten_2D2Osym]', got%is_approx(expected, tol=t), &
                    got%vals, expected%vals, N2SYM, t)
    end subroutine check_2D2Osym

    subroutine check_2D2O(passed, label, got, expected, tol)
        !! Compares two `ten_2D2O` (verdict: `is_approx`).
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        type(ten_2D2O), intent(in) :: got
            !! Obtained tensor.
        type(ten_2D2O), intent(in) :: expected
            !! Expected tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default: that of `.approx.`, `1E-12`).
        real(real64) :: t
        t = DEFAULT_TOL
        if (present(tol)) t = tol
        call report(passed, label // ' [ten_2D2O]', got%is_approx(expected, tol=t), &
                    got%vals, expected%vals, N2GEN, t)
    end subroutine check_2D2O

    subroutine check_2D4O3sym(passed, label, got, expected, tol)
        !! Compares two `ten_2D4O3sym` (verdict: `is_approx`).
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        type(ten_2D4O3sym), intent(in) :: got
            !! Obtained tensor.
        type(ten_2D4O3sym), intent(in) :: expected
            !! Expected tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default: that of `.approx.`, `1E-12`).
        real(real64) :: t
        t = DEFAULT_TOL
        if (present(tol)) t = tol
        call report(passed, label // ' [ten_2D4O3sym]', got%is_approx(expected, tol=t), &
                    got%vals, expected%vals, N2D4O3, t)
    end subroutine check_2D4O3sym

    subroutine check_2D4O2sym(passed, label, got, expected, tol)
        !! Compares two `ten_2D4O2sym` (verdict: `is_approx`).
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        type(ten_2D4O2sym), intent(in) :: got
            !! Obtained tensor.
        type(ten_2D4O2sym), intent(in) :: expected
            !! Expected tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default: that of `.approx.`, `1E-8`).
        real(real64) :: t
        t = DEFAULT_TOL_4O2SYM
        if (present(tol)) t = tol
        call report(passed, label // ' [ten_2D4O2sym]', got%is_approx(expected, tol=t), &
                    reshape(got%vals, [16]), reshape(expected%vals, [16]), pair_names(N2SYM), t)
    end subroutine check_2D4O2sym

    subroutine check_2D4O(passed, label, got, expected, tol)
        !! Compares two `ten_2D4O` (verdict: `is_approx`).
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        type(ten_2D4O), intent(in) :: got
            !! Obtained tensor.
        type(ten_2D4O), intent(in) :: expected
            !! Expected tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default: that of `.approx.`, `1E-12`).
        real(real64) :: t
        t = DEFAULT_TOL
        if (present(tol)) t = tol
        call report(passed, label // ' [ten_2D4O]', got%is_approx(expected, tol=t), &
                    reshape(got%vals, [25]), reshape(expected%vals, [25]), pair_names(N2GEN), t)
    end subroutine check_2D4O

    subroutine check_3D2Osym(passed, label, got, expected, tol)
        !! Compares two `ten_3D2Osym` (verdict: `is_approx`).
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        type(ten_3D2Osym), intent(in) :: got
            !! Obtained tensor.
        type(ten_3D2Osym), intent(in) :: expected
            !! Expected tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default: that of `.approx.`, `1E-12`).
        real(real64) :: t
        t = DEFAULT_TOL
        if (present(tol)) t = tol
        call report(passed, label // ' [ten_3D2Osym]', got%is_approx(expected, tol=t), &
                    got%vals, expected%vals, N3SYM, t)
    end subroutine check_3D2Osym

    subroutine check_3D2O(passed, label, got, expected, tol)
        !! Compares two `ten_3D2O` (verdict: `is_approx`).
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        type(ten_3D2O), intent(in) :: got
            !! Obtained tensor.
        type(ten_3D2O), intent(in) :: expected
            !! Expected tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default: that of `.approx.`, `1E-12`).
        real(real64) :: t
        t = DEFAULT_TOL
        if (present(tol)) t = tol
        call report(passed, label // ' [ten_3D2O]', got%is_approx(expected, tol=t), &
                    got%vals, expected%vals, N3GEN, t)
    end subroutine check_3D2O

    ! -------------------------------------------------------------------------
    ! Private helpers
    ! -------------------------------------------------------------------------

    subroutine compare(passed, label, got, expected, names, tol)
        !! Applies the component-wise relative criterion and reports on failure.
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        real(real64), intent(in) :: got(:)
            !! Obtained components.
        real(real64), intent(in) :: expected(:)
            !! Expected components.
        character(len=*), intent(in) :: names(:)
            !! Component labels.
        real(real64), intent(in) :: tol
            !! Relative tolerance.
        real(real64) :: scale
        scale = max(maxval(abs(got)), maxval(abs(expected)), EPS_ABS)
        call report(passed, label, all(abs(got - expected) <= tol*scale), got, expected, names, tol)
    end subroutine compare

    subroutine report(passed, label, ok, got, expected, names, tol)
        !! Updates `passed` and, if `ok` is false, prints the components outside the tolerance.
        logical, intent(inout) :: passed
            !! Accumulated test result.
        character(len=*), intent(in) :: label
            !! Operation under test.
        logical, intent(in) :: ok
            !! Verdict of the comparison.
        real(real64), intent(in) :: got(:)
            !! Obtained components.
        real(real64), intent(in) :: expected(:)
            !! Expected components.
        character(len=*), intent(in) :: names(:)
            !! Component labels.
        real(real64), intent(in) :: tol
            !! Relative tolerance used for the verdict.
        real(real64) :: scale, diff
        integer :: i, n_bad
        character(len=10) :: comp

        passed = passed .and. ok
        if (ok) return

        scale = max(maxval(abs(got)), maxval(abs(expected)), EPS_ABS)
        n_bad = count(abs(got - expected) > tol*scale)
        print '(a)', ' FAIL: ' // label
        call print_criterion(tol, scale, 'max(max|got|, max|expected|)')
        if (any(ieee_is_nan(got)) .or. any(ieee_is_nan(expected))) then
            print '(a)', '   (NaN present in the operands)'
        end if
        print '(3x, i0, a, i0, a)', n_bad, ' of ', size(got), ' component(s) outside tolerance:'
        print '(a)', '   comp        got                       expected                  |diff|     |diff|/scale'
        do i = 1, size(got)
            diff = abs(got(i) - expected(i))
            if (.not. (diff <= tol*scale)) then
                comp = names(i)
                print '(3x, a10, es24.15, 2x, es24.15, 2x, es9.2, 2x, es9.2)', &
                    comp, got(i), expected(i), diff, diff/scale
            end if
        end do
        flush(output_unit)
    end subroutine report

    subroutine print_criterion(tol, scale, scale_def)
        !! Prints the acceptance criterion with its numeric threshold.
        real(real64), intent(in) :: tol
            !! Tolerance.
        real(real64), intent(in) :: scale
            !! Scale multiplying the tolerance.
        character(len=*), intent(in) :: scale_def
            !! Definition of the scale.
        print '(a, es8.2, a, es8.2, a, es8.2)', &
            '   criterion: |got - expected| <= tol * scale = ', tol, ' * ', scale, ' = ', tol*scale
        print '(a)', '   scale = ' // scale_def
    end subroutine print_criterion

    pure function pair_names(base) result(names)
        !! Labels `ab:cd` of a matrix whose rows and columns are labelled by `base`
        !! (column-major, as the matrix is flattened).
        character(len=2), intent(in) :: base(:)
            !! Row/column labels.
        character(len=5) :: names(size(base)**2)
            !! Component labels.
        integer :: i, j
        do j = 1, size(base)
            do i = 1, size(base)
                names(i + (j-1)*size(base)) = base(i) // ':' // base(j)
            end do
        end do
    end function pair_names

end module muscle_test_utils
