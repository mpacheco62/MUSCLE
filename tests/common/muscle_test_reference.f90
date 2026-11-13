! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_test_reference
    !! Reference checks of double precision results against reference values stored in
    !! `tests/reference/<area>/<name>_dp.txt` from a validated version (D-021). The references were
    !! captured to migrate the code bit for bit (Phases 0-9); they are now regression checks of the
    !! behaviour, compared within the tolerance of their category:
    !!
    !! - `REFERENCE_TOL_ALGEBRA` (default): direct computations (tensors, laws, kinematics), rounding only;
    !! - `REFERENCE_TOL_FD`: finite differences, which amplify the rounding by eps/h;
    !! - `REFERENCE_TOL_SOLVER`: iterative return mappings (a different rounding may change the Newton path:
    !!   10 times NEWTON_TOL_dp);
    !! - with `status_iters`, the last two values of a record are the status (exact) and the number of
    !!   iterations (within 1).
    !!
    !! The error of a value is relative to the largest magnitude of the real values of its record. A NaN
    !! matches any NaN (sign and payload are not results). With the environment variable
    !! `MUSCLE_REFERENCE_STRICT=1` every value must match bit for bit (the migration check, only meaningful
    !! with the reference configuration: gfortran, Debug); with `MUSCLE_REFERENCE_REPORT=1` the largest
    !! relative error of every record that is not bit-identical is printed (to calibrate tolerances).
    !!
    !! A reference test program calls `reference_init(relpath[, tol])` (`tol`: the default tolerance of its
    !! records), then `reference_check(key, x[, tol][, status_iters])` for every result (a `real(dp)` array;
    !! integers are recorded converted to `real(dp)`, which is exact), and finally `reference_finish()`,
    !! which stops with the first failure code:
    !!
    !! - `stop 1`: a result is outside the tolerance of its category (every differing component is
    !!   printed with its relative error);
    !! - `stop 2`: a result has no reference record;
    !! - `stop 3`: a reference record was not produced by the test (with `reference_finish(partial=.true.,
    !!   expected=n)`: fewer or more than `n` distinct records matched);
    !! - `stop 4`: the reference cannot be read (run the test through ctest, which sets
    !!   `MUSCLE_TEST_REFERENCE_DIR`).
    !!
    !! Regeneration mode: `<test> --write FILE` writes the records to FILE instead of checking them
    !! (only through the `regenerate.sh` of the reference directory, on the reference code).
    !!
    !! Record format (one per line): `<key> | <n> | <IEEE-754 bits of each value, hex>`.
    use, intrinsic :: iso_fortran_env, only : int64, output_unit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use muscle_kinds, only : dp
    implicit none
    private

    public :: reference_init, reference_check, reference_finish
    public :: REFERENCE_TOL_ALGEBRA, REFERENCE_TOL_FD, REFERENCE_TOL_SOLVER

    integer, parameter :: MAXV = 512
        !! Largest number of values of a record.
    integer, parameter :: MAXREC = 4096
        !! Largest number of records of a reference file.
    integer, parameter :: KEYLEN = 160
        !! Largest length of a key.
    real(dp), parameter :: REFERENCE_TOL_ALGEBRA = 1.0e-13_dp
        !! Tolerance of direct computations (rounding only).
    real(dp), parameter :: REFERENCE_TOL_FD = 1.0e-9_dp
        !! Tolerance of finite differences (rounding amplified by eps/h).
    real(dp), parameter :: REFERENCE_TOL_SOLVER = 1.0e-9_dp
        !! Tolerance of iterative return mappings (10 times NEWTON_TOL_dp).

    real(dp) :: default_tol = REFERENCE_TOL_ALGEBRA
        !! Tolerance of the records checked without `tol` (set by `reference_init`).
    logical :: strict = .false.
        !! Bit-for-bit mode (`MUSCLE_REFERENCE_STRICT=1`).
    logical :: report = .false.
        !! Print the largest relative error of every record that is not bit-identical
        !! (`MUSCLE_REFERENCE_REPORT=1`).

    logical :: writing = .false.
        !! Regeneration mode (`--write FILE`).
    integer :: wunit = -1
        !! Unit of the reference being written.
    integer :: nref = 0
        !! Records read from the reference.
    character(len=KEYLEN), allocatable :: ref_key(:)
        !! Key of each record.
    integer, allocatable :: ref_n(:)
        !! Number of values of each record.
    integer(int64), allocatable :: ref_bits(:,:)
        !! Bit patterns of the values.
    logical, allocatable :: ref_used(:)
        !! Record already matched by a computed result.
    integer :: code = 0
        !! Exit code: first failure found.
    integer :: n_checked = 0, n_bad = 0
        !! Statistics.

contains

    subroutine reference_init(relpath, tol)
        !! Reads the reference `$MUSCLE_TEST_REFERENCE_DIR/<relpath>`, or opens the output file in
        !! regeneration mode (`--write FILE`).
        character(len=*), intent(in) :: relpath
            !! Reference path relative to `tests/reference`, e.g. `math/math_operations_dp.txt`.
        real(dp), intent(in), optional :: tol
            !! Default tolerance of the records of the test (default `REFERENCE_TOL_ALGEBRA`).
        character(len=4096) :: arg, dir
        character(len=16384) :: line
        integer :: u, ios, p1, p2, n

        if (present(tol)) default_tol = tol
        call get_environment_variable('MUSCLE_REFERENCE_STRICT', arg, status=ios)
        strict = ios == 0 .and. trim(arg) == '1'
        call get_environment_variable('MUSCLE_REFERENCE_REPORT', arg, status=ios)
        report = ios == 0 .and. trim(arg) == '1'
        allocate(ref_key(MAXREC), ref_n(MAXREC), ref_bits(MAXV, MAXREC), ref_used(MAXREC))
        ref_used = .false.
        if (command_argument_count() >= 2) then
            call get_command_argument(1, arg)
            if (trim(arg) == '--write') then
                call get_command_argument(2, arg)
                open(newunit=wunit, file=trim(arg), status='replace', action='write', iostat=ios)
                if (ios /= 0) then
                    print '(a)', ' FAIL: cannot open for writing ' // trim(arg)
                    stop 4
                end if
                writing = .true.
                write(wunit, '(a)') '# Reference values (double precision), written by the test with --write.'
                write(wunit, '(a)') '# <key> | <n> | <IEEE-754 bits of each value, hex>'
                return
            end if
        end if

        call get_environment_variable('MUSCLE_TEST_REFERENCE_DIR', dir, status=ios)
        if (ios /= 0) then
            print '(a)', ' FAIL: MUSCLE_TEST_REFERENCE_DIR is not set (run the test through ctest)'
            stop 4
        end if
        open(newunit=u, file=trim(dir) // '/' // relpath, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print '(a)', ' FAIL: cannot open reference ' // trim(dir) // '/' // relpath
            stop 4
        end if
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            ! The key may contain ' | ' itself: the last two separators delimit n and the values.
            p2 = index(line, ' | ', back=.true.)
            p1 = index(line(1:max(p2-1, 1)), ' | ', back=.true.)
            if (p1 == 0 .or. p2 == 0 .or. p1 > KEYLEN + 1 .or. nref == MAXREC) then
                print '(a)', ' FAIL: malformed reference line: ' // trim(line)
                stop 4
            end if
            read(line(p1+3:p2-1), *, iostat=ios) n
            if (ios /= 0 .or. n < 1 .or. n > MAXV) then
                print '(a)', ' FAIL: malformed value count: ' // trim(line)
                stop 4
            end if
            nref = nref + 1
            ref_key(nref) = line(1:p1-1)
            ref_n(nref) = n
            read(line(p2+3:), '(*(z16, :, 1x))', iostat=ios) ref_bits(1:n, nref)
            if (ios /= 0) then
                print '(a)', ' FAIL: malformed hex values: ' // trim(line)
                stop 4
            end if
        end do
        close(u)
        if (nref == 0) then
            print '(a)', ' FAIL: empty reference ' // relpath
            stop 4
        end if
    end subroutine reference_init

    subroutine reference_check(key, x, tol, status_iters)
        !! Writes (regeneration) or checks (test) the values of one result.
        character(len=*), intent(in) :: key
            !! Unique key of the result (routine, test function, case).
        real(dp), intent(in) :: x(:)
            !! Values of the result.
        real(dp), intent(in), optional :: tol
            !! Relative tolerance of the record (default: that of `reference_init`).
        logical, intent(in), optional :: status_iters
            !! The last two values are the status (exact) and the number of iterations (within 1).
        integer(int64) :: bits(size(x))
        integer :: i, j, nr
        real(dp) :: expected(size(x)), scale, rel(size(x)), maxrel, t
        logical :: same(size(x)), ok(size(x)), tail

        if (len_trim(key) > KEYLEN .or. size(x) > MAXV .or. size(x) < 1) then
            print '(a)', ' FAIL: key too long or bad record size: ' // trim(key)
            stop 4
        end if
        bits = transfer(x, bits)
        if (writing) then
            write(wunit, '(a, " | ", i0, " |", *(1x, z16.16))') trim(key), size(x), bits
            return
        end if

        n_checked = n_checked + 1
        j = 0
        do i = 1, nref
            if (ref_key(i) == key) then
                j = i
                exit
            end if
        end do
        if (j == 0) then
            print '(a)', ' FAIL: no reference for ' // trim(key)
            if (code == 0) code = 2
            n_bad = n_bad + 1
            return
        end if
        if (ref_used(j)) then
            print '(a)', ' FAIL: key produced twice: ' // trim(key)
            if (code == 0) code = 2
            n_bad = n_bad + 1
            return
        end if
        ref_used(j) = .true.
        if (ref_n(j) /= size(x)) then
            print '(a, i0, a, i0)', ' FAIL: ' // trim(key) // ': ', size(x), ' values, reference has ', ref_n(j)
            if (code == 0) code = 1
            n_bad = n_bad + 1
            return
        end if
        expected = transfer(ref_bits(1:size(x), j), expected)
        ! A NaN matches any NaN: its sign and payload are not a result (they depend on the order of
        ! the operands the compiler emits)
        same = bits == ref_bits(1:size(x), j) .or. (ieee_is_nan(x) .and. ieee_is_nan(expected))
        if (all(same)) return

        t = default_tol
        if (present(tol)) t = tol
        tail = .false.
        if (present(status_iters)) tail = status_iters .and. size(x) >= 3
        nr = size(x)
        if (tail) nr = size(x) - 2
        scale = max(maxval(abs(expected(1:nr)), mask=.not. ieee_is_nan(expected(1:nr))), tiny(1.0_dp))
        rel = 0.0_dp
        where (.not. same) rel = abs(x - expected)/scale
        ! A NaN on one side only is never within a tolerance
        where (.not. same .and. (ieee_is_nan(x) .neqv. ieee_is_nan(expected))) rel = huge(1.0_dp)
        if (strict) then
            ok = same
        else
            ok(1:nr) = same(1:nr) .or. rel(1:nr) <= t
            if (tail) then
                ok(nr+1) = same(nr+1)
                ok(nr+2) = abs(x(nr+2) - expected(nr+2)) <= 1.0_dp
            end if
        end if
        maxrel = maxval(rel(1:nr))
        if (report) print '(a, es9.2, a)', ' REPORT: ' // trim(key) // ': max relative error ', maxrel, &
            merge(' (status/iters differ)', '                      ', tail .and. .not. all(same(nr+1:)))
        if (all(ok)) return

        if (strict) then
            print '(a)', ' FAIL: ' // trim(key) // ' differs from the reference value (bit-for-bit check, ' // &
                'MUSCLE_REFERENCE_STRICT=1)'
        else
            print '(a, es9.2, a)', ' FAIL: ' // trim(key) // ' is outside the tolerance of its category (rel. ', &
                t, ')'
        end if
        print '(a)', '   comp   got                        expected                   |diff|/max|expected|'
        do i = 1, size(x)
            if (ok(i)) cycle
            print '(3x, i4, 2x, es25.17, 2x, es25.17, 2x, es9.2)', i, x(i), expected(i), rel(i)
        end do
        print '(a, es9.2)', '   max relative error ', maxrel
        flush(output_unit)
        if (code == 0) code = 1
        n_bad = n_bad + 1
    end subroutine reference_check

    subroutine reference_finish(partial, expected)
        !! Checks that every reference record was produced and stops with the first failure code.
        logical, intent(in), optional :: partial
            !! The test checks only a subset of the reference (e.g. a new API against the records of
            !! the cases it supports): reference records not produced are not a failure.
        integer, intent(in), optional :: expected
            !! With `partial`: number of distinct reference records the test must have matched (a
            !! skipped case is then a failure, code 3).
        integer :: i
        logical :: subset
        subset = .false.
        if (present(partial)) subset = partial
        if (.not. writing .and. subset .and. present(expected)) then
            if (count(ref_used(1:nref)) /= expected) then
                print '(a, i0, a, i0)', ' FAIL: reference records matched: ', count(ref_used(1:nref)), &
                    ', expected ', expected
                if (code == 0) code = 3
                n_bad = n_bad + 1
            end if
        end if
        if (writing) then
            close(wunit)
            return
        end if
        do i = 1, nref
            if (.not. ref_used(i) .and. .not. subset) then
                print '(a)', ' FAIL: reference record not produced by the test: ' // trim(ref_key(i))
                if (code == 0) code = 3
                n_bad = n_bad + 1
            end if
        end do
        print '(a, i0, a, i0, a, i0, a)', ' reference: ', n_checked, ' results checked, ', &
            nref, ' reference records, ', n_bad, ' failure(s)'
        select case (code)
        case (1)
            stop 1
        case (2)
            stop 2
        case (3)
            stop 3
        end select
    end subroutine reference_finish

end module muscle_test_reference
