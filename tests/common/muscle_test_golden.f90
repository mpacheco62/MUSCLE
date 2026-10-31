! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_test_golden
    !! Bit-for-bit golden (characterization) checks of double precision results against a
    !! reference file captured from the original code (`tests/reference/<area>/<name>_dp.txt`).
    !!
    !! A golden test program calls `golden_init(relpath)`, then `golden_check(key, x)` for every
    !! result (a `real(dp)` array; integers are recorded converted to `real(dp)`, which is exact),
    !! and finally `golden_finish()`, which stops with the first failure code:
    !!
    !! - `stop 1`: a result differs from its reference (every differing component is printed with
    !!   its relative error; <= 1e-14 is reported as approvable, testing.md, but still fails);
    !! - `stop 2`: a result has no reference record;
    !! - `stop 3`: a reference record was not produced by the test;
    !! - `stop 4`: the reference cannot be read (run the test through ctest, which sets
    !!   `MUSCLE_TEST_REFERENCE_DIR`).
    !!
    !! Regeneration mode: `<test> --write FILE` writes the records to FILE instead of checking them
    !! (only through the `regenerate.sh` of the reference directory, on the reference code).
    !!
    !! Record format (one per line): `<key> | <n> | <IEEE-754 bits of each value, hex>`.
    use, intrinsic :: iso_fortran_env, only : int64, output_unit
    use muscle_kinds, only : dp
    implicit none
    private

    public :: golden_init, golden_check, golden_finish

    integer, parameter :: MAXV = 512
        !! Largest number of values of a record.
    integer, parameter :: MAXREC = 2048
        !! Largest number of records of a reference file.
    integer, parameter :: KEYLEN = 160
        !! Largest length of a key.
    real(dp), parameter :: ACCEPTABLE_REL = 1.0e-14_dp
        !! Relative error reported as "within the approvable bound" of testing.md (still a failure).

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

    subroutine golden_init(relpath)
        !! Reads the reference `$MUSCLE_TEST_REFERENCE_DIR/<relpath>`, or opens the output file in
        !! regeneration mode (`--write FILE`).
        character(len=*), intent(in) :: relpath
            !! Reference path relative to `tests/reference`, e.g. `math/math_derivatives_dp.txt`.
        character(len=4096) :: arg, dir
        character(len=16384) :: line
        integer :: u, ios, p1, p2, n

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
                write(wunit, '(a)') '# Golden values (double precision), written by the test with --write.'
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
    end subroutine golden_init

    subroutine golden_check(key, x)
        !! Writes (regeneration) or checks bit for bit (test) the values of one result.
        character(len=*), intent(in) :: key
            !! Unique key of the result (routine, test function, case).
        real(dp), intent(in) :: x(:)
            !! Values of the result.
        integer(int64) :: bits(size(x))
        integer :: i, j
        real(dp) :: expected(size(x)), scale, rel, maxrel

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
        if (all(bits == ref_bits(1:size(x), j))) return

        expected = transfer(ref_bits(1:size(x), j), expected)
        scale = max(maxval(abs(expected)), tiny(1.0_dp))
        maxrel = 0.0_dp
        print '(a)', ' FAIL: ' // trim(key) // ' differs from the golden value (bit-for-bit check)'
        print '(a)', '   comp   got                        expected                   |diff|/max|expected|'
        do i = 1, size(x)
            if (bits(i) == ref_bits(i, j)) cycle
            rel = abs(x(i) - expected(i))/scale
            maxrel = max(maxrel, rel)
            print '(3x, i4, 2x, es25.17, 2x, es25.17, 2x, es9.2)', i, x(i), expected(i), rel
        end do
        if (maxrel <= ACCEPTABLE_REL) then
            print '(a, es9.2, a)', '   max relative error ', maxrel, &
                ' <= 1e-14: acceptable only with explicit approval (testing.md)'
        else
            print '(a, es9.2, a)', '   max relative error ', maxrel, ' > 1e-14: behaviour change'
        end if
        flush(output_unit)
        if (code == 0) code = 1
        n_bad = n_bad + 1
    end subroutine golden_check

    subroutine golden_finish()
        !! Checks that every reference record was produced and stops with the first failure code.
        integer :: i
        if (writing) then
            close(wunit)
            return
        end if
        do i = 1, nref
            if (.not. ref_used(i)) then
                print '(a)', ' FAIL: reference record not produced by the test: ' // trim(ref_key(i))
                if (code == 0) code = 3
                n_bad = n_bad + 1
            end if
        end do
        print '(a, i0, a, i0, a, i0, a)', ' golden: ', n_checked, ' results checked, ', &
            nref, ' reference records, ', n_bad, ' failure(s)'
        select case (code)
        case (1)
            stop 1
        case (2)
            stop 2
        case (3)
            stop 3
        end select
    end subroutine golden_finish

end module muscle_test_golden
