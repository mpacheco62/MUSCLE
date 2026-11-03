! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_test_solver_golden
    !! Shared pieces of the golden (characterization) tests of the return mapping solvers
    !! (`tests/reference/solvers/`, Phase 7 step 1, D-018), independent of the API under test: the
    !! strain sequence, the record layout and the mutation hook.
    !!
    !! Strain sequence of every case (multiples of the yield strain `a` of the case along the fixed
    !! directions D1 and D2, see `seq_strain`):
    !!
    !! - `seq step 0..7`: zero strain (zero-strain case, first increment from the zero state), elastic
    !!   step (0.5 D1), first plastic step (1.5 D1), plastic loading (2.0, 2.5, 3.0 D1), elastic
    !!   unloading (2.6 D1) and one large non-proportional increment (2.6 D1 + 25 D2), with the state
    !!   committed after every step;
    !! - `first`: a single plastic increment 4.0 D1 from the zero state.
    !!
    !! Records of a step (`<solver> | <case> | <seq step i | first> | <what>`):
    !! `result` = [strain(6), stress(6), ep, strain_p(6), status, iters] (status and iters converted
    !! to real(dp), exact), and `tangent` / `tangent_num` = the 6x6 Voigt matrix (column major).
    !!
    !! Mutation hook: `<test> --mutate 1` perturbs by one ulp the parameter named by the test (see
    !! `mutate_name`) when its component is built; the ctest `*_golden_mutation` checks that the
    !! golden comparison detects it.
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym, ten_3D4O2sym
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: NSTEP, FIRST_FACTOR, seq_strain, first_strain, step_key
    public :: record_result, record_tangent, compare_tangents, n_tangent_notes
    public :: mutation_requested, scales_requested, mutate_name, bump

    integer, parameter :: NSTEP = 8
        !! Steps of the sequence (0..NSTEP-1).
    real(dp), parameter :: D1(6) = [1.0_dp, -0.3_dp, -0.2_dp, 0.25_dp, -0.1_dp, 0.15_dp]
        !! Main loading direction (volumetric part, three shear components).
    real(dp), parameter :: D2(6) = [-0.2_dp, 0.8_dp, -0.4_dp, -0.3_dp, 0.2_dp, 0.1_dp]
        !! Direction of the large non-proportional increment.
    real(dp), parameter :: F1(0:NSTEP-1) = [0.0_dp, 0.5_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, 2.6_dp, 2.6_dp]
        !! Factor of D1 at each step.
    real(dp), parameter :: F2(0:NSTEP-1) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 25.0_dp]
        !! Factor of D2 at each step.
    real(dp), parameter :: FIRST_FACTOR = 4.0_dp
        !! Factor of D1 of the single first increment.
    real(dp), parameter :: TANGENT_NOTE_TOL = 1.0e-5_dp
        !! Relative difference analytic vs numerical tangent above which a note is printed
        !! (information only: the golden records both tangents bit for bit).

    integer :: n_tangent_notes = 0
        !! Notes printed by `compare_tangents`.
    character(len=32) :: mutate_name = ''
        !! Parameter set perturbed by one ulp (empty: none). Set by the test with `--mutate`.

contains

    pure function seq_strain(step, a) result(strain)
        !! Strain of step `step` of the sequence for a yield strain `a`.
        integer, intent(in) :: step
        real(dp), intent(in) :: a
        type(ten_3D2Osym) :: strain
        strain%vals = (F1(step)*a)*D1 + (F2(step)*a)*D2
    end function seq_strain

    pure function first_strain(a) result(strain)
        !! Strain of the single first increment from the zero state.
        real(dp), intent(in) :: a
        type(ten_3D2Osym) :: strain
        strain%vals = (FIRST_FACTOR*a)*D1
    end function first_strain

    function step_key(prefix, step) result(key)
        !! `<prefix> | seq step <step>`, or `<prefix> | first` for step < 0.
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: step
        character(len=:), allocatable :: key
        character(len=16) :: buf
        if (step < 0) then
            key = prefix // ' | first'
        else
            write(buf, '(i0)') step
            key = prefix // ' | seq step ' // trim(buf)
        end if
    end function step_key

    subroutine record_result(key, strain, stress, ep, strain_p, status, iters)
        !! Golden record of the result of a step.
        character(len=*), intent(in) :: key
        type(ten_3D2Osym), intent(in) :: strain, stress, strain_p
        real(dp), intent(in) :: ep
        integer, intent(in) :: status, iters
        call golden_check(key // ' | result', [strain%vals, stress%vals, ep, strain_p%vals, &
            real(status, dp), real(iters, dp)])
    end subroutine record_result

    subroutine record_tangent(key, what, t)
        !! Golden record of a tangent (6x6 Voigt matrix, column major).
        character(len=*), intent(in) :: key, what
        type(ten_3D4O2sym), intent(in) :: t
        call golden_check(key // ' | ' // what, reshape(t%vals, [36]))
    end subroutine record_tangent

    subroutine compare_tangents(key, ta, tn)
        !! Prints a note when the analytic and numerical tangents differ beyond TANGENT_NOTE_TOL
        !! (relative to the largest numerical component). Information for the report, not a check.
        character(len=*), intent(in) :: key
        type(ten_3D4O2sym), intent(in) :: ta, tn
        real(dp) :: rel
        rel = maxval(abs(ta%vals - tn%vals))/max(maxval(abs(tn%vals)), tiny(1.0_dp))
        if (.not. (rel <= TANGENT_NOTE_TOL)) then
            n_tangent_notes = n_tangent_notes + 1
            print '(a, es10.3)', ' NOTE: ' // key // ': analytic vs numerical tangent, max rel diff ', rel
        end if
    end subroutine compare_tangents

    logical function mutation_requested()
        !! `--mutate` given as first argument.
        character(len=64) :: arg
        mutation_requested = .false.
        if (command_argument_count() >= 1) then
            call get_command_argument(1, arg)
            mutation_requested = trim(arg) == '--mutate'
        end if
    end function mutation_requested

    logical function scales_requested()
        !! `--scales` given as first argument (prints the yield strain of every case, used once to
        !! write the A_Y tables of the tests; never run by ctest).
        character(len=64) :: arg
        scales_requested = .false.
        if (command_argument_count() >= 1) then
            call get_command_argument(1, arg)
            scales_requested = trim(arg) == '--scales'
        end if
    end function scales_requested

    pure function bump(name, x) result(y)
        !! `x` perturbed by one ulp when `name` is the mutated parameter set.
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: x
        real(dp) :: y
        y = x
        if (len_trim(mutate_name) > 0 .and. name == mutate_name) y = nearest(x, 1.0_dp)
    end function bump

end module muscle_test_solver_golden
