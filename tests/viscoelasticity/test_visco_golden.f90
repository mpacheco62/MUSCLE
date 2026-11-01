! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module visco_golden_cases
    !! Cases of the golden test of the viscoelasticity family: every public and type-bound procedure
    !! of `Base_viscoelasticity` (through `Generalized_maxwell` and `Kelvin_voigt`, called through the
    !! polymorphic `class(Base_viscoelasticity)` interface), of the frozen algorithmic elasticity
    !! `Elasticity_frozen` (through `class(Base_elasticity)`, as the return mapping solvers use it,
    !! 3D and 2D, transient overloads included) and of the transactional state
    !! `Viscoelastic_history` (init, commit, rollback, pack/unpack of the FEA array), along
    !! deterministic multi-step load histories. After every step the coefficients, the frozen
    !! elasticity, its stresses and tangents and the WHOLE state (all `MAX_VISCO_BRANCHES` slots)
    !! are recorded.
    use, intrinsic :: iso_fortran_env, only : int64
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym, ten_2D2Osym, ten_3D4O3sym, ten_2D4O3sym
    use muscle_elasticity_base, only : Base_elasticity, Base_internal_elasticity_variables
    use muscle_visco_state, only : Viscoelastic_history, Viscoelastic_state, MAX_VISCO_BRANCHES
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Base_viscoelasticity, Visco_coefficients
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: golden_internal, record_coefs, record_state, record_frozen, record_to_string, &
        run_history, history_names, n_alt_bad, NHIST, NHSV

    integer, parameter :: NHIST = 7
        !! Number of load histories. (The eighth, dt = -0.01 then 0.1, was removed with D-017.13b: a
        !! negative time step is an error in every model; see `check_kernel_status`.)
    integer, parameter :: NHSV = 6*(MAX_VISCO_BRANCHES + 1)
        !! Length of the FEA history array (room for every branch).

    real(dp), parameter :: D1(6) = [1.0_dp, -0.4_dp, 0.2_dp, 0.6_dp, -0.2_dp, 0.1_dp]*1.0e-3_dp
        !! First strain direction (normal and shear, volumetric and deviatoric parts).
    real(dp), parameter :: D2(6) = [0.3_dp, 0.5_dp, -0.7_dp, -0.2_dp, 0.4_dp, 0.3_dp]*1.0e-3_dp
        !! Second strain direction (non-proportional paths).

    type, extends(Base_internal_elasticity_variables) :: golden_internal
        !! Concrete internal-variable container of the transient overloads (recorded after
        !! `stress_and_update` to detect any modification).
        real(dp) :: q(6) = [1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp, 5.0_dp, -6.0_dp]
    end type golden_internal

    integer :: n_alt_bad = 0
        !! Alternative paths (freeze_dt/update_dt, FEA wrapper, rollback, commit) that differ.

contains

    function history_names(ih) result(name)
        !! Name of the load history `ih`.
        integer, intent(in) :: ih
        character(len=:), allocatable :: name
        select case (ih)
        case (1); name = 'relaxation'
        case (2); name = 'constant rate'
        case (3); name = 'cycle with rollback'
        case (4); name = 'dt zero'
        case (5); name = 'dt large'
        case (6); name = 'dt tiny'
        case default; name = 'prestrain'
        end select
    end function history_names

    subroutine history_data(ih, eps, dts, rollback_at)
        !! Elastic strain at the end of each step, time step of each step, steps with a trial update
        !! that is rolled back before the real one.
        integer, intent(in) :: ih
        real(dp), allocatable, intent(out) :: eps(:, :), dts(:)
        logical, allocatable, intent(out) :: rollback_at(:)
        integer :: i
        select case (ih)
        case (1) ! step strain in a very short step, then hold with growing steps
            dts = [1.0e-12_dp, 0.01_dp, 0.1_dp, 1.0_dp, 10.0_dp, 100.0_dp]
            allocate(eps(6, size(dts)))
            do i = 1, size(dts)
                eps(:, i) = D1
            end do
        case (2) ! constant strain rate
            dts = [(0.05_dp, i = 1, 6)]
            allocate(eps(6, size(dts)))
            do i = 1, size(dts)
                eps(:, i) = 0.5_dp*real(i, dp)*D1
            end do
        case (3) ! non-proportional load-unload cycle through zero and compression
            dts = [(0.1_dp, i = 1, 8)]
            eps = reshape([D1, 2.0_dp*D1, D1 + D2, D2, 0.0_dp*D1, -D1, -2.0_dp*D1 + D2, 0.0_dp*D1], [6, 8])
        case (4) ! dt = 0 with and without a strain change
            dts = [0.1_dp, 0.0_dp, 0.0_dp, 0.1_dp]
            eps = reshape([D1, 2.0_dp*D1, 2.0_dp*D1, 2.0_dp*D1], [6, 4])
        case (5) ! dt much larger than every relaxation time
            dts = [1.0e6_dp, 1.0e6_dp, 1.0e6_dp]
            eps = reshape([D1, D1, D2], [6, 3])
        case (6) ! dt much smaller; the last two around the series switch x = 1e-5 of tau = 1
            dts = [1.0e-14_dp, 1.0e-14_dp, 9.99e-6_dp, 1.001e-5_dp]
            eps = reshape([D1, 1.001_dp*D1, D2, D2], [6, 4])
        case default ! (7) history initialised with a nonzero elastic strain (D2)
            dts = [0.2_dp, 0.2_dp, 0.2_dp]
            eps = reshape([D2, D1, D1], [6, 3])
        end select
        allocate(rollback_at(size(dts)))
        rollback_at = .false.
        if (ih == 3) rollback_at([3, 6]) = .true.
    end subroutine history_data

    function coefs_vals(c) result(v)
        !! Every component of the step coefficients (all branch slots).
        type(Visco_coefficients), intent(in) :: c
        real(dp) :: v(3 + 2*MAX_VISCO_BRANCHES)
        v = [c%dt, c%k_alg, c%g_alg, c%decay, c%weight]
    end function coefs_vals

    function state_vals(s) result(v)
        !! Every component of a state: n_branches, strain_e and all history slots.
        type(Viscoelastic_state), intent(in) :: s
        real(dp) :: v(7 + 6*MAX_VISCO_BRANCHES)
        v = [real(s%n_branches, dp), s%strain_e%vals, reshape(s%h, [6*MAX_VISCO_BRANCHES])]
    end function state_vals

    function state_active(s) result(v)
        !! Recorded part of a state: n_branches, strain_e and the active history tensors. The
        !! inactive slots are checked to be zero (`check_inactive`), not recorded.
        type(Viscoelastic_state), intent(in) :: s
        real(dp), allocatable :: v(:)
        integer :: n
        n = min(max(s%n_branches, 0), MAX_VISCO_BRANCHES)
        v = [real(s%n_branches, dp), s%strain_e%vals, reshape(s%h(:, 1:n), [6*n])]
    end function state_active

    subroutine check_inactive(key, s)
        !! The history slots beyond n_branches stay zero (no model writes them).
        character(len=*), intent(in) :: key
        type(Viscoelastic_state), intent(in) :: s
        integer :: n
        n = min(max(s%n_branches, 0), MAX_VISCO_BRANCHES)
        if (any(s%h(:, n+1:) /= 0.0_dp)) then
            print '(a)', ' FAIL: inactive history slots written: ' // key
            n_alt_bad = n_alt_bad + 1
        end if
    end subroutine check_inactive

    subroutine bits_equal(key, got, expected)
        !! Bitwise check of an alternative path against the reference path.
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: got(:), expected(:)
        if (size(got) /= size(expected)) then
            print '(a)', ' FAIL: alternative path differs (size): ' // key
            n_alt_bad = n_alt_bad + 1
        else if (any(transfer(got, 1_int64, size(got)) /= transfer(expected, 1_int64, size(expected)))) then
            print '(a)', ' FAIL: alternative path differs: ' // key
            n_alt_bad = n_alt_bad + 1
        end if
    end subroutine bits_equal

    subroutine record_coefs(key, c)
        character(len=*), intent(in) :: key
        type(Visco_coefficients), intent(in) :: c
        call golden_check(key, coefs_vals(c))
    end subroutine record_coefs

    subroutine record_state(key, h)
        !! Both states of a history.
        character(len=*), intent(in) :: key
        type(Viscoelastic_history), intent(in) :: h
        call golden_check(key // ' | state_n', state_active(h%state_n))
        call golden_check(key // ' | state_np1', state_active(h%state_np1))
        call check_inactive(key // ' state_n', h%state_n)
        call check_inactive(key // ' state_np1', h%state_np1)
    end subroutine record_state

    subroutine record_to_string(key, el)
        !! `to_string` (character codes) through the polymorphic base.
        character(len=*), intent(in) :: key
        class(Base_elasticity), intent(in) :: el
        character(len=275) :: str
        integer :: i
        str = el%to_string()
        call golden_check(key // ' | to_string', [(real(iachar(str(i:i)), dp), i = 1, len(str))])
    end subroutine record_to_string

    subroutine record_frozen(key, elas, eps, eps_n, dt, tangents, transient, string)
        !! Frozen elasticity at the elastic strain `eps` of the step: own bindings (bulk, shear,
        !! history_stress) and, through `class(Base_elasticity)`, stress and tangents in 3D and 2D
        !! (2D strain = in-plane components 11, 22, 33, 12 of `eps`); optionally the tangents, the
        !! transient overloads (increment eps - eps_n, `dt`) with `stress_and_update`, and
        !! `to_string`.
        character(len=*), intent(in) :: key
        type(Elasticity_frozen), intent(in) :: elas
        real(dp), intent(in) :: eps(6), eps_n(6), dt
        logical, intent(in) :: tangents, transient, string
        class(Base_elasticity), allocatable :: el
        type(ten_3D2Osym) :: e3, d3, s3, h3, t3
        type(ten_2D2Osym) :: e2, d2, s2, sd2, t2, td2
        type(ten_3D4O3sym) :: c3, cd3
        type(ten_2D4O3sym) :: c2
        type(golden_internal) :: internal
        real(dp), allocatable :: r(:)

        el = elas
        call e3%init(eps)
        call d3%init(eps - eps_n)
        call e2%init(eps(1:4))
        call d2%init(eps(1:4) - eps_n(1:4))
        h3 = elas%history_stress()
        s3 = el%stress(e3)
        t3 = el%stress_dev(e3)
        s2 = el%stress(e2)
        sd2 = el%stress_dev(e2)
        call golden_check(key // ' | frozen stress', [elas%bulk(), elas%shear(), h3%vals, s3%vals, &
            t3%vals, s2%vals, sd2%vals])
        c3 = el%dstress_dstrain(e3)
        cd3 = el%dstress_dstrain_dev(e3)
        c2 = el%dstress_dstrain(e2)
        if (tangents) call golden_check(key // ' | frozen tangents', [c3%vals, cd3%vals, c2%vals])
        if (string) call record_to_string(key // ' | frozen', el)
        if (.not. transient) return
        s3 = el%stress(e3, d3, dt, internal)
        t3 = el%stress_dev(e3, d3, dt, internal)
        r = [s3%vals, t3%vals]
        c3 = el%dstress_dstrain(e3, d3, dt, internal)
        cd3 = el%dstress_dstrain_dev(e3, d3, dt, internal)
        r = [r, c3%vals, cd3%vals]
        t2 = el%stress(e2, d2, dt, internal)
        td2 = el%stress_dev(e2, d2, dt, internal)
        c2 = el%dstress_dstrain(e2, d2, dt, internal)
        r = [r, t2%vals, td2%vals, c2%vals]
        call el%stress_and_update(s3, e3, d3, dt, internal)
        r = [r, s3%vals, internal%q]
        call golden_check(key // ' | frozen transient', r)
    end subroutine record_frozen

    subroutine run_history(key, vm, ih)
        !! Runs the load history `ih` on the model `vm`, as a solver does (coefficients once per step,
        !! freeze from state_n, update with the converged elastic strain, commit), recording after
        !! every step. In parallel: the convenience path (freeze_dt/update_dt) and the FEA wrapper
        !! path (init_history + unpack_from_fea + freeze + update + pack_to_fea on a new local
        !! history every step), both checked bit for bit (the FEA array also against its layout and
        !! recorded at the last step).
        character(len=*), intent(in) :: key
        class(Base_viscoelasticity), intent(in) :: vm
        integer, intent(in) :: ih
        real(dp), allocatable :: eps(:, :), dts(:)
        logical, allocatable :: rollback_at(:)
        type(Viscoelastic_history) :: hist, hist_dt, hist_fea
        type(Visco_coefficients) :: coefs
        type(Elasticity_frozen) :: elas, elas_dt, elas_fea
        type(ten_3D2Osym) :: e3, trial, s_a, s_b
        type(Viscoelastic_state) :: state_before
        real(dp) :: hsv(NHSV), expected_hsv(NHSV), eps_n(6)
        integer :: i, n
        character(len=:), allocatable :: kk, ks
        character(len=160) :: buf

        kk = key // ' | ' // history_names(ih)
        call history_data(ih, eps, dts, rollback_at)
        hsv = 0.0_dp
        if (ih == 7) then
            call e3%init(eps(:, 1))
            call vm%init_history(hist, e3)
            hsv(1:6) = eps(:, 1)
        else
            call vm%init_history(hist)
        end if
        hist_dt = hist
        call record_state(kk // ' | step 0', hist)

        do i = 1, size(dts)
            write(buf, '(a, " | step ", i0)') kk, i
            ks = trim(buf)
            eps_n = hist%state_n%strain_e%vals
            call e3%init(eps(:, i))
            coefs = vm%coefficients(dts(i))
            call golden_check(ks // ' | dstrain_e, dt, coefficients', [eps(:, i) - eps_n, dts(i), &
                coefs_vals(coefs)])

            elas = vm%freeze(coefs, hist)
            call record_frozen(ks, elas, eps(:, i), eps_n, dts(i), tangents=(i == 1 .or. i == size(dts)), &
                transient=(i == 1), string=(i == size(dts) .and. ih == 3))
            ! Convenience path, own state.
            elas_dt = vm%freeze_dt(dts(i), hist_dt)
            s_a = elas%stress(e3)
            s_b = elas_dt%stress(e3)
            call bits_equal(ks // ' freeze_dt stress', s_b%vals, s_a%vals)
            ! FEA wrapper path: everything rebuilt from the array.
            call vm%init_history(hist_fea)
            call hist_fea%unpack_from_fea(hsv)
            elas_fea = vm%freeze(coefs, hist_fea)
            s_b = elas_fea%stress(e3)
            call bits_equal(ks // ' FEA wrapper stress', s_b%vals, s_a%vals)

            if (rollback_at(i)) then
                ! A failed global attempt: trial update, then rollback.
                call trial%init(3.0_dp*eps(:, i) + D2)
                call vm%update(coefs, trial, hist)
                call golden_check(ks // ' | trial state_np1', state_active(hist%state_np1))
                call hist%rollback()
                call bits_equal(ks // ' rollback state_np1', state_vals(hist%state_np1), &
                    state_vals(hist%state_n))
            end if

            state_before = hist%state_n
            call vm%update(coefs, e3, hist)
            call golden_check(ks // ' | state_np1 after update', state_active(hist%state_np1))
            call check_inactive(ks // ' state_np1 after update', hist%state_np1)
            call bits_equal(ks // ' update keeps state_n', state_vals(hist%state_n), state_vals(state_before))
            call vm%update_dt(dts(i), e3, hist_dt)
            call bits_equal(ks // ' update_dt state_np1', state_vals(hist_dt%state_np1), &
                state_vals(hist%state_np1))
            call vm%update(coefs, e3, hist_fea)
            call hist_fea%pack_to_fea(hsv)
            ! FEA layout: strain_e, then the history tensor of each active branch; the rest of the
            ! array untouched (zero here, as initialised).
            n = hist%state_np1%n_branches
            expected_hsv = 0.0_dp
            expected_hsv(1:6) = hist%state_np1%strain_e%vals
            expected_hsv(7:6*(n + 1)) = reshape(hist%state_np1%h(:, 1:n), [6*n])
            call bits_equal(ks // ' FEA array layout', hsv, expected_hsv)
            call bits_equal(ks // ' FEA wrapper state_np1', state_vals(hist_fea%state_np1), &
                state_vals(hist%state_np1))
            if (i == size(dts)) call golden_check(ks // ' | FEA array', hsv)

            call hist%commit()
            call hist_dt%commit()
            call bits_equal(ks // ' commit state_n', state_vals(hist%state_n), state_vals(hist%state_np1))
        end do
    end subroutine run_history

end module visco_golden_cases

program test_visco_golden
    !! Golden (characterization) test of the viscoelasticity family in double precision, bit for
    !! bit against tests/reference/viscoelasticity/viscoelasticity_dp.txt, captured from the
    !! original OO code (regenerate.sh there). Uses the public OO API as its users do
    !! (`set_parameters`, `class(Base_viscoelasticity)`/`class(Base_elasticity)` allocatables with
    !! intrinsic assignment, `Viscoelastic_history`, the FEA array).
    !!
    !! Exit codes: 1-4 from `golden_finish` (muscle_test_golden); 5: an alternative path
    !! (freeze_dt/update_dt, FEA wrapper, rollback, commit) differs from the reference path, or a
    !! status of the kernel API is not the expected one.
    !!
    !! D-017.13b (approved change): a negative time step is an error in every model; the OO adapters
    !! stop (`error stop`), so the Maxwell records with dt < 0 (the former load history 8 and the
    !! dt = -0.01 entry of the coefficient tables) can no longer be produced in process. They were
    !! removed from the reference and replaced by checks of the kernel status (`check_kernel_status`,
    !! not recorded): dt < 0 gives STATUS_INVALID_TIME_STEP in both models, dt = 0 is valid in Maxwell
    !! and an error in Kelvin-Voigt.
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym
    use muscle_visco_state, only : Viscoelastic_history, MAX_VISCO_BRANCHES
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Base_viscoelasticity, Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    use muscle_visco_kelvin, only : Kelvin_voigt
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use visco_golden_cases, only : record_coefs, record_state, record_frozen, run_history, &
        n_alt_bad, NHIST
    use muscle_status, only : STATUS_OK, STATUS_INVALID_TIME_STEP
    use muscle_visco_types, only : visco_coefs_dp
    use muscle_visco_maxwell, only : visco_maxwell_dp, visco_maxwell_init, visco_step_coefficients_maxwell_dp
    use muscle_visco_kelvin, only : visco_kelvin_dp, visco_kelvin_init, visco_step_coefficients_kelvin_dp
    implicit none

    real(dp), parameter :: DT_LIST(10) = [0.0_dp, 1.0e-14_dp, 9.99e-6_dp, 1.0e-5_dp, 1.001e-5_dp, &
        1.0e-3_dp, 1.0_dp, 1.0e3_dp, 1.0e12_dp, -0.01_dp]
        !! Time steps of the coefficient tables (Maxwell; Kelvin only the positive ones).
    real(dp), parameter :: T_LIST(8) = [0.0_dp, 1.0e-6_dp, 0.1_dp, 1.0_dp, 10.0_dp, 1.0e6_dp, &
        1.0e300_dp, -1.0_dp]
        !! Times of the relaxation modulus.
    real(dp), parameter :: KV(4, 5) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &           ! placeholder of the unset object (not used)
        2000.0_dp, 100.0_dp, 50.0_dp, 20.0_dp, &    ! K, G, zeta, mu_v
        2000.0_dp, 100.0_dp, 0.0_dp, 0.0_dp, &      ! no dashpot: elastic
        2000.0_dp, 100.0_dp, 1000.0_dp, 0.0_dp, &   ! volumetric dashpot only
        0.0_dp, 0.0_dp, 30.0_dp, 5.0_dp], [4, 5])   ! no spring: pure dashpot
        !! Kelvin-Voigt parameter sets (bulk, shear, bulk_viscosity, shear_viscosity).

    type(Generalized_maxwell) :: gm, gm_unset
    type(Kelvin_voigt) :: kv_model, kv_unset
    class(Base_viscoelasticity), allocatable :: vm
    type(Viscoelastic_history) :: hist_default
    type(Visco_coefficients) :: coefs_default
    type(Elasticity_frozen) :: frozen
    type(ten_3D2Osym) :: hs
    real(dp), allocatable :: r(:)
    real(dp) :: tau10(MAX_VISCO_BRANCHES), g10(MAX_VISCO_BRANCHES)
    character(len=48) :: key
    integer :: im, ih, i

    call golden_init('viscoelasticity/viscoelasticity_dp.txt')

    ! Default-initialised containers and frozen elasticity.
    call record_state('defaults | history', hist_default)
    call record_coefs('defaults | coefficients', coefs_default)
    call record_frozen('defaults | frozen', frozen, D_1(), [(0.0_dp, i = 1, 6)], 0.1_dp, .true., .true., .true.)
    ! Elasticity_frozen%set directly (as freeze does), then every binding.
    call hs%init([10.0_dp, -20.0_dp, 30.0_dp, -4.0_dp, 5.0_dp, -6.0_dp])
    call frozen%set(1500.0_dp, 250.0_dp, hs)
    call record_frozen('set | frozen', frozen, D_1(), [(0.0_dp, i = 1, 6)], 0.1_dp, .true., .true., .true.)
    call frozen%set(0.0_dp, 0.0_dp, hs)
    call record_frozen('set zero moduli | frozen', frozen, D_1(), [(0.0_dp, i = 1, 6)], 0.1_dp, .true., .true., .true.)

    ! Generalized Maxwell: parameter sets.
    do i = 1, MAX_VISCO_BRANCHES
        tau10(i) = 10.0_dp**(i - 5)
        g10(i) = 10.0_dp*real(i, dp)
    end do
    do im = 0, 9
        select case (im)
        case (0) ! default object (set_parameters never called)
            gm = gm_unset
        case (1) ! no branches: elastic
            call gm%set_parameters(bulk=2000.0_dp, shear_inf=100.0_dp, shear=[real(dp) ::], tau=[real(dp) ::])
        case (2) ! standard linear solid
            call gm%set_parameters(bulk=2000.0_dp, shear_inf=100.0_dp, shear=[300.0_dp], tau=[1.0_dp])
        case (3) ! single Maxwell element (fluid), G_inf = 0
            call gm%set_parameters(bulk=2000.0_dp, shear_inf=0.0_dp, shear=[500.0_dp], tau=[0.5_dp])
        case (4) ! two branches
            call gm%set_parameters(bulk=2000.0_dp, shear_inf=100.0_dp, shear=[200.0_dp, 300.0_dp], &
                tau=[0.1_dp, 1.0_dp])
        case (5) ! three branches with very different relaxation times, G_inf = 0
            call gm%set_parameters(bulk=2000.0_dp, shear_inf=0.0_dp, shear=[100.0_dp, 1000.0_dp, 10.0_dp], &
                tau=[1.0e-6_dp, 1.0_dp, 1.0e6_dp])
        case (6) ! MAX_VISCO_BRANCHES branches, tau = 1e-4 .. 1e5
            call gm%set_parameters(bulk=2000.0_dp, shear_inf=50.0_dp, shear=g10, tau=tau10)
        case (7) ! normalized (Abaqus) Prony series
            call gm%set_parameters_normalized(young=1000.0_dp, poisson=0.3_dp, g=[0.3_dp, 0.5_dp], &
                tau=[0.2_dp, 5.0_dp])
        case (8) ! normalized, sum(g) = 1: Maxwell fluid
            call gm%set_parameters_normalized(young=1000.0_dp, poisson=0.3_dp, g=[0.25_dp, 0.75_dp], &
                tau=[1.0_dp, 10.0_dp])
        case default ! re-initialisation with fewer branches: the unused slots must be reset
            call gm%set_parameters(bulk=2000.0_dp, shear_inf=50.0_dp, shear=g10, tau=tau10)
            call gm%set_parameters(bulk=3000.0_dp, shear_inf=100.0_dp, shear=[300.0_dp], tau=[2.0_dp])
        end select
        write(key, '(a, i0)') 'maxwell ', im
        vm = gm
        r = [real(vm%n_branches(), dp), gm%shear_instantaneous()]
        do i = 1, size(T_LIST)
            r = [r, gm%relaxation_shear(T_LIST(i))]
        end do
        call golden_check(trim(key) // ' | n_branches, G0, G(t)', r)
        r = [real(dp) ::]
        do i = 1, size(DT_LIST)
            if (DT_LIST(i) < 0.0_dp) cycle ! an error since D-017.13b (check_kernel_status)
            r = [r, coefs_table(vm%coefficients(DT_LIST(i)))]
        end do
        call golden_check(trim(key) // ' | coefficients table', r)
        do ih = 1, NHIST
            call run_history(trim(key), vm, ih)
        end do
        deallocate(vm)
    end do
    call check_kernel_status()

    ! Kelvin-Voigt: parameter sets (dt = 0 stops with an error: history 4 not run).
    do im = 1, size(KV, 2)
        if (im == 1) then
            kv_model = kv_unset
        else
            call kv_model%set_parameters(bulk=KV(1, im), shear=KV(2, im), bulk_viscosity=KV(3, im), &
                shear_viscosity=KV(4, im))
        end if
        write(key, '(a, i0)') 'kelvin ', im - 1
        vm = kv_model
        r = [real(vm%n_branches(), dp)]
        do i = 1, size(DT_LIST)
            if (DT_LIST(i) > 0.0_dp) r = [r, coefs_table(vm%coefficients(DT_LIST(i)))]
        end do
        call golden_check(trim(key) // ' | n_branches, coefficients table', r)
        do ih = 1, NHIST
            if (ih == 4) cycle ! dt = 0
            call run_history(trim(key), vm, ih)
        end do
        deallocate(vm)
    end do
    ! Re-initialisation of a Kelvin-Voigt object.
    call kv_model%set_parameters(bulk=1.0_dp, shear=2.0_dp, bulk_viscosity=3.0_dp, shear_viscosity=4.0_dp)
    call kv_model%set_parameters(bulk=KV(1, 2), shear=KV(2, 2), bulk_viscosity=KV(3, 2), &
        shear_viscosity=KV(4, 2))
    vm = kv_model
    call run_history('kelvin reinit', vm, 3)

    call golden_finish()
    if (n_alt_bad > 0) stop 5

contains

    subroutine check_kernel_status()
        !! D-017.13b: status of `visco_step_coefficients` with dt < 0 (the time steps of the former
        !! load history 8 and of the coefficient tables) and dt = 0, for the default parameters and a parameter set of
        !! each model; on an error the coefficients keep their default values.
        type(visco_maxwell_dp) :: pm(2)
        type(visco_kelvin_dp) :: pk(2)
        type(visco_coefs_dp) :: c, c_default
        integer :: status, j

        call visco_maxwell_init(pm(2), bulk=2000.0_dp, shear_inf=100.0_dp, shear=[300.0_dp], tau=[1.0_dp], &
            status=status)
        call expect('maxwell init', status, STATUS_OK)
        call visco_kelvin_init(pk(2), bulk=KV(1, 2), shear=KV(2, 2), bulk_viscosity=KV(3, 2), &
            shear_viscosity=KV(4, 2), status=status)
        call expect('kelvin init', status, STATUS_OK)
        do j = 1, 2
            call visco_step_coefficients_maxwell_dp(pm(j), -0.01_dp, c, status)
            call expect('maxwell dt < 0', status, STATUS_INVALID_TIME_STEP)
            if (any(coefs_table(c) /= coefs_table(c_default))) call expect('maxwell dt < 0 coefs', 1, 0)
            call visco_step_coefficients_maxwell_dp(pm(j), 0.0_dp, c, status)
            call expect('maxwell dt = 0', status, STATUS_OK)
            call visco_step_coefficients_kelvin_dp(pk(j), -0.01_dp, c, status)
            call expect('kelvin dt < 0', status, STATUS_INVALID_TIME_STEP)
            if (any(coefs_table(c) /= coefs_table(c_default))) call expect('kelvin dt < 0 coefs', 1, 0)
            call visco_step_coefficients_kelvin_dp(pk(j), 0.0_dp, c, status)
            call expect('kelvin dt = 0', status, STATUS_INVALID_TIME_STEP)
        end do
    end subroutine check_kernel_status

    subroutine expect(what, got, expected)
        !! Counts a kernel status that is not the expected one as a failed alternative path.
        character(len=*), intent(in) :: what
        integer, intent(in) :: got, expected
        if (got /= expected) then
            print '(a, i0, a, i0)', ' FAIL: kernel status of ' // what // ': ', got, ', expected ', expected
            n_alt_bad = n_alt_bad + 1
        end if
    end subroutine expect

    function D_1() result(v)
        !! First strain direction of the load histories.
        real(dp) :: v(6)
        v = [1.0_dp, -0.4_dp, 0.2_dp, 0.6_dp, -0.2_dp, 0.1_dp]*1.0e-3_dp
    end function D_1

    function coefs_table(c) result(v)
        type(Visco_coefficients), intent(in) :: c
        real(dp) :: v(3 + 2*MAX_VISCO_BRANCHES)
        v = [c%dt, c%k_alg, c%g_alg, c%decay, c%weight]
    end function coefs_table

end program test_visco_golden
