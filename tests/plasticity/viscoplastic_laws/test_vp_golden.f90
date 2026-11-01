! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module vp_golden_cases
    !! Cases of the golden test of the viscoplastic laws: every public type-bound procedure of
    !! every law through the polymorphic bases used by the solvers (`Base_viscoplastic_law`:
    !! `flow_stress`, `dstress_dep`; `Base_overstress_law`: `flow_rate`, `dflow_dsigma`,
    !! `dflow_dep`), plus the inherited finite-difference defaults of both bases, evaluated on
    !! every law through a test-local extension that only provides the deferred procedure.
    use, intrinsic :: iso_fortran_env, only : int64
    use muscle_kinds, only : dp
    use muscle_hard_base, only : Base_hardening_laws
    use muscle_vp_base, only : Base_viscoplastic_law
    use muscle_vp_overstress_base, only : Base_overstress_law
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: vp_sweep, os_sweep, bits_equal, n_direct_bad, DT_VALUES

    real(dp), parameter :: EP_SWEEP(10) = [0.0_dp, 1.0e-10_dp, 1.0e-7_dp, 1.0e-5_dp, 1.0e-3_dp, &
        1.0e-2_dp, 5.0e-2_dp, 0.2_dp, 1.0_dp, 5.0_dp]
        !! Equivalent plastic strains: zero (first increment), around the 1e-10 and 1e-7 floors of
        !! NNL and VA, small, intermediate and large.
    real(dp), parameter :: EPD_SWEEP(17) = [-1.0e-3_dp, 0.0_dp, 1.0e-12_dp, 1.0e-10_dp, 1.0e-9_dp, &
        1.0e-7_dp, 1.0e-6_dp, 1.0e-5_dp, 1.0e-4_dp, 1.0e-3_dp, 1.0e-2_dp, 0.1_dp, 1.0_dp, 10.0_dp, &
        1.0e2_dp, 1.0e3_dp, 1.0e4_dp]
        !! Equivalent plastic strain rates (1/s): negative and zero (first increment), the 1e-12,
        !! 1e-7 and 1e-6 floors of the laws, quasi-static to dynamic rates. Each law adds its own
        !! thresholds (reference/limit rates).
    real(dp), parameter :: SIG_FACTORS(11) = [0.0_dp, 0.5_dp, 0.999_dp, 1.0_dp, 1.0_dp + 1.0e-9_dp, &
        1.001_dp, 1.01_dp, 1.1_dp, 1.5_dp, 2.0_dp, 5.0_dp]
        !! Equivalent stress as a multiple of the static yield stress Y(ep): elastic, exactly at
        !! the threshold, just above it and well above it.
    real(dp), parameter :: DT_VALUES(2) = [1.0e-3_dp, 1.0_dp]
        !! Time steps of `dstress_dep` (two values separate d/dep from d/depd).

    type, extends(Base_viscoplastic_law) :: Numeric_vp
        !! Only `flow_stress` (delegated): `dstress_dep` is the inherited finite difference.
        class(Base_viscoplastic_law), allocatable :: inner
    contains
        procedure :: flow_stress => numeric_vp_flow
    end type Numeric_vp

    type, extends(Base_overstress_law) :: Numeric_os
        !! Only `flow_rate` (delegated): `dflow_dsigma`/`dflow_dep` are the inherited differences.
        class(Base_overstress_law), allocatable :: inner
    contains
        procedure :: flow_rate => numeric_os_flow
    end type Numeric_os

    integer :: n_direct_bad = 0
        !! Direct (non-polymorphic) calls that differ from the polymorphic ones.

contains

    pure function numeric_vp_flow(self, ep, epd) result(res)
        class(Numeric_vp), intent(in) :: self
        real(dp), intent(in) :: ep, epd
        real(dp) :: res
        res = self%inner%flow_stress(ep, epd)
    end function numeric_vp_flow

    pure function numeric_os_flow(self, sigma_eq, ep) result(res)
        class(Numeric_os), intent(in) :: self
        real(dp), intent(in) :: sigma_eq, ep
        real(dp) :: res
        res = self%inner%flow_rate(sigma_eq, ep)
    end function numeric_os_flow

    subroutine vp_sweep(key, v, epd_extra, ep_g, epd_g, flow, dflow)
        !! Golden records on the grid EP_SWEEP x [EPD_SWEEP, epd_extra] (ep fastest) of the
        !! inputs, `flow_stress`, `dstress_dep` at every DT_VALUES, and the inherited numeric
        !! `dstress_dep` of the base, all through the polymorphic base interface.
        character(len=*), intent(in) :: key
            !! Key prefix (law, parameter set, hardening law).
        class(Base_viscoplastic_law), intent(in) :: v
            !! Viscoplastic law.
        real(dp), intent(in) :: epd_extra(:)
            !! Law-specific rates (thresholds).
        real(dp), allocatable, intent(out) :: ep_g(:), epd_g(:)
            !! Grid.
        real(dp), allocatable, intent(out) :: flow(:), dflow(:)
            !! `flow_stress` and `dstress_dep` (first time step) on the grid.
        type(Numeric_vp) :: w
        real(dp), allocatable :: epd(:), d(:)
        character(len=16) :: dts
        integer :: i, j, n, it

        epd = [EPD_SWEEP, epd_extra]
        n = size(EP_SWEEP)*size(epd)
        allocate(ep_g(n), epd_g(n), flow(n), dflow(n), d(n))
        do j = 1, size(epd)
            do i = 1, size(EP_SWEEP)
                ep_g(i + (j - 1)*size(EP_SWEEP)) = EP_SWEEP(i)
                epd_g(i + (j - 1)*size(EP_SWEEP)) = epd(j)
            end do
        end do
        call golden_check(key // ' | ep', ep_g)
        call golden_check(key // ' | epd', epd_g)
        do i = 1, n
            flow(i) = v%flow_stress(ep_g(i), epd_g(i))
        end do
        call golden_check(key // ' | flow_stress', flow)
        allocate(w%inner, source=v)
        do it = 1, size(DT_VALUES)
            write(dts, '(es8.1)') DT_VALUES(it)
            do i = 1, n
                d(i) = v%dstress_dep(ep_g(i), epd_g(i), DT_VALUES(it))
            end do
            if (it == 1) dflow = d
            call golden_check(key // ' | dstress_dep dt=' // trim(adjustl(dts)), d)
            do i = 1, n
                d(i) = w%dstress_dep(ep_g(i), epd_g(i), DT_VALUES(it))
            end do
            call golden_check(key // ' | base numeric dstress_dep dt=' // trim(adjustl(dts)), d)
        end do
    end subroutine vp_sweep

    subroutine os_sweep(key, o, h, ep_g, sig_g, g)
        !! Golden records on the grid EP_SWEEP x SIG_FACTORS*Y(ep) (ep fastest) of the inputs,
        !! `flow_rate`, `dflow_dsigma`, `dflow_dep` and the inherited numeric derivatives of the
        !! base, through the polymorphic base interface.
        character(len=*), intent(in) :: key
            !! Key prefix (law, parameter set, hardening law).
        class(Base_overstress_law), intent(in) :: o
            !! Overstress law.
        class(Base_hardening_laws), intent(in) :: h
            !! The hardening law composed by `o` (to place the stresses around Y(ep)).
        real(dp), allocatable, intent(out) :: ep_g(:), sig_g(:)
            !! Grid.
        real(dp), allocatable, intent(out) :: g(:)
            !! `flow_rate` on the grid.
        type(Numeric_os) :: w
        real(dp), allocatable :: r(:)
        integer :: i, j, n

        n = size(EP_SWEEP)*size(SIG_FACTORS)
        allocate(ep_g(n), sig_g(n), g(n), r(n))
        do j = 1, size(SIG_FACTORS)
            do i = 1, size(EP_SWEEP)
                ep_g(i + (j - 1)*size(EP_SWEEP)) = EP_SWEEP(i)
                sig_g(i + (j - 1)*size(EP_SWEEP)) = SIG_FACTORS(j)*h%stress(EP_SWEEP(i))
            end do
        end do
        call golden_check(key // ' | ep', ep_g)
        call golden_check(key // ' | sigma_eq', sig_g)
        do i = 1, n
            g(i) = o%flow_rate(sig_g(i), ep_g(i))
        end do
        call golden_check(key // ' | flow_rate', g)
        do i = 1, n
            r(i) = o%dflow_dsigma(sig_g(i), ep_g(i))
        end do
        call golden_check(key // ' | dflow_dsigma', r)
        do i = 1, n
            r(i) = o%dflow_dep(sig_g(i), ep_g(i))
        end do
        call golden_check(key // ' | dflow_dep', r)
        allocate(w%inner, source=o)
        do i = 1, n
            r(i) = w%dflow_dsigma(sig_g(i), ep_g(i))
        end do
        call golden_check(key // ' | base numeric dflow_dsigma', r)
        do i = 1, n
            r(i) = w%dflow_dep(sig_g(i), ep_g(i))
        end do
        call golden_check(key // ' | base numeric dflow_dep', r)
    end subroutine os_sweep

    subroutine bits_equal(key, got, expected)
        !! Bitwise check of direct (concrete type) calls against the polymorphic ones.
        character(len=*), intent(in) :: key
            !! Description of the result.
        real(dp), intent(in) :: got(:), expected(:)
            !! Direct and polymorphic results.
        if (any(transfer(got, 1_int64, size(got)) /= transfer(expected, 1_int64, size(expected)))) then
            print '(a)', ' FAIL: direct call differs from polymorphic: ' // key
            n_direct_bad = n_direct_bad + 1
        end if
    end subroutine bits_equal

end module vp_golden_cases

program test_vp_golden
    !! Golden (characterization) test of the viscoplastic laws in double precision, bit for bit
    !! against tests/reference/plasticity/viscoplastic_laws/vp_laws_dp.txt, captured from the
    !! original OO code (regenerate.sh there). Uses the public OO API as today's solvers do
    !! (keyword structure constructors, hardening law assigned to the `hard_law` component,
    !! polymorphic `class(Base_viscoplastic_law)` / `class(Base_overstress_law)` calls), so it
    !! stays valid with the OO adapters of Phase 6 (D-017.1). Includes ep = 0 and epd = 0 (first
    !! plastic increment, whose numeric derivatives changed in Phase 5, D-016.4).
    !!
    !! Exit codes: 1-4 from `golden_finish` (muscle_test_golden); 5: a direct call on the
    !! concrete type differs from the polymorphic call.
    use muscle_kinds, only : dp
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_hard_swift, only : Swift_hardening
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_jc, only : JC_viscoplastic
    use muscle_vp_jcfull, only : JohnsonCook_full_visco_hardening
    use muscle_vp_mrk, only : MRK_viscoplastic
    use muscle_vp_rk, only : RK_viscoplastic
    use muscle_vp_nnl, only : NNL_viscoplastic
    use muscle_vp_va, only : VA_viscoplastic
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_vp_chaboche, only : Chaboche_overstress
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use vp_golden_cases, only : vp_sweep, os_sweep, bits_equal, n_direct_bad, DT_VALUES
    implicit none

    character(len=*), parameter :: HNAME(3) = [character(len=8) :: 'bilinear', 'swift', 'vocemod']
        !! Composed hardening laws (static yield stress Y(ep), MPa).

    real(dp), parameter :: JCP(3, 3) = reshape([ &
        0.00946101_dp, 0.32_dp, 1.0e-7_dp, &   ! fit of test_muscle_vp_models, default epmin
        0.014_dp, 1.0_dp, 1.0e-7_dp, &         ! classic JC steel, reference rate 1/s
        0.05_dp, 1.0e-3_dp, 1.0e-5_dp], [3, 3])  ! strong rate sensitivity, larger floor
        !! JC: C, epdmax (reference rate), epmin (rate floor).
    real(dp), parameter :: JCFP(6, 6) = reshape([ &
        0.014_dp, 1.0_dp, 1.03_dp, 293.15_dp, 1793.0_dp, 293.15_dp, &  ! T = Troom: Tstar = 0
        0.014_dp, 1.0_dp, 1.03_dp, 293.15_dp, 1793.0_dp, 700.0_dp, &   ! intermediate Tstar
        0.014_dp, 1.0_dp, 1.03_dp, 293.15_dp, 1793.0_dp, 2000.0_dp, &  ! T > Tmelt: Tstar = 1
        0.022_dp, 1.0e-3_dp, 0.8_dp, 293.15_dp, 1356.0_dp, 200.0_dp, & ! T < Troom: Tstar = 0
        0.014_dp, 1.0_dp, 0.0_dp, 293.15_dp, 1793.0_dp, 700.0_dp, &    ! m = 0 branch
        0.014_dp, 1.0_dp, 1.0_dp, 293.15_dp, 200.0_dp, 700.0_dp], [6, 6])  ! Tmelt < Troom
        !! JC full: C, epdmax, m, Troom, Tmelt, T_current.

    type(Bilinear_hardening) :: bil
    type(Swift_hardening) :: swi
    type(Voce_modified_hardening) :: voc
    type(JC_viscoplastic) :: jc
    type(JohnsonCook_full_visco_hardening) :: jcf
    type(MRK_viscoplastic) :: mrk
    type(RK_viscoplastic) :: rk
    type(NNL_viscoplastic) :: nnl
    type(VA_viscoplastic) :: va
    type(Perzyna_overstress) :: per
    type(Chaboche_overstress) :: cha
    real(dp), allocatable :: ep(:), epd(:), flow(:), dflow(:), sig(:)
    character(len=64) :: key
    integer :: i, ih, n, i_

    call golden_init('plasticity/viscoplastic_laws/vp_laws_dp.txt')

    bil = Bilinear_hardening(y0=250.0_dp, K=1000.0_dp)
    swi = Swift_hardening(k=500.0_dp, e0=1.0e-3_dp, n=0.2_dp)
    voc = Voce_modified_hardening(Sy=272.075784_dp, k=0.04230655_dp, q=297.603923_dp, n=1.48375781_dp)

    ! ---- Johnson-Cook (rate factor only), composed with each hardening law.
    do i = 1, size(JCP, 2)
        jc = JC_viscoplastic(C=JCP(1, i), epdmax=JCP(2, i))
        if (i == 3) jc%epmin = JCP(3, i)
        write(key, '(a, i0)') 'jc | set ', i
        call golden_check(trim(key) // ' | constructor', [jc%C, jc%epdmax, jc%epmin])
        do ih = 1, 3
            call set_hard(ih, jc_law=jc)
            call vp_sweep(trim(key) // ' | ' // trim(HNAME(ih)), jc, &
                [jc%epmin, jc%epmin*(1.0_dp + 1.0e-9_dp), jc%epdmax], ep, epd, flow, dflow)
            n = size(ep)
            call bits_equal(trim(key) // ' flow_stress', [(jc%flow_stress(ep(i_), epd(i_)), i_ = 1, n)], flow)
            call bits_equal(trim(key) // ' dstress_dep', &
                [(jc%dstress_dep(ep(i_), epd(i_), DT_VALUES(1)), i_ = 1, n)], dflow)
        end do
    end do

    ! ---- Johnson-Cook full (rate and temperature; inherited numeric dstress_dep).
    do i = 1, size(JCFP, 2)
        jcf = JohnsonCook_full_visco_hardening(C=JCFP(1, i), epdmax=JCFP(2, i), m=JCFP(3, i), &
            Troom=JCFP(4, i), Tmelt=JCFP(5, i), T_current=JCFP(6, i))
        write(key, '(a, i0)') 'jcfull | set ', i
        call golden_check(trim(key) // ' | constructor', &
            [jcf%C, jcf%epdmax, jcf%m, jcf%Troom, jcf%Tmelt, jcf%T_current])
        do ih = 1, 3
            if (ih == 1 .and. i /= 2) cycle   ! bilinear only with the intermediate temperature
            if (ih == 3 .and. i /= 2) cycle
            call set_hard(ih, jcf_law=jcf)
            call vp_sweep(trim(key) // ' | ' // trim(HNAME(ih)), jcf, [jcf%epdmax, tiny(1.0_dp)], &
                ep, epd, flow, dflow)
            n = size(ep)
            call bits_equal(trim(key) // ' flow_stress', [(jcf%flow_stress(ep(i_), epd(i_)), i_ = 1, n)], flow)
            call bits_equal(trim(key) // ' dstress_dep', &
                [(jcf%dstress_dep(ep(i_), epd(i_), DT_VALUES(1)), i_ = 1, n)], dflow)
        end do
    end do
    ! Default components (only the hardening law given).
    jcf = JohnsonCook_full_visco_hardening()
    call golden_check('jcfull | defaults | constructor', &
        [jcf%C, jcf%epdmax, jcf%m, jcf%Troom, jcf%Tmelt, jcf%T_current])
    call set_hard(2, jcf_law=jcf)
    call vp_sweep('jcfull | defaults | swift', jcf, [real(dp) ::], ep, epd, flow, dflow)

    ! ---- Modified Rusinek-Klepaczko (own hardening).
    do i = 1, 2
        if (i == 1) then   ! fit of test_muscle_vp_models
            mrk = MRK_viscoplastic(B01=100.020701_dp, B02=100.0_dp, epdmax=0.01_dp, &
                nu1=0.01439459_dp, nu2=0.14313556_dp, n0=3.95295384_dp, D2=5.6375e-12_dp, &
                epdmin=1.0e-5_dp, sig_u=253.772550_dp, chi1=9.1503e-4_dp, chi2=0.01942062_dp)
        else               ! dynamic range: epdmax = 1e4, rate-dependent n and saturation
            mrk = MRK_viscoplastic(B01=200.0_dp, B02=150.0_dp, epdmax=1.0e4_dp, &
                nu1=0.05_dp, nu2=0.1_dp, n0=5.0_dp, D2=0.05_dp, &
                epdmin=1.0e-4_dp, sig_u=150.0_dp, chi1=0.05_dp, chi2=0.5_dp)
        end if
        write(key, '(a, i0)') 'mrk | set ', i
        call golden_check(trim(key) // ' | constructor', [mrk%B01, mrk%B02, mrk%epdmax, mrk%nu1, &
            mrk%nu2, mrk%n0, mrk%D2, mrk%epdmin, mrk%chi1, mrk%chi2, mrk%sig_u])
        call vp_sweep(trim(key), mrk, [2.0e-12_dp, mrk%epdmin, 0.5_dp*mrk%epdmax, 0.99_dp*mrk%epdmax, &
            mrk%epdmax, 2.0_dp*mrk%epdmax], ep, epd, flow, dflow)
        n = size(ep)
        call bits_equal(trim(key) // ' flow_stress', [(mrk%flow_stress(ep(i_), epd(i_)), i_ = 1, n)], flow)
        call bits_equal(trim(key) // ' dstress_dep', &
            [(mrk%dstress_dep(ep(i_), epd(i_), DT_VALUES(1)), i_ = 1, n)], dflow)
    end do

    ! ---- Rusinek-Klepaczko (own hardening).
    do i = 1, 2
        if (i == 1) then   ! fit of test_muscle_vp_models
            rk = RK_viscoplastic(B0=18.5805384_dp, epdmax=0.01_dp, nu=0.02014578_dp, &
                ep0=2.36959316_dp, n0=3.05066530_dp, D2=1.5302e-8_dp, D1=0.41425354_dp, &
                epdmin=1.0e-6_dp, sig0=1.1693e-5_dp, m=0.94914871_dp)
        else               ! steel-like: epdmax = 1e7, factor_n reaches 0 at epd = 1e15*epdmin
            rk = RK_viscoplastic(B0=600.0_dp, epdmax=1.0e7_dp, nu=0.05_dp, &
                ep0=0.01_dp, n0=0.3_dp, D2=0.25_dp, D1=0.2_dp, &
                epdmin=1.0e-5_dp, sig0=200.0_dp, m=2.0_dp)
            rk%epmin = 1.0e-8_dp
        end if
        write(key, '(a, i0)') 'rk | set ', i
        call golden_check(trim(key) // ' | constructor', [rk%B0, rk%ep0, rk%n0, rk%D2, rk%D1, rk%nu, &
            rk%epdmin, rk%epdmax, rk%sig0, rk%m, rk%epmin, rk%epmax])
        call vp_sweep(trim(key), rk, [rk%epmin, rk%epdmin, 0.5_dp*rk%epdmax, 0.99_dp*rk%epdmax, &
            rk%epdmax, 2.0_dp*rk%epdmax], ep, epd, flow, dflow)
        n = size(ep)
        call bits_equal(trim(key) // ' flow_stress', [(rk%flow_stress(ep(i_), epd(i_)), i_ = 1, n)], flow)
        call bits_equal(trim(key) // ' dstress_dep', &
            [(rk%dstress_dep(ep(i_), epd(i_), DT_VALUES(1)), i_ = 1, n)], dflow)
    end do

    ! ---- Nemat-Nasser-Li (own hardening).
    do i = 1, 2
        if (i == 1) then   ! fit of test_muscle_vp_models
            nnl = NNL_viscoplastic(sig_a=334.776175_dp, sig_0=5.88318003_dp, KG0=0.21906996_dp, &
                n1=0.06301105_dp, epd0=0.01425422_dp, at=22.8752986_dp, &
                n0=0.97427303_dp, q=2.0_dp, p=0.6666667_dp)
        else               ! large epd0: X clamped to 0.999 at low rates
            nnl = NNL_viscoplastic(sig_a=350.0_dp, sig_0=600.0_dp, KG0=0.04_dp, &
                n1=0.08_dp, epd0=1.0e8_dp, at=10.0_dp, n0=0.5_dp, q=1.5_dp, p=0.5_dp)
        end if
        write(key, '(a, i0)') 'nnl | set ', i
        call golden_check(trim(key) // ' | constructor', [nnl%sig_a, nnl%n1, nnl%sig_0, nnl%KG0, &
            nnl%epd0, nnl%at, nnl%n0, nnl%q, nnl%p, nnl%eps_log, nnl%eps_base])
        call vp_sweep(trim(key), nnl, [nnl%epd0], ep, epd, flow, dflow)
        n = size(ep)
        call bits_equal(trim(key) // ' flow_stress', [(nnl%flow_stress(ep(i_), epd(i_)), i_ = 1, n)], flow)
        call bits_equal(trim(key) // ' dstress_dep', &
            [(nnl%dstress_dep(ep(i_), epd(i_), DT_VALUES(1)), i_ = 1, n)], dflow)
    end do

    ! ---- VA (own hardening).
    do i = 1, 3
        select case (i)
        case (1)   ! fit of test_muscle_vp_models
            va = VA_viscoplastic(B=374.525856_dp, B1=0.00331794_dp, B2=0.50280316_dp, &
                n=0.26412620_dp, m=0.03766668_dp, sig_u=195.233128_dp)
        case (2)   ! m = 0 branch (rate independent)
            va = VA_viscoplastic(B=374.525856_dp, B1=0.00331794_dp, B2=0.50280316_dp, &
                n=0.26412620_dp, m=0.0_dp, sig_u=195.233128_dp)
        case (3)
            va = VA_viscoplastic(B=500.0_dp, B1=0.02_dp, B2=0.1_dp, n=0.3_dp, m=5.0_dp, sig_u=100.0_dp)
        end select
        write(key, '(a, i0)') 'va | set ', i
        call golden_check(trim(key) // ' | constructor', [va%B, va%n, va%B1, va%B2, va%m, va%sig_u])
        call vp_sweep(trim(key), va, [real(dp) ::], ep, epd, flow, dflow)
        n = size(ep)
        call bits_equal(trim(key) // ' flow_stress', [(va%flow_stress(ep(i_), epd(i_)), i_ = 1, n)], flow)
        call bits_equal(trim(key) // ' dstress_dep', &
            [(va%dstress_dep(ep(i_), epd(i_), DT_VALUES(1)), i_ = 1, n)], dflow)
    end do

    ! ---- Perzyna overstress, composed with each hardening law.
    do i = 1, 4
        select case (i)
        case (1)
            per = Perzyna_overstress(eta=1.0e-2_dp, n=1.0_dp)
        case (2)
            per = Perzyna_overstress(eta=1.0_dp, n=0.5_dp)   ! n < 1: infinite slope at the threshold
        case (3)
            per = Perzyna_overstress(eta=1.0e-4_dp, n=5.0_dp)
        case (4)
            per = Perzyna_overstress()                     ! defaults
        end select
        write(key, '(a, i0)') 'perzyna | set ', i
        call golden_check(trim(key) // ' | constructor', [per%eta, per%n])
        do ih = 1, 3
            if (i == 4 .and. ih /= 1) cycle
            call set_hard(ih, per_law=per)
            call os_sweep(trim(key) // ' | ' // trim(HNAME(ih)), per, per%hard_law, ep, sig, flow)
            n = size(ep)
            call bits_equal(trim(key) // ' flow_rate', [(per%flow_rate(sig(i_), ep(i_)), i_ = 1, n)], flow)
        end do
    end do

    ! ---- Chaboche overstress, composed with each hardening law.
    do i = 1, 4
        select case (i)
        case (1)
            cha = Chaboche_overstress(k=80.0_dp, n=2.3_dp)
        case (2)
            cha = Chaboche_overstress(k=50.0_dp, n=1.0_dp)
        case (3)
            cha = Chaboche_overstress(k=200.0_dp, n=0.5_dp)   ! n < 1: infinite slope at the threshold
        case (4)
            cha = Chaboche_overstress()                       ! defaults
        end select
        write(key, '(a, i0)') 'chaboche | set ', i
        call golden_check(trim(key) // ' | constructor', [cha%k, cha%n])
        do ih = 1, 3
            if (i == 4 .and. ih /= 1) cycle
            call set_hard(ih, cha_law=cha)
            call os_sweep(trim(key) // ' | ' // trim(HNAME(ih)), cha, cha%hard_law, ep, sig, flow)
            n = size(ep)
            call bits_equal(trim(key) // ' flow_rate', [(cha%flow_rate(sig(i_), ep(i_)), i_ = 1, n)], flow)
        end do
    end do

    call golden_finish()
    if (n_direct_bad > 0) stop 5

contains

    subroutine set_hard(ih, jc_law, jcf_law, per_law, cha_law)
        !! Assigns hardening law `ih` (1 bilinear, 2 swift, 3 vocemod) to the `hard_law`
        !! component of the given law, as the users do (intrinsic assignment).
        integer, intent(in) :: ih
        type(JC_viscoplastic), intent(inout), optional :: jc_law
        type(JohnsonCook_full_visco_hardening), intent(inout), optional :: jcf_law
        type(Perzyna_overstress), intent(inout), optional :: per_law
        type(Chaboche_overstress), intent(inout), optional :: cha_law
        ! Deallocate first: gfortran 15 does not reallocate a polymorphic allocatable on
        ! intrinsic assignment when only the dynamic type changes (heap overflow).
        if (present(jc_law)) then
            if (allocated(jc_law%hard_law)) deallocate(jc_law%hard_law)
        end if
        if (present(jcf_law)) then
            if (allocated(jcf_law%hard_law)) deallocate(jcf_law%hard_law)
        end if
        if (present(per_law)) then
            if (allocated(per_law%hard_law)) deallocate(per_law%hard_law)
        end if
        if (present(cha_law)) then
            if (allocated(cha_law%hard_law)) deallocate(cha_law%hard_law)
        end if
        select case (ih)
        case (1)
            if (present(jc_law)) jc_law%hard_law = bil
            if (present(jcf_law)) jcf_law%hard_law = bil
            if (present(per_law)) per_law%hard_law = bil
            if (present(cha_law)) cha_law%hard_law = bil
        case (2)
            if (present(jc_law)) jc_law%hard_law = swi
            if (present(jcf_law)) jcf_law%hard_law = swi
            if (present(per_law)) per_law%hard_law = swi
            if (present(cha_law)) cha_law%hard_law = swi
        case (3)
            if (present(jc_law)) jc_law%hard_law = voc
            if (present(jcf_law)) jcf_law%hard_law = voc
            if (present(per_law)) per_law%hard_law = voc
            if (present(cha_law)) cha_law%hard_law = voc
        end select
    end subroutine set_hard

end program test_vp_golden
