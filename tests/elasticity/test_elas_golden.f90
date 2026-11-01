! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module elas_golden_cases
    !! Cases of the golden test of the elasticity family: every public and type-bound procedure of
    !! `Base_elasticity` (deferred ones through `Elasticity_linear`, the only concrete model of the
    !! family) called through the polymorphic `class(Base_elasticity)` interface used by the
    !! return mapping solvers, in 3D and 2D, time-independent and transient (`dtime`, `dstrain`,
    !! `internal`), plus the constructor (`set_parameters` and its cache flags) via `to_string`.
    use, intrinsic :: iso_fortran_env, only : int64
    use muscle_kinds, only : dp, TOL_APPROX_dp
    use muscle_tensors, only : ten_3D2Osym, ten_2D2Osym, ten_3D4O3sym, ten_2D4O3sym, &
        operator(.ddot.), operator(*)
    use muscle_elasticity_base, only : Base_elasticity, Base_internal_elasticity_variables
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: base_defaults, golden_internal, run_3D, run_2D, run_tangents, record_to_string, check_embedding
    public :: n_direct_bad, n_embed_bad, MAGS, DTS, STRAINS_3D, STRAINS_2D

    type, extends(Base_internal_elasticity_variables) :: golden_internal
        !! Concrete internal-variable container (purely elastic models ignore it; recorded after
        !! `stress_and_update` to detect any modification).
        real(dp) :: q(6) = [1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp, 5.0_dp, -6.0_dp]
    end type golden_internal

    type, extends(Base_elasticity) :: base_defaults
        !! Test model that implements only the deferred procedures (delegating to a linear model)
        !! and inherits every default of `Base_elasticity`, so that the defaults overridden by
        !! `Elasticity_linear` (`stress_dev_3D`, `dstress_dstrain_dev_3D`: `.dev.` and projector)
        !! are also covered, as a user model such as `Elasticity_frozen` would inherit them.
        type(Elasticity_linear) :: lin
    contains
        procedure :: stress_3D => defaults_stress_3D
        procedure :: stress_2D => defaults_stress_2D
        procedure :: dstress_dstrain_3D => defaults_dstress_dstrain_3D
        procedure :: dstress_dstrain_2D => defaults_dstress_dstrain_2D
        procedure :: to_string => defaults_to_string
    end type base_defaults

    real(dp), parameter :: MAGS(3) = [1.0e-6_dp, 1.0e-3_dp, 1.0e-1_dp]
        !! Strain magnitudes.
    real(dp), parameter :: DTS(3) = [0.0_dp, 1.0e-3_dp, 10.0_dp]
        !! Time steps of the transient overloads (dstrain = 0 with dt = 0).
    character(len=*), parameter :: STRAINS_3D(5) = [character(len=10) :: 'zero', 'uniaxial', &
        'shear', 'volumetric', 'rotated']
        !! 3D strain states.
    character(len=*), parameter :: STRAINS_2D(6) = [character(len=10) :: 'zero', 'uniaxial', &
        'shear', 'volumetric', 'plane', 'axisym']
        !! 2D strain states (11, 22, 33, 12): plane = general in-plane with e33 = 0, axisym = with
        !! e33 /= 0.

    integer, parameter :: RQ(3, 3) = reshape([-20, 20, 10, 4, -10, 28, 22, 20, 4], [3, 3])
        !! 30 * rotation matrix of the quaternion (1, 2, 3, 4)/sqrt(30) (exact, no libm).

    integer :: n_direct_bad = 0
        !! Specific bindings that differ from the generic call.
    integer :: n_embed_bad = 0
        !! 2D results that differ from the 3D results of the embedded strain.

contains

    pure function defaults_stress_3D(self, strain) result(res)
        class(base_defaults), intent(in) :: self
        class(ten_3D2Osym), intent(in) :: strain
        type(ten_3D2Osym) :: res
        res = self%lin%stress_3D(strain)
    end function defaults_stress_3D

    pure function defaults_stress_2D(self, strain) result(res)
        class(base_defaults), intent(in) :: self
        class(ten_2D2Osym), intent(in) :: strain
        type(ten_2D2Osym) :: res
        res = self%lin%stress_2D(strain)
    end function defaults_stress_2D

    pure function defaults_dstress_dstrain_3D(self, strain) result(res)
        class(base_defaults), intent(in) :: self
        class(ten_3D2Osym), intent(in) :: strain
        type(ten_3D4O3sym) :: res
        res = self%lin%dstress_dstrain_3D(strain)
    end function defaults_dstress_dstrain_3D

    pure function defaults_dstress_dstrain_2D(self, strain) result(res)
        class(base_defaults), intent(in) :: self
        class(ten_2D2Osym), intent(in) :: strain
        type(ten_2D4O3sym) :: res
        res = self%lin%dstress_dstrain_2D(strain)
    end function defaults_dstress_dstrain_2D

    pure function defaults_to_string(self) result(str)
        class(base_defaults), intent(in) :: self
        character(len=275) :: str
        str = self%lin%to_string()
    end function defaults_to_string

    function strain_3D(name, a) result(v)
        !! Voigt (11, 22, 33, 12, 23, 13) components of the 3D strain state `name` with magnitude `a`.
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a
        real(dp) :: v(6), r(3, 3), e(3, 3), dg(3)
        integer :: i, j
        select case (trim(name))
        case ('zero')
            v = 0.0_dp
        case ('uniaxial')
            v = [a, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        case ('shear')
            v = [0.0_dp, 0.0_dp, 0.0_dp, a, 0.0_dp, 0.0_dp]
        case ('volumetric')
            v = [a, a, a, 0.0_dp, 0.0_dp, 0.0_dp]
        case default ! rotated: R diag(a, -0.3a, 0.55a) R^T
            r = real(RQ, dp)/30.0_dp
            dg = [a, -0.3_dp*a, 0.55_dp*a]
            do j = 1, 3
                do i = 1, 3
                    e(i, j) = sum(r(i, :)*dg*r(j, :))
                end do
            end do
            v = [e(1, 1), e(2, 2), e(3, 3), e(1, 2), e(2, 3), e(1, 3)]
        end select
    end function strain_3D

    function strain_2D(name, a) result(v)
        !! Components (11, 22, 33, 12) of the 2D strain state `name` with magnitude `a`.
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: a
        real(dp) :: v(4), c, s, d1, d2
        c = 0.6_dp
        s = 0.8_dp
        d1 = a
        d2 = -0.4_dp*a
        select case (trim(name))
        case ('zero')
            v = 0.0_dp
        case ('uniaxial')
            v = [a, 0.0_dp, 0.0_dp, 0.0_dp]
        case ('shear')
            v = [0.0_dp, 0.0_dp, 0.0_dp, a]
        case ('volumetric')
            v = [a, a, a, 0.0_dp]
        case ('plane') ! in-plane rotation (3-4-5) of diag(a, -0.4a), e33 = 0
            v = [c*c*d1 + s*s*d2, s*s*d1 + c*c*d2, 0.0_dp, c*s*(d1 - d2)]
        case default ! axisym: same with e33 = 0.25a
            v = [c*c*d1 + s*s*d2, s*s*d1 + c*c*d2, 0.25_dp*a, c*s*(d1 - d2)]
        end select
    end function strain_2D

    integer function nmag(name)
        !! Number of magnitudes of a strain state (the zero strain has one).
        character(len=*), intent(in) :: name
        nmag = 3
        if (trim(name) == 'zero') nmag = 1
    end function nmag

    subroutine bits_equal(key, got, expected)
        !! Bitwise check of a specific binding against the generic call.
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: got(:), expected(:)
        if (any(transfer(got, 1_int64, size(got)) /= transfer(expected, 1_int64, size(expected)))) then
            print '(a)', ' FAIL: specific binding differs from the generic call: ' // key
            n_direct_bad = n_direct_bad + 1
        end if
    end subroutine bits_equal

    subroutine run_3D(key, el)
        !! Golden records of every 3D procedure for every 3D strain state and magnitude.
        character(len=*), intent(in) :: key
        class(Base_elasticity), intent(in) :: el
        type(ten_3D2Osym) :: eps, deps, s, sb
        type(ten_3D4O3sym) :: c, cb
        type(golden_internal) :: internal
        real(dp), allocatable :: rin(:), rs(:), rsd(:), rc(:), rcd(:), rdin(:), &
            ts(:), tsd(:), tc(:), tcd(:), us(:), uq(:)
        integer :: is, im, it
        character(len=:), allocatable :: kk

        do is = 1, size(STRAINS_3D)
            kk = key // ' | 3D | ' // trim(STRAINS_3D(is))
            rin = [real(dp) ::]; rs = rin; rsd = rin; rc = rin; rcd = rin
            rdin = rin; ts = rin; tsd = rin; tc = rin; tcd = rin; us = rin; uq = rin
            do im = 1, nmag(STRAINS_3D(is))
                call eps%init(strain_3D(STRAINS_3D(is), MAGS(im)))
                rin = [rin, eps%vals]
                s = el%stress(eps)
                rs = [rs, s%vals]
                sb = el%stress_3D(eps)
                call bits_equal(kk // ' stress_3D', sb%vals, s%vals)
                s = el%stress_dev(eps)
                rsd = [rsd, s%vals]
                sb = el%stress_dev_3D(eps)
                call bits_equal(kk // ' stress_dev_3D', sb%vals, s%vals)
                c = el%dstress_dstrain(eps)
                rc = [rc, c%vals]
                cb = el%dstress_dstrain_3D(eps)
                call bits_equal(kk // ' dstress_dstrain_3D', cb%vals, c%vals)
                c = el%dstress_dstrain_dev(eps)
                rcd = [rcd, c%vals]
                cb = el%dstress_dstrain_dev_3D(eps)
                call bits_equal(kk // ' dstress_dstrain_dev_3D', cb%vals, c%vals)
                do it = 1, size(DTS)
                    ! dt = 0 with a zero increment; dt > 0 with increments 0.1*eps and -eps.
                    select case (it)
                    case (1)
                        call deps%init(spread(0.0_dp, 1, 6))
                    case (2)
                        call deps%init(0.1_dp*eps%vals)
                    case default
                        call deps%init(-eps%vals)
                    end select
                    rdin = [rdin, deps%vals, DTS(it)]
                    s = el%stress(eps, deps, DTS(it), internal)
                    ts = [ts, s%vals]
                    sb = el%stress_3D_transient(eps, deps, DTS(it), internal)
                    call bits_equal(kk // ' stress_3D_transient', sb%vals, s%vals)
                    s = el%stress_dev(eps, deps, DTS(it), internal)
                    tsd = [tsd, s%vals]
                    sb = el%stress_dev_3D_transient(eps, deps, DTS(it), internal)
                    call bits_equal(kk // ' stress_dev_3D_transient', sb%vals, s%vals)
                    c = el%dstress_dstrain(eps, deps, DTS(it), internal)
                    tc = [tc, c%vals]
                    cb = el%dstress_dstrain_3D_transient(eps, deps, DTS(it), internal)
                    call bits_equal(kk // ' dstress_dstrain_3D_transient', cb%vals, c%vals)
                    c = el%dstress_dstrain_dev(eps, deps, DTS(it), internal)
                    tcd = [tcd, c%vals]
                    cb = el%dstress_dstrain_dev_3D_transient(eps, deps, DTS(it), internal)
                    call bits_equal(kk // ' dstress_dstrain_dev_3D_transient', cb%vals, c%vals)
                    call el%stress_and_update(s, eps, deps, DTS(it), internal)
                    us = [us, s%vals]
                    uq = [uq, internal%q]
                    call el%stress_and_update_3D(s, eps, deps, DTS(it), internal)
                    call bits_equal(kk // ' stress_and_update_3D', s%vals, us(size(us)-5:))
                end do
            end do
            call golden_check(kk // ' | strain', rin)
            call golden_check(kk // ' | stress', rs)
            call golden_check(kk // ' | stress_dev', rsd)
            call golden_check(kk // ' | dstress_dstrain', rc)
            call golden_check(kk // ' | dstress_dstrain_dev', rcd)
            call golden_check(kk // ' | transient dstrain, dtime', rdin)
            call golden_check(kk // ' | stress transient', ts)
            call golden_check(kk // ' | stress_dev transient', tsd)
            call golden_check(kk // ' | dstress_dstrain transient', tc)
            call golden_check(kk // ' | dstress_dstrain_dev transient', tcd)
            call golden_check(kk // ' | stress_and_update stress', us)
            call golden_check(kk // ' | stress_and_update internal', uq)
        end do
    end subroutine run_3D

    subroutine run_2D(key, el)
        !! Golden records of every 2D procedure for every 2D strain state and magnitude.
        character(len=*), intent(in) :: key
        class(Base_elasticity), intent(in) :: el
        type(ten_2D2Osym) :: eps, deps, s, sb
        type(ten_2D4O3sym) :: c, cb
        type(golden_internal) :: internal
        real(dp), allocatable :: rin(:), rs(:), rsd(:), rc(:), rdin(:), ts(:), tsd(:), tc(:)
        integer :: is, im, it
        character(len=:), allocatable :: kk

        do is = 1, size(STRAINS_2D)
            kk = key // ' | 2D | ' // trim(STRAINS_2D(is))
            rin = [real(dp) ::]; rs = rin; rsd = rin; rc = rin
            rdin = rin; ts = rin; tsd = rin; tc = rin
            do im = 1, nmag(STRAINS_2D(is))
                call eps%init(strain_2D(STRAINS_2D(is), MAGS(im)))
                rin = [rin, eps%vals]
                s = el%stress(eps)
                rs = [rs, s%vals]
                sb = el%stress_2D(eps)
                call bits_equal(kk // ' stress_2D', sb%vals, s%vals)
                s = el%stress_dev(eps)
                rsd = [rsd, s%vals]
                sb = el%stress_dev_2D(eps)
                call bits_equal(kk // ' stress_dev_2D', sb%vals, s%vals)
                c = el%dstress_dstrain(eps)
                rc = [rc, c%vals]
                cb = el%dstress_dstrain_2D(eps)
                call bits_equal(kk // ' dstress_dstrain_2D', cb%vals, c%vals)
                do it = 1, size(DTS)
                    select case (it)
                    case (1)
                        call deps%init(spread(0.0_dp, 1, 4))
                    case (2)
                        call deps%init(0.1_dp*eps%vals)
                    case default
                        call deps%init(-eps%vals)
                    end select
                    rdin = [rdin, deps%vals, DTS(it)]
                    s = el%stress(eps, deps, DTS(it), internal)
                    ts = [ts, s%vals]
                    sb = el%stress_2D_transient(eps, deps, DTS(it), internal)
                    call bits_equal(kk // ' stress_2D_transient', sb%vals, s%vals)
                    s = el%stress_dev(eps, deps, DTS(it), internal)
                    tsd = [tsd, s%vals]
                    sb = el%stress_dev_2D_transient(eps, deps, DTS(it), internal)
                    call bits_equal(kk // ' stress_dev_2D_transient', sb%vals, s%vals)
                    c = el%dstress_dstrain(eps, deps, DTS(it), internal)
                    tc = [tc, c%vals]
                    cb = el%dstress_dstrain_2D_transient(eps, deps, DTS(it), internal)
                    call bits_equal(kk // ' dstress_dstrain_2D_transient', cb%vals, c%vals)
                end do
            end do
            call golden_check(kk // ' | strain', rin)
            call golden_check(kk // ' | stress', rs)
            call golden_check(kk // ' | stress_dev', rsd)
            call golden_check(kk // ' | dstress_dstrain', rc)
            call golden_check(kk // ' | transient dstrain, dtime', rdin)
            call golden_check(kk // ' | stress transient', ts)
            call golden_check(kk // ' | stress_dev transient', tsd)
            call golden_check(kk // ' | dstress_dstrain transient', tc)
        end do
    end subroutine run_2D

    subroutine run_tangents(key, el)
        !! Golden records of the tangents only (cache flags of `set_parameters`), at the rotated
        !! 3D strain and the axisymmetric 2D strain of magnitude 1e-3, dt = 1e-3.
        character(len=*), intent(in) :: key
        class(Base_elasticity), intent(in) :: el
        type(ten_3D2Osym) :: e3
        type(ten_2D2Osym) :: e2
        type(golden_internal) :: internal
        type(ten_3D4O3sym) :: c3
        type(ten_2D4O3sym) :: c2
        call e3%init(strain_3D('rotated', 1.0e-3_dp))
        call e2%init(strain_2D('axisym', 1.0e-3_dp))
        c3 = el%dstress_dstrain(e3)
        call golden_check(key // ' | 3D | dstress_dstrain', c3%vals)
        c3 = el%dstress_dstrain_dev(e3)
        call golden_check(key // ' | 3D | dstress_dstrain_dev', c3%vals)
        c3 = el%dstress_dstrain(e3, 0.1_dp*e3, 1.0e-3_dp, internal)
        call golden_check(key // ' | 3D | dstress_dstrain transient', c3%vals)
        c3 = el%dstress_dstrain_dev(e3, 0.1_dp*e3, 1.0e-3_dp, internal)
        call golden_check(key // ' | 3D | dstress_dstrain_dev transient', c3%vals)
        c2 = el%dstress_dstrain(e2)
        call golden_check(key // ' | 2D | dstress_dstrain', c2%vals)
        c2 = el%dstress_dstrain(e2, 0.1_dp*e2, 1.0e-3_dp, internal)
        call golden_check(key // ' | 2D | dstress_dstrain transient', c2%vals)
    end subroutine run_tangents

    subroutine record_to_string(key, el)
        !! Golden record of `to_string` (character codes): internal state after construction
        !! (Lame parameters, e1, cache flags), formatted with ES14.6.
        character(len=*), intent(in) :: key
        class(Base_elasticity), intent(in) :: el
        character(len=275) :: str
        integer :: i
        str = el%to_string()
        call golden_check(key // ' | to_string', [(real(iachar(str(i:i)), dp), i = 1, len(str))])
    end subroutine record_to_string

    subroutine close_enough(key, got, expected, scale)
        !! Relative check (TOL_APPROX_dp of the largest expected component, or of `scale`) of an
        !! embedding.
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: got(:), expected(:)
        real(dp), intent(in), optional :: scale
            !! Magnitude of the full result whose part is compared (a deviator that cancels to
            !! rounding noise, e.g. of a volumetric strain, is compared at the scale of the stress).
        real(dp) :: ref
        ref = max(maxval(abs(expected)), tiny(1.0_dp))
        if (present(scale)) ref = max(ref, scale)
        if (any(abs(got - expected) > TOL_APPROX_dp*ref)) then
            print '(a)', ' FAIL: 2D result differs from the 3D result of the embedded strain: ' // key
            print '(a, *(es25.17))', '   2D: ', got
            print '(a, *(es25.17))', '   3D: ', expected
            n_embed_bad = n_embed_bad + 1
        end if
    end subroutine close_enough

    subroutine check_embedding(key, el)
        !! 2D vs 3D (dims-2d-3d.md): every 2D stress/tangent equals the 3D one of the embedded
        !! strain (11, 22, 33, 12, 0, 0); the out-of-plane shears 23, 13 of the 3D result are zero.
        character(len=*), intent(in) :: key
        class(Base_elasticity), intent(in) :: el
        type(ten_2D2Osym) :: e2, d2, s2, p2
        type(ten_3D2Osym) :: e3, d3, s3, p3
        type(golden_internal) :: internal
        real(dp), parameter :: PROBE(4) = [0.7_dp, -1.3_dp, 0.4_dp, 0.9_dp]
            !! 2D tensor to which the tangents are applied (.ddot.).
        real(dp) :: sc
        integer :: is, im
        character(len=:), allocatable :: kk
        call p2%init(PROBE)
        call p3%init([PROBE, 0.0_dp, 0.0_dp])
        do is = 1, size(STRAINS_2D)
            do im = 1, nmag(STRAINS_2D(is))
                kk = key // ' | ' // trim(STRAINS_2D(is))
                call e2%init(strain_2D(STRAINS_2D(is), MAGS(im)))
                call e3%init([e2%vals, 0.0_dp, 0.0_dp])
                call d2%init(0.1_dp*e2%vals)
                call d3%init(0.1_dp*e3%vals)
                s2 = el%stress(e2)
                s3 = el%stress(e3)
                call close_enough(kk // ' stress', [s2%vals, 0.0_dp, 0.0_dp], s3%vals)
                sc = maxval(abs(s3%vals))
                s2 = el%stress_dev(e2)
                s3 = el%stress_dev(e3)
                call close_enough(kk // ' stress_dev', [s2%vals, 0.0_dp, 0.0_dp], s3%vals, sc)
                s2 = el%stress(e2, d2, 1.0e-3_dp, internal)
                s3 = el%stress(e3, d3, 1.0e-3_dp, internal)
                call close_enough(kk // ' stress transient', [s2%vals, 0.0_dp, 0.0_dp], s3%vals)
                s2 = el%stress_dev(e2, d2, 1.0e-3_dp, internal)
                s3 = el%stress_dev(e3, d3, 1.0e-3_dp, internal)
                call close_enough(kk // ' stress_dev transient', [s2%vals, 0.0_dp, 0.0_dp], s3%vals, sc)
                s2 = el%dstress_dstrain(e2) .ddot. p2
                s3 = el%dstress_dstrain(e3) .ddot. p3
                call close_enough(kk // ' dstress_dstrain : probe', [s2%vals, 0.0_dp, 0.0_dp], s3%vals)
                s2 = el%dstress_dstrain(e2, d2, 1.0e-3_dp, internal) .ddot. p2
                s3 = el%dstress_dstrain(e3, d3, 1.0e-3_dp, internal) .ddot. p3
                call close_enough(kk // ' dstress_dstrain transient : probe', [s2%vals, 0.0_dp, 0.0_dp], &
                    s3%vals)
            end do
        end do
    end subroutine check_embedding

end module elas_golden_cases

program test_elas_golden
    !! Golden (characterization) test of the elasticity family in double precision, bit for bit
    !! against tests/reference/elasticity/elasticity_dp.txt, captured from the original OO code
    !! (regenerate.sh there). Uses the public OO API as today's solvers do
    !! (`set_parameters`, `class(Base_elasticity), allocatable` with intrinsic assignment,
    !! generic and specific type-bound calls), so it stays valid with the OO adapters of Phase 6
    !! (D-017.1).
    !!
    !! Exit codes: 1-4 from `golden_finish` (muscle_test_golden); 5: a specific binding differs
    !! from the generic call; 6: a 2D result differs from the 3D result of the embedded strain.
    use muscle_kinds, only : dp
    use muscle_elasticity_base, only : Base_elasticity
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use elas_golden_cases, only : base_defaults, run_3D, run_2D, run_tangents, record_to_string, check_embedding, &
        n_direct_bad, n_embed_bad
    implicit none

    real(dp), parameter :: PRM(2, 6) = reshape([ &
        210000.0_dp, 0.3_dp, &     ! steel
        70000.0_dp, 0.0_dp, &      ! nu = 0: lambda = 0
        1000.0_dp, 0.49_dp, &      ! nearly incompressible
        200.0_dp, 0.4999_dp, &     ! lambda/mu ~ 5e3
        1.0_dp, -0.5_dp, &         ! auxetic, unit modulus
        210000.0_dp, 0.5_dp], [2, 6])  ! incompressible limit: lambda = Inf (anomaly, recorded)
        !! Young's modulus E, Poisson's ratio nu.
    logical, parameter :: FLAGS(2, 4) = reshape([.true., .true., .false., .false., &
        .true., .false., .false., .true.], [2, 4])
        !! Explicit (calc_tan_3D, calc_tan_2D) of `set_parameters`.

    class(Base_elasticity), allocatable :: el
    type(Elasticity_linear) :: lin, unset
    type(base_defaults) :: dflt
    character(len=48) :: key
    integer :: i, j

    call golden_init('elasticity/elasticity_dp.txt')

    ! Default-initialised object (set_parameters never called): every modulus is zero.
    call record_to_string('unset', unset)
    el = unset
    call run_3D('unset', el)
    call run_2D('unset', el)
    deallocate(el)

    do i = 1, size(PRM, 2)
        write(key, '(a, i0)') 'set ', i
        call golden_check(trim(key) // ' | E, nu', PRM(:, i))
        ! Default flags (both stiffness tensors cached), full cases through the polymorphic base.
        call lin%set_parameters(young=PRM(1, i), poisson=PRM(2, i))
        call record_to_string(trim(key), lin)
        el = lin
        call run_3D(trim(key), el)
        call run_2D(trim(key), el)
        if (i < size(PRM, 2)) call check_embedding(trim(key), el)
        deallocate(el)
        ! Same parameters through the inherited defaults of Base_elasticity (3D: the 2D procedures
        ! of Elasticity_linear already are the defaults).
        dflt%lin = lin
        el = dflt
        call run_3D('base defaults | ' // trim(key), el)
        if (i < size(PRM, 2)) call check_embedding('base defaults | ' // trim(key), el)
        deallocate(el)
        ! Explicit cache flags (cached and assembled on the fly): tangents (the only results that
        ! read the cache) and state.
        do j = 1, size(FLAGS, 2)
            write(key, '(a, i0, a, l1, l1)') 'set ', i, ' | calc_tan ', FLAGS(:, j)
            call lin%set_parameters(young=PRM(1, i), poisson=PRM(2, i), calc_tan_3D=FLAGS(1, j), &
                calc_tan_2D=FLAGS(2, j))
            call record_to_string(trim(key), lin)
            el = lin
            call run_tangents(trim(key), el)
            deallocate(el)
        end do
    end do

    ! Re-initialisation of an object: the second call must replace parameters and caches.
    call lin%set_parameters(young=PRM(1, 1), poisson=PRM(2, 1))
    call lin%set_parameters(young=PRM(1, 3), poisson=PRM(2, 3), calc_tan_3D=.false., calc_tan_2D=.false.)
    call record_to_string('reinit cached -> uncached', lin)
    el = lin
    call run_tangents('reinit cached -> uncached', el)
    deallocate(el)
    call lin%set_parameters(young=PRM(1, 1), poisson=PRM(2, 1), calc_tan_3D=.false., calc_tan_2D=.false.)
    call lin%set_parameters(young=PRM(1, 3), poisson=PRM(2, 3))
    call record_to_string('reinit uncached -> cached', lin)
    el = lin
    call run_tangents('reinit uncached -> cached', el)
    deallocate(el)

    call golden_finish()
    if (n_direct_bad > 0) stop 5
    if (n_embed_bad > 0) stop 6
end program test_elas_golden
