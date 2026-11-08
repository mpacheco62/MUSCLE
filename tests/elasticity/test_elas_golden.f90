! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module elas_golden_cases
    !! Cases of the golden test of the elasticity family: every function of the family (`elas_stress`,
    !! `elas_stress_dev`, `elas_dstress_dstrain`, `elas_dstress_dstrain_dev`) of the linear model
    !! through the kernel API, in 3D and 2D.
    !!
    !! Removed with the OO layer (D-018.8): the records of `to_string`, of the cache flags of
    !! `set_parameters` and of the defaults of Base_elasticity (its transient overloads, also those
    !! that Elasticity_linear overrode, `stress_and_update`, and the test model that inherited every
    !! default).
    use muscle_kinds, only : dp, TOL_APPROX_dp
    use muscle_tensors, only : ten_3D2Osym, ten_2D2Osym, ten_3D4O3sym, ten_2D4O3sym, operator(.ddot.)
    use muscle_elasticity_laws, only : elas_stress, elas_stress_dev, elas_dstress_dstrain, &
        elas_dstress_dstrain_dev, elas_linear_dp
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: run_3D, run_2D, check_embedding
    public :: n_embed_bad, MAGS, STRAINS_3D, STRAINS_2D

    real(dp), parameter :: MAGS(3) = [1.0e-6_dp, 1.0e-3_dp, 1.0e-1_dp]
        !! Strain magnitudes.
    character(len=*), parameter :: STRAINS_3D(5) = [character(len=10) :: 'zero', 'uniaxial', &
        'shear', 'volumetric', 'rotated']
        !! 3D strain states.
    character(len=*), parameter :: STRAINS_2D(6) = [character(len=10) :: 'zero', 'uniaxial', &
        'shear', 'volumetric', 'plane', 'axisym']
        !! 2D strain states (11, 22, 33, 12): plane = general in-plane with e33 = 0, axisym = with
        !! e33 /= 0.

    integer, parameter :: RQ(3, 3) = reshape([-20, 20, 10, 4, -10, 28, 22, 20, 4], [3, 3])
        !! 30 * rotation matrix of the quaternion (1, 2, 3, 4)/sqrt(30) (exact, no libm).

    integer :: n_embed_bad = 0
        !! 2D results that differ from the 3D results of the embedded strain.

contains

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

    subroutine run_3D(key, el)
        !! Golden records of every 3D function for every 3D strain state and magnitude.
        character(len=*), intent(in) :: key
        type(elas_linear_dp), intent(in) :: el
        type(ten_3D2Osym) :: eps, s
        type(ten_3D4O3sym) :: c
        real(dp), allocatable :: rin(:), rs(:), rsd(:), rc(:), rcd(:)
        integer :: is, im
        character(len=:), allocatable :: kk

        do is = 1, size(STRAINS_3D)
            kk = key // ' | 3D | ' // trim(STRAINS_3D(is))
            rin = [real(dp) ::]; rs = rin; rsd = rin; rc = rin; rcd = rin
            do im = 1, nmag(STRAINS_3D(is))
                call eps%init(strain_3D(STRAINS_3D(is), MAGS(im)))
                rin = [rin, eps%vals]
                s = elas_stress(el, eps)
                rs = [rs, s%vals]
                s = elas_stress_dev(el, eps)
                rsd = [rsd, s%vals]
                c = elas_dstress_dstrain(el, eps)
                rc = [rc, c%vals]
                c = elas_dstress_dstrain_dev(el, eps)
                rcd = [rcd, c%vals]
            end do
            call golden_check(kk // ' | strain', rin)
            call golden_check(kk // ' | stress', rs)
            call golden_check(kk // ' | stress_dev', rsd)
            call golden_check(kk // ' | dstress_dstrain', rc)
            call golden_check(kk // ' | dstress_dstrain_dev', rcd)
        end do
    end subroutine run_3D

    subroutine run_2D(key, el)
        !! Golden records of every 2D function for every 2D strain state and magnitude.
        character(len=*), intent(in) :: key
        type(elas_linear_dp), intent(in) :: el
        type(ten_2D2Osym) :: eps, s
        type(ten_2D4O3sym) :: c
        real(dp), allocatable :: rin(:), rs(:), rsd(:), rc(:)
        integer :: is, im
        character(len=:), allocatable :: kk

        do is = 1, size(STRAINS_2D)
            kk = key // ' | 2D | ' // trim(STRAINS_2D(is))
            rin = [real(dp) ::]; rs = rin; rsd = rin; rc = rin
            do im = 1, nmag(STRAINS_2D(is))
                call eps%init(strain_2D(STRAINS_2D(is), MAGS(im)))
                rin = [rin, eps%vals]
                s = elas_stress(el, eps)
                rs = [rs, s%vals]
                s = elas_stress_dev(el, eps)
                rsd = [rsd, s%vals]
                c = elas_dstress_dstrain(el, eps)
                rc = [rc, c%vals]
            end do
            call golden_check(kk // ' | strain', rin)
            call golden_check(kk // ' | stress', rs)
            call golden_check(kk // ' | stress_dev', rsd)
            call golden_check(kk // ' | dstress_dstrain', rc)
        end do
    end subroutine run_2D

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
        type(elas_linear_dp), intent(in) :: el
        type(ten_2D2Osym) :: e2, s2, p2
        type(ten_3D2Osym) :: e3, s3, p3
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
                s2 = elas_stress(el, e2)
                s3 = elas_stress(el, e3)
                call close_enough(kk // ' stress', [s2%vals, 0.0_dp, 0.0_dp], s3%vals)
                sc = maxval(abs(s3%vals))
                s2 = elas_stress_dev(el, e2)
                s3 = elas_stress_dev(el, e3)
                call close_enough(kk // ' stress_dev', [s2%vals, 0.0_dp, 0.0_dp], s3%vals, sc)
                s2 = elas_dstress_dstrain(el, e2) .ddot. p2
                s3 = elas_dstress_dstrain(el, e3) .ddot. p3
                call close_enough(kk // ' dstress_dstrain : probe', [s2%vals, 0.0_dp, 0.0_dp], s3%vals)
            end do
        end do
    end subroutine check_embedding

end module elas_golden_cases

program test_elas_golden
    !! Golden (characterization) test of the elasticity family in double precision, within the tolerance of its category (D-021)
    !! against tests/reference/elasticity/elasticity_dp.txt, captured from the original OO code
    !! (regenerate.sh there) and reproduced by the kernel API (D-018.18f): `elas_linear_init` and the
    !! family generics of muscle_elasticity_laws.
    !!
    !! Exit codes: 1-4 from `golden_finish` (muscle_test_golden); 6: a 2D result differs from the 3D
    !! result of the embedded strain.
    use muscle_kinds, only : dp
    use muscle_elasticity_laws, only : elas_linear_dp
    use muscle_elasticity_linear, only : elas_linear_init
    use muscle_status, only : STATUS_OK
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use elas_golden_cases, only : run_3D, run_2D, check_embedding, n_embed_bad
    implicit none

    real(dp), parameter :: PRM(2, 6) = reshape([ &
        210000.0_dp, 0.3_dp, &     ! steel
        70000.0_dp, 0.0_dp, &      ! nu = 0: lambda = 0
        1000.0_dp, 0.49_dp, &      ! nearly incompressible
        200.0_dp, 0.4999_dp, &     ! lambda/mu ~ 5e3
        1.0_dp, -0.5_dp, &         ! auxetic, unit modulus
        210000.0_dp, 0.5_dp], [2, 6])  ! incompressible limit: lambda = Inf (anomaly, recorded)
        !! Young's modulus E, Poisson's ratio nu.

    type(elas_linear_dp) :: lin, unset
    character(len=48) :: key
    integer :: i, status

    call golden_init('elasticity/elasticity_dp.txt')

    ! Default-initialised parameters (elas_linear_init never called): every modulus is zero.
    call run_3D('unset', unset)
    call run_2D('unset', unset)

    do i = 1, size(PRM, 2)
        write(key, '(a, i0)') 'set ', i
        call golden_check(trim(key) // ' | E, nu', PRM(:, i))
        call elas_linear_init(lin, young=PRM(1, i), poisson=PRM(2, i), status=status)
        if (status /= STATUS_OK) error stop 'elas_linear_init failed'
        call run_3D(trim(key), lin)
        call run_2D(trim(key), lin)
        if (i < size(PRM, 2)) call check_embedding(trim(key), lin)
    end do

    call golden_finish()
    if (n_embed_bad > 0) stop 6
end program test_elas_golden
