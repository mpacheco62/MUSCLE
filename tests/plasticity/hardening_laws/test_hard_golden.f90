! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module hard_golden_cases
    !! Cases of the golden test of the hardening laws: every public type-bound procedure of every
    !! law (`stress`, `dstress_dep`, `ddstress_ddep`) through the polymorphic `Base_hardening_laws`
    !! interface used by the solvers and viscoplastic laws, over an `ep` sweep from 0 to 10.
    use, intrinsic :: iso_fortran_env, only : int64
    use muscle_kinds, only : dp
    use muscle_hard_base, only : Base_hardening_laws
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: EP_SWEEP, sweep, record_bits_equal, n_direct_bad

    real(dp), parameter :: EP_SWEEP(21) = [0.0_dp, 1.0e-12_dp, 1.0e-10_dp, 1.0e-9_dp, 1.0e-6_dp, &
        1.0e-5_dp, 1.0e-4_dp, 1.0e-3_dp, 2.0e-3_dp, 5.0e-3_dp, 1.0e-2_dp, 2.0e-2_dp, 5.0e-2_dp, &
        0.1_dp, 0.2_dp, 0.3_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp]
        !! Common equivalent plastic strain sweep: zero, around the 1e-10 offset of Ludwik's
        !! derivatives, small (1e-6..1e-3), intermediate and large (>= 1) strains.

    integer :: n_direct_bad = 0
        !! Direct (non-polymorphic) calls that differ from the polymorphic ones.

contains

    subroutine sweep(key, h, ep)
        !! Golden records of the input strains and of `stress`, `dstress_dep` and `ddstress_ddep`
        !! at every strain, called through the polymorphic base interface (as the solvers do).
        character(len=*), intent(in) :: key
            !! Key prefix (law and parameter set).
        class(Base_hardening_laws), intent(in) :: h
            !! Hardening law.
        real(dp), intent(in) :: ep(:)
            !! Equivalent plastic strains.
        real(dp) :: s(size(ep)), ds(size(ep)), dds(size(ep))
        integer :: i
        do i = 1, size(ep)
            s(i) = h%stress(ep(i))
            ds(i) = h%dstress_dep(ep(i))
            dds(i) = h%ddstress_ddep(ep(i))
        end do
        call golden_check(key // ' | ep', ep)
        call golden_check(key // ' | stress', s)
        call golden_check(key // ' | dstress_dep', ds)
        call golden_check(key // ' | ddstress_ddep', dds)
    end subroutine sweep

    subroutine record_bits_equal(key, got, expected)
        !! Bitwise check of a direct (concrete type) call against the polymorphic one.
        character(len=*), intent(in) :: key
            !! Description of the result.
        real(dp), intent(in) :: got, expected
            !! Direct and polymorphic results.
        if (transfer(got, 1_int64) /= transfer(expected, 1_int64)) then
            print '(a, 2es25.17)', ' FAIL: direct call differs from polymorphic: ' // key, got, expected
            n_direct_bad = n_direct_bad + 1
        end if
    end subroutine record_bits_equal

end module hard_golden_cases

program test_hard_golden
    !! Golden (characterization) test of the hardening laws in double precision, bit for bit
    !! against tests/reference/plasticity/hardening_laws/hard_laws_dp.txt, captured from the
    !! original OO code (regenerate.sh there). Uses the public OO API as today's solvers do
    !! (keyword structure constructors, `class(Base_hardening_laws), allocatable` with intrinsic
    !! assignment, type-bound calls), so it stays valid with the OO adapters of Phase 6 (D-017.1).
    !!
    !! Exit codes: 1-4 from `golden_finish` (muscle_test_golden); 5: a direct call on the
    !! concrete type differs from the polymorphic call.
    use muscle_kinds, only : dp
    use muscle_hard_base, only : Base_hardening_laws
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_hard_ludwik, only : ludwik_hardening
    use muscle_hard_swift, only : Swift_hardening
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use hard_golden_cases, only : EP_SWEEP, sweep, record_bits_equal, n_direct_bad
    implicit none

    ! Parameter sets (realistic metal values, MPa): one row per set.
    real(dp), parameter :: BIL(2, 4) = reshape([ &
        250.0_dp, 1000.0_dp, &   ! mild steel, linear hardening
        355.0_dp, 0.0_dp, &      ! perfect plasticity
        70.0_dp, 690.0_dp, &     ! aluminium
        480.0_dp, -150.0_dp], [2, 4])  ! softening
        !! Bilinear: y0, K.
    real(dp), parameter :: LUD(3, 4) = reshape([ &
        250.0_dp, 500.0_dp, 0.2_dp, &
        300.0_dp, 800.0_dp, 0.5_dp, &
        200.0_dp, 1000.0_dp, 1.0_dp, &   ! n = 1: linear
        100.0_dp, 300.0_dp, 1.5_dp], [3, 4])  ! n > 1: d2 > 0
        !! Ludwik: sigma0, k, n.
    real(dp), parameter :: SWI(3, 4) = reshape([ &
        500.0_dp, 1.0e-3_dp, 0.2_dp, &
        800.0_dp, 1.0e-2_dp, 0.15_dp, &
        1200.0_dp, 2.0e-3_dp, 0.3_dp, &
        600.0_dp, 0.0_dp, 0.25_dp], [3, 4])   ! e0 = 0 (Hollomon): derivatives are Inf at ep = 0
        !! Swift: k, e0, n.
    real(dp), parameter :: VOC(4, 4) = reshape([ &
        250.0_dp, 100.0_dp, 150.0_dp, 20.0_dp, &
        300.0_dp, 0.0_dp, 200.0_dp, 10.0_dp, &   ! pure Voce (k = 0)
        200.0_dp, 500.0_dp, 100.0_dp, 50.0_dp, &
        180.0_dp, 50.0_dp, 400.0_dp, 2.0_dp], [4, 4])
        !! Modified Voce: Sy, k, q, n.

    class(Base_hardening_laws), allocatable :: h
    type(Bilinear_hardening) :: bil_law
    type(ludwik_hardening) :: lud_law
    type(Swift_hardening) :: swi_law
    type(Voce_modified_hardening) :: voc_law
    real(dp), allocatable :: ep(:)
    character(len=32) :: key
    integer :: i, j

    call golden_init('plasticity/hardening_laws/hard_laws_dp.txt')

    do i = 1, size(BIL, 2)
        write(key, '(a, i0)') 'bilinear | set ', i
        bil_law = Bilinear_hardening(y0=BIL(1, i), K=BIL(2, i))
        call golden_check(trim(key) // ' | constructor', [bil_law%y0, bil_law%K])
        h = bil_law
        call sweep(trim(key), h, EP_SWEEP)
        do j = 1, size(EP_SWEEP)
            call record_bits_equal(trim(key) // ' stress', bil_law%stress(EP_SWEEP(j)), h%stress(EP_SWEEP(j)))
            call record_bits_equal(trim(key) // ' dstress', bil_law%dstress_dep(EP_SWEEP(j)), &
                h%dstress_dep(EP_SWEEP(j)))
            call record_bits_equal(trim(key) // ' ddstress', bil_law%ddstress_ddep(EP_SWEEP(j)), &
                h%ddstress_ddep(EP_SWEEP(j)))
        end do
        deallocate(h)
    end do

    do i = 1, size(LUD, 2)
        write(key, '(a, i0)') 'ludwik | set ', i
        lud_law = ludwik_hardening(sigma0=LUD(1, i), k=LUD(2, i), n=LUD(3, i))
        call golden_check(trim(key) // ' | constructor', [lud_law%sigma0, lud_law%k, lud_law%n])
        h = lud_law
        ! Extra points on both sides of the 1e-10 offset of the derivatives.
        ep = [EP_SWEEP, 5.0e-11_dp, 2.0e-10_dp, 1.0e-8_dp]
        call sweep(trim(key), h, ep)
        do j = 1, size(ep)
            call record_bits_equal(trim(key) // ' stress', lud_law%stress(ep(j)), h%stress(ep(j)))
            call record_bits_equal(trim(key) // ' dstress', lud_law%dstress_dep(ep(j)), h%dstress_dep(ep(j)))
            call record_bits_equal(trim(key) // ' ddstress', lud_law%ddstress_ddep(ep(j)), h%ddstress_ddep(ep(j)))
        end do
        deallocate(h)
    end do

    do i = 1, size(SWI, 2)
        write(key, '(a, i0)') 'swift | set ', i
        swi_law = Swift_hardening(k=SWI(1, i), e0=SWI(2, i), n=SWI(3, i))
        call golden_check(trim(key) // ' | constructor', [swi_law%k, swi_law%e0, swi_law%n])
        h = swi_law
        ! Extra points around the offset e0.
        ep = [EP_SWEEP, 0.5_dp*SWI(2, i), SWI(2, i), 2.0_dp*SWI(2, i)]
        call sweep(trim(key), h, ep)
        do j = 1, size(ep)
            call record_bits_equal(trim(key) // ' stress', swi_law%stress(ep(j)), h%stress(ep(j)))
            call record_bits_equal(trim(key) // ' dstress', swi_law%dstress_dep(ep(j)), h%dstress_dep(ep(j)))
            call record_bits_equal(trim(key) // ' ddstress', swi_law%ddstress_ddep(ep(j)), h%ddstress_ddep(ep(j)))
        end do
        deallocate(h)
    end do

    do i = 1, size(VOC, 2)
        write(key, '(a, i0)') 'vocemod | set ', i
        voc_law = Voce_modified_hardening(Sy=VOC(1, i), k=VOC(2, i), q=VOC(3, i), n=VOC(4, i))
        call golden_check(trim(key) // ' | constructor', [voc_law%Sy, voc_law%k, voc_law%q, voc_law%n])
        h = voc_law
        ! Saturation: ep = 1/n, 5/n, 40/n (exp ~ 4e-18) and 800/n (exp underflows to 0).
        ep = [EP_SWEEP, 1.0_dp/VOC(4, i), 5.0_dp/VOC(4, i), 40.0_dp/VOC(4, i), 800.0_dp/VOC(4, i)]
        call sweep(trim(key), h, ep)
        do j = 1, size(ep)
            call record_bits_equal(trim(key) // ' stress', voc_law%stress(ep(j)), h%stress(ep(j)))
            call record_bits_equal(trim(key) // ' dstress', voc_law%dstress_dep(ep(j)), h%dstress_dep(ep(j)))
            call record_bits_equal(trim(key) // ' ddstress', voc_law%ddstress_ddep(ep(j)), h%ddstress_ddep(ep(j)))
        end do
        deallocate(h)
    end do

    call golden_finish()
    if (n_direct_bad > 0) stop 5
end program test_hard_golden
