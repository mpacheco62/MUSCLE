! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module yield_golden_cases
    !! Cases of the golden test of the yield criteria: every public type-bound procedure of
    !! every criterion through the polymorphic base used by the return-mapping solvers
    !! (`class(Base_yield_critera)`: `stress_eq`, `dstressEq_dstress`, `ddstressEq_ddstress`,
    !! `stressEq_and_derivatives`, `joint_evaluation`) and the finite-difference procedures of
    !! the base (`dstressEq_dstress_numeric`, `ddstressEq_ddstress_numeric`), on a fixed list of
    !! 3D stress states.
    !!
    !! Records per criterion, parameter set and stress state `sNN` (layout of the values):
    !!
    !! - `analytic`: [stress_eq, dstressEq_dstress%vals(1:6), ddstressEq_ddstress%vals(1:21)];
    !! - `base numeric`: [dstressEq_dstress_numeric%vals(1:6), ddstressEq_ddstress_numeric%vals(1:21)];
    !! - `stressEq_and_derivatives` (only when `joint_evaluation()` is true, i.e. the criterion
    !!   overrides it): [f, df%vals(1:6), ddf%vals(1:21)]. When it is false the inherited default
    !!   must give the same bits as the three bindings (checked bit for bit, exit code 5).
    !!
    !! A NaN is recorded as one canonical quiet NaN: the golden fixes where the original code
    !! gives NaN (e.g. the derivatives at a zero equivalent stress), not its sign or payload,
    !! which depend on the order of the operations that produce it.
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym, ten_3D4O3sym
    use muscle_yield_base, only : Base_yield_critera
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: yield_sweep, record_stresses, n_joint_bad, NSTRESS

    integer, parameter :: NBASE = 22
        !! Stress states of unit magnitude.
    integer, parameter :: NSCALED = 5
        !! States repeated at other magnitudes.
    real(dp), parameter :: MAGNITUDES(3) = [1.0e-3_dp, 250.0_dp, 1.0e5_dp]
        !! Magnitudes of the scaled states (the unit ones are magnitude 1).
    integer, parameter :: NSTRESS = NBASE + NSCALED*size(MAGNITUDES)
        !! Number of stress states.
    real(dp), parameter :: T3 = 1.0_dp/3.0_dp
        !! One third (components of the rotated states).
    real(dp), parameter :: BASE(6, NBASE) = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &           ! 01 zero tensor
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &           ! 02 uniaxial tension x
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &           ! 03 uniaxial tension y
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &           ! 04 uniaxial tension z
        -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &          ! 05 uniaxial compression x
        0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, &           ! 06 uniaxial at 45 deg in xy
        0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp, &           ! 07 uniaxial at 45 deg in yz
        1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &           ! 08 equibiaxial (repeated)
        1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &           ! 09 biaxial 2:1
        1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &          ! 10 tension-compression
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &           ! 11 pure shear xy
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &           ! 12 pure shear yz
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &           ! 13 pure shear xz
        1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &           ! 14 hydrostatic tension
        -1.0_dp, -1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &        ! 15 hydrostatic compression
        1.0_dp, 1.0_dp, 1.0_dp, 1.0e-10_dp, 0.0_dp, 0.0_dp, &       ! 16 hydrostatic + tiny shear
        3.0_dp, -1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &          ! 17 triaxial, distinct
        2.0_dp, -1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &         ! 18 axisymmetric (repeated)
        1.205_dp, -0.3525_dp, 0.6075_dp, 0.425_dp, -0.28_dp, 0.3325_dp, & ! 19 general
        T3, T3, T3, T3, T3, T3, &                                   ! 20 uniaxial along (1,1,1)
        2.0_dp*T3, 2.0_dp*T3, 2.0_dp*T3, -T3, -T3, -T3, &           ! 21 equibiaxial normal (1,1,1)
        0.2_dp, -0.6_dp, 0.1_dp, 0.9_dp, -0.4_dp, 0.7_dp], [6, NBASE])   ! 22 general, shear dominated
        !! Unit stress states (Voigt tensorial xx, yy, zz, xy, yz, xz). 08, 14-16, 18, 20 and 21
        !! have repeated principal values (CPB06/CPB06ex2/Yld2004 decompose spectrally).
    integer, parameter :: SCALED(NSCALED) = [2, 8, 11, 14, 19]
        !! Unit states repeated at each of MAGNITUDES.

    integer :: n_joint_bad = 0
        !! Sets whose inherited stressEq_and_derivatives differs from the three bindings.

contains

    subroutine check(key, x)
        !! `golden_check` of `x` with every NaN replaced by the canonical quiet NaN.
        character(len=*), intent(in) :: key
            !! Key of the record.
        real(dp), intent(in) :: x(:)
            !! Values.
        call golden_check(key, canon(x))
    end subroutine check

    function canon(x) result(res)
        !! `x` with every NaN replaced by the canonical quiet NaN.
        real(dp), intent(in) :: x(:)
            !! Values.
        real(dp) :: res(size(x))
        res = merge(ieee_value(1.0_dp, ieee_quiet_nan), x, ieee_is_nan(x))
    end function canon

    pure function stress_case(j) result(s)
        !! j-th stress state: the unit states, then each scaled state at each magnitude.
        integer, intent(in) :: j
        type(ten_3D2Osym) :: s
        integer :: m, i
        if (j <= NBASE) then
            s%vals = BASE(:, j)
        else
            m = (j - NBASE - 1)/NSCALED + 1
            i = mod(j - NBASE - 1, NSCALED) + 1
            s%vals = MAGNITUDES(m)*BASE(:, SCALED(i))
        end if
    end function stress_case

    subroutine record_stresses()
        !! Golden records of the input stress states.
        integer :: j
        type(ten_3D2Osym) :: s
        character(len=8) :: tag
        do j = 1, NSTRESS
            s = stress_case(j)
            write(tag, '(a, i2.2)') 's', j
            call golden_check('stress | ' // trim(tag), s%vals)
        end do
    end subroutine record_stresses

    subroutine yield_sweep(key, y, states)
        !! Golden records of every procedure of `y` on the stress states (all when `states` is
        !! absent), through the polymorphic base.
        character(len=*), intent(in) :: key
            !! Key prefix (criterion, parameter set).
        class(Base_yield_critera), intent(in) :: y
            !! Yield criterion.
        integer, intent(in), optional :: states(:)
            !! Subset of stress states.
        type(ten_3D2Osym) :: s, df, ndf, jdf
        type(ten_3D4O3sym) :: ddf, nddf, jddf
        real(dp) :: f, jf
        logical :: joint
        character(len=8) :: tag
        integer :: j, ij, nj
        integer, allocatable :: list(:)

        if (present(states)) then
            list = states
        else
            list = [(j, j = 1, NSTRESS)]
        end if
        nj = size(list)
        joint = y%joint_evaluation()
        call check(key // ' | joint_evaluation', [merge(1.0_dp, 0.0_dp, joint)])
        do ij = 1, nj
            j = list(ij)
            s = stress_case(j)
            write(tag, '(a, i2.2)') 's', j
            f = y%stress_eq(s)
            df = y%dstressEq_dstress(s)
            ddf = y%ddstressEq_ddstress(s)
            call check(key // ' | ' // trim(tag) // ' | analytic', [f, df%vals, ddf%vals])
            ndf = y%dstressEq_dstress_numeric(s)
            nddf = y%ddstressEq_ddstress_numeric(s)
            call check(key // ' | ' // trim(tag) // ' | base numeric', [ndf%vals, nddf%vals])
            call y%stressEq_and_derivatives(s, jf, jdf, jddf)
            if (joint) then
                call check(key // ' | ' // trim(tag) // ' | stressEq_and_derivatives', &
                    [jf, jdf%vals, jddf%vals])
            else if (any(transfer(canon([jf, jdf%vals, jddf%vals]), 1_int64, 28) /= &
                         transfer(canon([f, df%vals, ddf%vals]), 1_int64, 28))) then
                print '(a)', ' FAIL: inherited stressEq_and_derivatives differs from the bindings: ' &
                    // key // ' | ' // trim(tag)
                n_joint_bad = n_joint_bad + 1
            end if
        end do
    end subroutine yield_sweep

end module yield_golden_cases

program test_yield_golden
    !! Golden (characterization) test of the yield criteria in double precision, bit for bit
    !! against tests/reference/plasticity/yield_criteria/yield_criteria_dp.txt, captured from
    !! the original OO code (regenerate.sh there). Uses the public OO API as today's solvers do
    !! (structure constructor of Hill48, `init` of the others, polymorphic
    !! `class(Base_yield_critera)` calls), so it stays valid with the OO adapters of Phase 6
    !! (D-017.1). Constructor records hold the public components after construction
    !! (DruckerPrager's are private: only its behaviour is recorded, including the default
    !! components of an uninitialised object).
    !!
    !! Exit codes: 1-4 from `golden_finish` (muscle_test_golden); 5: the inherited
    !! `stressEq_and_derivatives` differs from the three bindings.
    use muscle_kinds, only : dp
    use muscle_yield_vonmises, only : VonMises
    use muscle_yield_hill48, only : Hill48
    use muscle_yield_druckerprager, only : DruckerPrager, DP_HARDENING_TENSION, &
        DP_HARDENING_COMPRESSION, DP_HARDENING_SHEAR
    use muscle_yield_cpb06, only : CPB06
    use muscle_yield_cpb06ex2, only : CPB06ex2
    use muscle_yield_yld2004, only : Yld2004
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use yield_golden_cases, only : yield_sweep, record_stresses, n_joint_bad
    implicit none

    real(dp), parameter :: HILLP(6, 3) = reshape([ &
        0.5_dp, 0.5_dp, 0.5_dp, 1.5_dp, 1.5_dp, 1.5_dp, &     ! isotropic (von Mises)
        0.25_dp, 0.4_dp, 0.6_dp, 1.3_dp, 1.7_dp, 1.45_dp, &   ! anisotropic
        0.7_dp, 0.2_dp, 0.35_dp, 2.0_dp, 0.9_dp, 1.1_dp], [6, 3])  ! strongly anisotropic
        !! Hill48: f, g, h, l, m, n.
    real(dp), parameter :: DPP(2, 5) = reshape([ &
        20.0_dp, 0.8_dp, &      ! tension calibration
        20.0_dp, 0.8_dp, &      ! compression calibration
        30.0_dp, 1.0_dp, &      ! shear calibration, K = 1
        0.0_dp, 1.0_dp, &       ! beta = 0, K = 1: von Mises
        45.0_dp, 0.778_dp], [2, 5])  ! large friction angle, lowest K
        !! Drucker-Prager: beta (deg), K.
    integer, parameter :: DPMODE(5) = [DP_HARDENING_TENSION, DP_HARDENING_COMPRESSION, &
        DP_HARDENING_SHEAR, DP_HARDENING_TENSION, DP_HARDENING_TENSION]
        !! Drucker-Prager calibration mode of each set.
    real(dp), parameter :: CPBP(14, 6) = reshape([ &
        ! 1: isotropic limit (C = I, k = 0, a = 2), equal to von Mises
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, &
        ! 2: Ti-6Al-4V (Tuninetti et al., 2013, Table 2), tension-compression asymmetry
        1.0_dp, -2.373_dp, -2.364_dp, -2.373_dp, -1.838_dp, 1.196_dp, -2.364_dp, 1.196_dp, -2.444_dp, &
        3.607_dp, 3.607_dp, 3.607_dp, -0.136_dp, 2.0_dp, &
        ! 3: non-symmetric C1, three different shear factors, positive k
        1.0_dp, -2.373_dp, -2.364_dp, -1.5_dp, -1.838_dp, 1.196_dp, -2.364_dp, 0.8_dp, -2.444_dp, &
        3.607_dp, 2.9_dp, 4.1_dp, 0.25_dp, 2.0_dp, &
        ! 4: material 2 with a = 8
        1.0_dp, -2.373_dp, -2.364_dp, -2.373_dp, -1.838_dp, 1.196_dp, -2.364_dp, 1.196_dp, -2.444_dp, &
        3.607_dp, 3.607_dp, 3.607_dp, -0.136_dp, 8.0_dp, &
        ! 5: isotropic, k = 1
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
        ! 6: isotropic, k = -1
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 2.0_dp], [14, 6])
        !! CPB06: c11, c12, c13, c21, c22, c23, c31, c32, c33, c44, c55, c66, k, a.
    real(dp), parameter :: EX2P(27, 5) = reshape([ &
        ! 1: steel (Plunkett et al., 2008, Table 2)
        1.283_dp, 0.208_dp, -0.421_dp, 0.208_dp, 0.535_dp, -0.791_dp, -0.421_dp, -0.791_dp, 0.746_dp, &
        1.156_dp, 1.418_dp, 1.329_dp, 0.030_dp, &
        1.283_dp, 0.380_dp, -0.001_dp, 0.380_dp, 2.145_dp, 0.792_dp, -0.001_dp, 0.792_dp, 0.218_dp, &
        1.721_dp, 1.418_dp, 1.329_dp, -0.027_dp, 5.0_dp, &
        ! 2: 6111-T4 (Table 8), a = 12
        0.250_dp, -1.399_dp, -0.971_dp, -1.399_dp, -0.375_dp, -2.254_dp, -0.971_dp, -2.254_dp, -0.831_dp, &
        1.316_dp, 1.403_dp, 1.101_dp, 0.0_dp, &
        0.250_dp, -0.840_dp, -1.330_dp, -0.840_dp, 0.780_dp, -0.565_dp, -1.330_dp, -0.565_dp, 0.479_dp, &
        1.413_dp, 1.094_dp, 1.448_dp, 0.0_dp, 12.0_dp, &
        ! 3: AZ31B (Table 1), strong asymmetry
        0.771_dp, -0.914_dp, -1.312_dp, -0.914_dp, 0.654_dp, -1.209_dp, -1.312_dp, -1.209_dp, -0.293_dp, &
        1.617_dp, 1.0_dp, 1.0_dp, 0.295_dp, &
        0.771_dp, 0.337_dp, -0.764_dp, 0.337_dp, 0.646_dp, -0.401_dp, -0.764_dp, -0.401_dp, -0.022_dp, &
        0.499_dp, 1.0_dp, 1.0_dp, 0.900_dp, 4.0_dp, &
        ! 4: two equal non-symmetric transformations (reduces to CPB06), a = 4
        1.0_dp, -0.2_dp, 0.3_dp, 0.15_dp, 0.8_dp, -0.4_dp, -0.1_dp, 0.2_dp, 1.2_dp, &
        1.3_dp, 0.7_dp, 1.1_dp, 0.25_dp, &
        1.0_dp, -0.2_dp, 0.3_dp, 0.15_dp, 0.8_dp, -0.4_dp, -0.1_dp, 0.2_dp, 1.2_dp, &
        1.3_dp, 0.7_dp, 1.1_dp, 0.25_dp, 4.0_dp, &
        ! 5: isotropic, k = kp = 0, a = 2 (von Mises)
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp], [27, 5])
        !! CPB06ex2: c11..c33, c44, c55, c66, k, cp11..cp33, cp44, cp55, cp66, kp, a.
    real(dp), parameter :: YLDP(19, 4) = reshape([ &
        ! 1: isotropic (all coefficients 1), a = 2
        1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
        ! 2: AA2090-T3 (Barlat et al., 2005), a = 8
        -0.069888_dp, 0.936408_dp, 0.079143_dp, 1.003060_dp, 0.524741_dp, 1.363180_dp, &
        0.954322_dp, 1.023770_dp, 1.069060_dp, &
        0.981171_dp, 0.476741_dp, 0.575316_dp, 0.866827_dp, 1.145010_dp, -0.079294_dp, &
        1.404620_dp, 1.051660_dp, 1.147100_dp, 8.0_dp, &
        ! 3: C' isotropic, C'' of AA2090-T3 with c''12 = c''13 = 1, a = 2 (zero differences)
        1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 0.575316_dp, 0.866827_dp, 1.145010_dp, -0.079294_dp, &
        1.404620_dp, 1.051660_dp, 1.147100_dp, 2.0_dp, &
        ! 4: AA2090-T3 with a = 6
        -0.069888_dp, 0.936408_dp, 0.079143_dp, 1.003060_dp, 0.524741_dp, 1.363180_dp, &
        0.954322_dp, 1.023770_dp, 1.069060_dp, &
        0.981171_dp, 0.476741_dp, 0.575316_dp, 0.866827_dp, 1.145010_dp, -0.079294_dp, &
        1.404620_dp, 1.051660_dp, 1.147100_dp, 6.0_dp], [19, 4])
        !! Yld2004: cp12, cp13, cp21, cp23, cp31, cp32, cp44, cp55, cp66, cpp12..cpp66, a.

    type(VonMises) :: vm
    type(Hill48) :: h48
    type(DruckerPrager) :: dpr, dp0
    type(CPB06) :: cpb
    type(CPB06ex2) :: ex2
    type(Yld2004) :: yld
    character(len=48) :: key
    integer :: i
    real(dp) :: c(27)

    call golden_init('plasticity/yield_criteria/yield_criteria_dp.txt')
    call record_stresses()

    ! ---- von Mises (no parameters).
    call yield_sweep('vonmises', vm)

    ! ---- Hill48 (keyword structure constructor).
    do i = 1, size(HILLP, 2)
        h48 = Hill48(f=HILLP(1, i), g=HILLP(2, i), h=HILLP(3, i), l=HILLP(4, i), m=HILLP(5, i), &
                     n=HILLP(6, i))
        write(key, '(a, i0)') 'hill48 | set ', i
        call golden_check(trim(key) // ' | constructor', [h48%f, h48%g, h48%h, h48%l, h48%m, h48%n])
        call yield_sweep(trim(key), h48)
    end do

    ! ---- Drucker-Prager (init; private components).
    do i = 1, size(DPP, 2)
        call dpr%init(beta_deg=DPP(1, i), K=DPP(2, i), hardening_mode=DPMODE(i))
        write(key, '(a, i0)') 'druckerprager | set ', i
        call yield_sweep(trim(key), dpr)
    end do
    ! Default components (no init): f_factor = I1_factor = K = 0; only the zero and hydrostatic
    ! states, where the result does not involve 1/K.
    call yield_sweep('druckerprager | defaults', dp0, [1, 14, 15, 16])

    ! ---- CPB06 (init).
    do i = 1, size(CPBP, 2)
        c(1:14) = CPBP(:, i)
        call cpb%init(c11=c(1), c12=c(2), c13=c(3), c21=c(4), c22=c(5), c23=c(6), &
                      c31=c(7), c32=c(8), c33=c(9), c44=c(10), c55=c(11), c66=c(12), k=c(13), a=c(14))
        write(key, '(a, i0)') 'cpb06 | set ', i
        call golden_check(trim(key) // ' | constructor', &
            [reshape(cpb%C1, [9]), cpb%C2, cpb%k, cpb%a, cpb%B, reshape(cpb%L, [9])])
        call yield_sweep(trim(key), cpb)
    end do

    ! ---- CPB06ex2 (init, positional as its users).
    do i = 1, size(EX2P, 2)
        c = EX2P(:, i)
        call ex2%init(c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10), c(11), c(12), &
                      c(13), c(14), c(15), c(16), c(17), c(18), c(19), c(20), c(21), c(22), c(23), &
                      c(24), c(25), c(26), c(27))
        write(key, '(a, i0)') 'cpb06ex2 | set ', i
        call golden_check(trim(key) // ' | constructor', [reshape(ex2%L, [9]), reshape(ex2%Lp, [9]), &
            ex2%c_shear, ex2%cp_shear, ex2%k, ex2%kp, ex2%a, ex2%B])
        call yield_sweep(trim(key), ex2)
    end do

    ! ---- Yld2004-18p (init; overrides stressEq_and_derivatives).
    do i = 1, size(YLDP, 2)
        c(1:19) = YLDP(:, i)
        call yld%init(cp12=c(1), cp13=c(2), cp21=c(3), cp23=c(4), cp31=c(5), cp32=c(6), &
                      cp44=c(7), cp55=c(8), cp66=c(9), cpp12=c(10), cpp13=c(11), cpp21=c(12), &
                      cpp23=c(13), cpp31=c(14), cpp32=c(15), cpp44=c(16), cpp55=c(17), cpp66=c(18), &
                      a=c(19))
        write(key, '(a, i0)') 'yld2004 | set ', i
        call golden_check(trim(key) // ' | constructor', [reshape(yld%Lp, [9]), yld%cp_shear, &
            reshape(yld%Lpp, [9]), yld%cpp_shear, yld%a])
        call yield_sweep(trim(key), yld)
    end do

    if (n_joint_bad > 0) then
        call golden_finish()
        stop 5
    end if
    call golden_finish()
end program test_yield_golden
