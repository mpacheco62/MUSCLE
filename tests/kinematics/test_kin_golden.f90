! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module kin_golden_cases
    !! Cases of the golden test of the kinematics family: every public and type-bound procedure of
    !! `Base_kinematics` (deferred bindings and `rotate_to/from_corotational`, called through the
    !! polymorphic `class(Base_kinematics)`), of `Base_F_kinematics` (update by F, incremental update,
    !! `update_F_base`, the cached state, the strain getters and the ten stress conversions, through
    !! `class(Base_F_kinematics)`) and of the six concrete kinematics (Small_strain, Total_lagrangian,
    !! Updated_lagrangian, Material_logarithmic, Spatial_logarithmic, Corotational), on deterministic
    !! deformation gradients (identity, rotations up to 180 degrees, stretches with repeated and
    !! distinct eigenvalues of C, simple shear, strong compression, nearly singular F, det F < 0) and
    !! incremental sequences (F_n, F_{n+1}) with small and large rotations. The kinematics family is
    !! 3D only (no 2D type exists).
    use, intrinsic :: iso_fortran_env, only : int64
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2O, ten_3D2Osym, ten_3D4O2sym
    use muscle_kinematics, only : Base_kinematics, Base_F_kinematics, Small_strain_kinematics, &
        Total_lagrangian_kinematics, Updated_lagrangian_kinematics, Corotational_kinematics, &
        Material_logarithmic_kinematics, Spatial_logarithmic_kinematics
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: NCASES, case_name, case_F, rotation, record_base, record_finite, bits_equal, &
        n_alt_bad, input_S, input_S2, input_P, input_C, input_R, kin_name, new_finite, NFINITE

    integer, parameter :: NCASES = 20
        !! Number of deformation gradients.
    integer, parameter :: NFINITE = 4
        !! Number of concrete F-driven kinematics.
    real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp

    integer :: n_alt_bad = 0
        !! Alternative paths (optional arguments, generic vs specific binding) that differ.

contains

    function kin_name(ik) result(name)
        !! Name of the F-driven kinematics `ik`.
        integer, intent(in) :: ik
        character(len=:), allocatable :: name
        select case (ik)
        case (1); name = 'total_lagrangian'
        case (2); name = 'updated_lagrangian'
        case (3); name = 'material_logarithmic'
        case default; name = 'spatial_logarithmic'
        end select
    end function kin_name

    subroutine new_finite(ik, kin)
        !! Allocates a fresh F-driven kinematics `ik`.
        integer, intent(in) :: ik
        class(Base_F_kinematics), allocatable, intent(out) :: kin
        select case (ik)
        case (1); allocate(Total_lagrangian_kinematics :: kin)
        case (2); allocate(Updated_lagrangian_kinematics :: kin)
        case (3); allocate(Material_logarithmic_kinematics :: kin)
        case default; allocate(Spatial_logarithmic_kinematics :: kin)
        end select
    end subroutine new_finite

    function rotation(axis, deg) result(Q)
        !! Rotation of `deg` degrees about `axis` (Rodrigues), column-major 3x3 values.
        real(dp), intent(in) :: axis(3), deg
        real(dp) :: Q(9)
        real(dp) :: n(3), c, s, t, m(3, 3)
        n = axis/sqrt(sum(axis**2))
        c = cos(deg*PI/180.0_dp)
        s = sin(deg*PI/180.0_dp)
        t = 1.0_dp - c
        m(1, :) = [t*n(1)*n(1) + c, t*n(1)*n(2) - s*n(3), t*n(1)*n(3) + s*n(2)]
        m(2, :) = [t*n(1)*n(2) + s*n(3), t*n(2)*n(2) + c, t*n(2)*n(3) - s*n(1)]
        m(3, :) = [t*n(1)*n(3) - s*n(2), t*n(2)*n(3) + s*n(1), t*n(3)*n(3) + c]
        Q = reshape(m, [9])
    end function rotation

    function matmul9(a, b) result(c)
        !! Product of two column-major 3x3 value arrays.
        real(dp), intent(in) :: a(9), b(9)
        real(dp) :: c(9)
        c = reshape(matmul(reshape(a, [3, 3]), reshape(b, [3, 3])), [9])
    end function matmul9

    function diag9(d) result(m)
        real(dp), intent(in) :: d(3)
        real(dp) :: m(9)
        m = 0.0_dp
        m([1, 5, 9]) = d
    end function diag9

    function case_name(ic) result(name)
        !! Name of the deformation gradient `ic`.
        integer, intent(in) :: ic
        character(len=:), allocatable :: name
        select case (ic)
        case (1); name = 'identity'
        case (2); name = 'rotation z 30'
        case (3); name = 'rotation 111 120'
        case (4); name = 'rotation x 180'
        case (5); name = 'rotation 122 180'
        case (6); name = 'rotation y 1e-6 deg'
        case (7); name = 'rotation z 90'
        case (8); name = 'uniaxial 1.5'
        case (9); name = 'biaxial 1.2 1.2 0.7'
        case (10); name = 'simple shear 0.5'
        case (11); name = 'hydrostatic compression 0.5'
        case (12); name = 'strong compression 0.05'
        case (13); name = 'general distinct'
        case (14); name = 'rotated uniaxial'
        case (15); name = 'general large stretch rotated'
        case (16); name = 'nearly singular diagonal 1e-9'
        case (17); name = 'nearly singular general'
        case (18); name = 'small strain'
        case (19); name = 'reflection det < 0'
        case default; name = 'rotated biaxial 180'
        end select
    end function case_name

    function case_F(ic) result(F)
        !! Deformation gradient `ic` (column-major values F11, F21, F31, F12, ...).
        integer, intent(in) :: ic
        real(dp) :: F(9)
        select case (ic)
        case (1); F = diag9([1.0_dp, 1.0_dp, 1.0_dp])
        case (2); F = rotation([0.0_dp, 0.0_dp, 1.0_dp], 30.0_dp)
        case (3); F = rotation([1.0_dp, 1.0_dp, 1.0_dp], 120.0_dp)
        case (4); F = rotation([1.0_dp, 0.0_dp, 0.0_dp], 180.0_dp)
        case (5); F = rotation([1.0_dp, 2.0_dp, 2.0_dp], 180.0_dp)
        case (6); F = rotation([0.0_dp, 1.0_dp, 0.0_dp], 1.0e-6_dp)
        case (7); F = rotation([0.0_dp, 0.0_dp, 1.0_dp], 90.0_dp)
        case (8); F = diag9([1.5_dp, 1.0_dp, 1.0_dp])
        case (9); F = diag9([1.2_dp, 1.2_dp, 0.7_dp])
        case (10)
            F = diag9([1.0_dp, 1.0_dp, 1.0_dp])
            F(4) = 0.5_dp ! F12
        case (11); F = diag9([0.5_dp, 0.5_dp, 0.5_dp])
        case (12); F = diag9([0.05_dp, 1.1_dp, 1.05_dp])
        case (13); F = [1.1_dp, 0.05_dp, 0.15_dp, 0.2_dp, 0.9_dp, -0.2_dp, -0.1_dp, 0.3_dp, 1.2_dp]
        case (14)
            F = matmul9(rotation([1.0_dp, -1.0_dp, 0.5_dp], 75.0_dp), diag9([1.3_dp, 0.9_dp, 0.9_dp]))
        case (15)
            F = matmul9(rotation([0.3_dp, 0.4_dp, -1.0_dp], 140.0_dp), &
                [3.0_dp, 0.4_dp, -0.2_dp, 0.4_dp, 0.5_dp, 0.1_dp, -0.2_dp, 0.1_dp, 0.8_dp])
        case (16); F = diag9([1.0_dp, 1.0_dp, 1.0e-9_dp])
        case (17)
            ! Third column = first + second columns + 1e-8 e3: det = 1e-8
            F = [1.0_dp, 0.2_dp, 0.3_dp, 0.1_dp, 0.9_dp, -0.4_dp, 1.1_dp, 1.1_dp, -0.1_dp + 1.0e-8_dp]
        case (18); F = [1.002_dp, 0.0_dp, 0.0_dp, 0.004_dp, 0.997_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.001_dp]
        case (19); F = diag9([-1.0_dp, 1.2_dp, 0.9_dp])
        case default
            F = matmul9(rotation([0.0_dp, 1.0_dp, 1.0_dp], 180.0_dp), diag9([1.2_dp, 1.2_dp, 0.7_dp]))
        end select
    end function case_F

    function input_S() result(S)
        !! Symmetric stress applied to every transformation.
        type(ten_3D2Osym) :: S
        call S%init(xx=120.0_dp, yy=-35.0_dp, zz=60.0_dp, xy=25.0_dp, yz=-12.0_dp, xz=8.0_dp)
    end function input_S

    function input_S2() result(S)
        !! Second symmetric stress (optional argument of the tangent transformations).
        type(ten_3D2Osym) :: S
        call S%init(xx=-7.0_dp, yy=11.0_dp, zz=3.0_dp, xy=-2.0_dp, yz=5.0_dp, xz=13.0_dp)
    end function input_S2

    function input_P() result(P)
        !! Non-symmetric first Piola-Kirchhoff-like stress.
        type(ten_3D2O) :: P
        P%vals = [90.0_dp, -14.0_dp, 22.0_dp, 31.0_dp, -48.0_dp, 6.0_dp, -17.0_dp, 9.0_dp, 55.0_dp]
    end function input_P

    function input_C() result(C)
        !! 6x6 tangent with minor symmetries and no major symmetry.
        type(ten_3D4O2sym) :: C
        integer :: i, j
        do j = 1, 6
            do i = 1, 6
                C%vals(i, j) = real(7*i - 3*j, dp) + 0.25_dp*real(i*j, dp)
            end do
            C%vals(j, j) = C%vals(j, j) + 1000.0_dp
        end do
    end function input_C

    function input_R() result(R)
        !! Rotation passed to rotate_to/from_corotational.
        type(ten_3D2O) :: R
        R%vals = rotation([2.0_dp, -1.0_dp, 0.5_dp], 50.0_dp)
    end function input_R

    subroutine bits_equal(key, got, expected)
        !! Bitwise check of an alternative path against the recorded path.
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

    subroutine record_base(key, kin)
        !! Every binding of Base_kinematics through the polymorphic interface.
        character(len=*), intent(in) :: key
        class(Base_kinematics), intent(in) :: kin
        type(ten_3D2Osym) :: s, t
        type(ten_3D4O2sym) :: c, c2
        real(dp) :: lfin

        lfin = merge(1.0_dp, 0.0_dp, kin%is_finite_strain())
        s = kin%get_strain()
        call golden_check(key // ' | jacobian, is_finite_strain, get_strain', [kin%jacobian(), lfin, s%vals])
        s = kin%push_forward_stress(input_S())
        t = kin%pull_back_stress(input_S())
        call golden_check(key // ' | push_forward_stress, pull_back_stress', [s%vals, t%vals])
        s = kin%to_Cauchy(input_S())
        call golden_check(key // ' | to_Cauchy', s%vals)
        c = kin%push_forward_tangent(input_C())
        c2 = kin%push_forward_tangent(input_C(), input_S2())
        call bits_equal(key // ' push_forward_tangent with stress', reshape(c2%vals, [36]), reshape(c%vals, [36]))
        call golden_check(key // ' | push_forward_tangent', reshape(c%vals, [36]))
        c = kin%pull_back_tangent(input_C())
        c2 = kin%pull_back_tangent(input_C(), input_S2())
        call bits_equal(key // ' pull_back_tangent with stress', reshape(c2%vals, [36]), reshape(c%vals, [36]))
        call golden_check(key // ' | pull_back_tangent', reshape(c%vals, [36]))
        c = kin%to_spatial_tangent(input_C())
        c2 = kin%to_spatial_tangent(input_C(), input_S2())
        call bits_equal(key // ' to_spatial_tangent with stress', reshape(c2%vals, [36]), reshape(c%vals, [36]))
        call golden_check(key // ' | to_spatial_tangent', reshape(c%vals, [36]))
        s = kin%rotate_to_corotational(input_S(), input_R())
        t = kin%rotate_from_corotational(input_S(), input_R())
        call golden_check(key // ' | rotate_to/from_corotational', [s%vals, t%vals])
        s = kin%rotate_to_corotational(input_S())
        t = kin%rotate_from_corotational(input_S())
        call bits_equal(key // ' rotate_to_corotational without R', s%vals, input_S_vals())
        call bits_equal(key // ' rotate_from_corotational without R', t%vals, input_S_vals())
    end subroutine record_base

    function input_S_vals() result(v)
        real(dp) :: v(6)
        type(ten_3D2Osym) :: s
        s = input_S()
        v = s%vals
    end function input_S_vals

    subroutine record_finite(key, kin, full)
        !! Cached state, getters and conversions of Base_F_kinematics (if `full`), then
        !! `record_base`.
        character(len=*), intent(in) :: key
        class(Base_F_kinematics), intent(in) :: kin
        logical, intent(in) :: full
        type(ten_3D2Osym) :: a, b, c, d, e, f
        type(ten_3D2O) :: P, Q

        call golden_check(key // ' | state F, F_inv, R_rot, J, C_mat, b_spat, b_inv', [kin%F%vals, &
            kin%F_inv%vals, kin%R_rot%vals, kin%J, kin%C_mat%vals, kin%b_spat%vals, kin%b_inv%vals])
        if (full) then
            a = kin%get_RightCauchyGreen()
            b = kin%get_LeftCauchyGreen_inv()
            c = kin%get_GreenLagrange()
            d = kin%get_EulerAlmansi()
            e = kin%get_Material_Logarithmic()
            f = kin%get_Spatial_Logarithmic()
            P = kin%get_Rotation()
            call golden_check(key // ' | getters C, b_inv, E_GL, e_EA, H_mat, h_spat, R', [a%vals, b%vals, &
                c%vals, d%vals, e%vals, f%vals, P%vals])
            a = kin%convert_PK2_to_Cauchy(input_S())
            b = kin%convert_Cauchy_to_PK2(input_S())
            c = kin%convert_PK2_to_Kirchhoff(input_S())
            d = kin%convert_Kirchhoff_to_PK2(input_S())
            e = kin%convert_Cauchy_to_Kirchhoff(input_S())
            f = kin%convert_Kirchhoff_to_Cauchy(input_S())
            call golden_check(key // ' | convert PK2>Cauchy, Cauchy>PK2, PK2>Kirchhoff, Kirchhoff>PK2,' // &
                ' Cauchy>Kirchhoff, Kirchhoff>Cauchy', [a%vals, b%vals, c%vals, d%vals, e%vals, f%vals])
            P = kin%convert_PK2_to_PK1(input_S())
            a = kin%convert_PK1_to_PK2(input_P())
            Q = kin%convert_Cauchy_to_PK1(input_S())
            b = kin%convert_PK1_to_Cauchy(input_P())
            call golden_check(key // ' | convert PK2>PK1, PK1>PK2, Cauchy>PK1, PK1>Cauchy', [P%vals, a%vals, &
                Q%vals, b%vals])
        end if
        call record_base(key, kin)
    end subroutine record_finite

end module kin_golden_cases


program test_kin_golden
    !! Golden (characterization) test of the kinematics family, double precision, bit for bit against
    !! tests/reference/kinematics/kinematics_dp.txt (captured from the original OO code; see
    !! regenerate.sh there). Exit codes: those of `golden_finish` (1-4); 5 if an alternative path
    !! (optional stress argument, update generic vs specific binding, dt present or absent, R absent)
    !! differs bit for bit.
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2O, ten_3D2Osym
    use muscle_kinematics, only : Base_kinematics, Base_F_kinematics, Small_strain_kinematics, &
        Total_lagrangian_kinematics, Corotational_kinematics, Material_logarithmic_kinematics
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use kin_golden_cases
    implicit none

    class(Base_F_kinematics), allocatable :: fk, fk2
    class(Base_kinematics), allocatable :: bk
    type(Small_strain_kinematics) :: ss
    type(Corotational_kinematics) :: cr
    type(Total_lagrangian_kinematics) :: tl
    type(Material_logarithmic_kinematics) :: ml
    type(ten_3D2O) :: F, R
    type(ten_3D2Osym) :: eps
    integer :: ic, ik
    character(len=:), allocatable :: key

    call golden_init('kinematics/kinematics_dp.txt')

    ! --- Defaults of freshly declared objects (only the default-initialised components are defined) ---
    block
        type(Small_strain_kinematics) :: s0
        type(Corotational_kinematics) :: c0
        type(Total_lagrangian_kinematics) :: t0
        type(ten_3D2Osym) :: sig
        sig = c0%to_Cauchy(input_S())
        call golden_check('defaults | small jacobian, corotational J/has_rotation/jacobian/to_Cauchy, F jacobian', &
            [s0%jacobian(), c0%J, merge(1.0_dp, 0.0_dp, c0%has_rotation), c0%jacobian(), sig%vals, t0%jacobian()])
    end block

    ! --- Small strain: update by F and by strain ---
    do ic = 1, NCASES
        key = 'small_strain | ' // case_name(ic)
        F%vals = case_F(ic)
        call ss%update(F)
        call record_base(key, ss)
    end do
    call eps%init(xx=1.0e-3_dp, yy=-4.0e-4_dp, zz=2.0e-4_dp, xy=6.0e-4_dp, yz=-2.0e-4_dp, xz=1.0e-4_dp)
    call ss%update(eps)
    call record_base('small_strain | update by strain', ss)

    ! --- F-driven kinematics: update by F (generic through the polymorphic base) ---
    do ik = 1, NFINITE
        do ic = 1, NCASES
            key = kin_name(ik) // ' | ' // case_name(ic)
            call new_finite(ik, fk)
            F%vals = case_F(ic)
            call fk%update(F)
            ! Only the first kinematics records the getters and conversions of Base_F_kinematics (they
            ! do not depend on the concrete type; the state record checks the others share it).
            call record_finite(key, fk, ik == 1)
            ! Alternative path: the specific binding update_F gives the same object.
            call new_finite(ik, fk2)
            call fk2%update_F(F)
            call bits_equal(key // ' update_F', state_of(fk2), state_of(fk))
        end do
    end do

    ! --- update_F_base: updates the cached state but not the strain of the concrete type ---
    F%vals = case_F(13)
    call tl%update(F)
    F%vals = case_F(15)
    call tl%update_F_base(F)
    call record_finite('total_lagrangian | update_F_base after update', tl, .false.)
    F%vals = case_F(13)
    call ml%update(F)
    F%vals = case_F(14)
    call ml%update_F_base(F)
    call record_finite('material_logarithmic | update_F_base after update', ml, .false.)

    ! --- Incremental updates F_{n+1} = (drot + dstrain) F_n ---
    do ik = 1, NFINITE
        call run_sequence(ik, 'small rotations', 1)
        call run_sequence(ik, 'large rotations', 2)
        call run_sequence(ik, 'pure rotations', 3)
        call run_sequence(ik, 'from general F', 4)
    end do

    ! --- Corotational: strain and R of each F (H = ln U and R of the polar decomposition) ---
    do ic = 1, NCASES
        key = 'corotational | ' // case_name(ic)
        F%vals = case_F(ic)
        call ml%update(F)
        eps = ml%get_strain()
        R = ml%get_Rotation()
        call cr%update(eps, R)
        call golden_check(key // ' | state strain, R_rot, J, has_rotation', [cr%strain%vals, cr%R_rot%vals, &
            cr%J, merge(1.0_dp, 0.0_dp, cr%has_rotation)])
        call record_base(key, cr)
        ! Without R: the stored R_rot is kept but not used.
        call cr%update(eps)
        call golden_check(key // ' no R | state strain, J, has_rotation', [cr%strain%vals, cr%J, &
            merge(1.0_dp, 0.0_dp, cr%has_rotation)])
        call record_base(key // ' no R', cr)
    end do
    ! Through the polymorphic base (allocated as Corotational).
    allocate(bk, source=cr)
    call record_base('corotational | polymorphic copy', bk)

    call golden_finish()
    if (n_alt_bad > 0) stop 5

contains

    function state_of(k) result(v)
        !! Cached state and strain of an F-driven kinematics.
        class(Base_F_kinematics), intent(in) :: k
        real(dp) :: v(52)
        type(ten_3D2Osym) :: s
        s = k%get_strain()
        v = [k%F%vals, k%F_inv%vals, k%R_rot%vals, k%J, k%C_mat%vals, k%b_spat%vals, k%b_inv%vals, s%vals]
    end function state_of

    subroutine run_sequence(ik, name, iseq)
        !! Incremental sequence `iseq` on the F-driven kinematics `ik`: every step calls
        !! update(dstrain, drot, dt) (generic, polymorphic); a twin object without dt must be identical.
        integer, intent(in) :: ik, iseq
        character(len=*), intent(in) :: name
        type(ten_3D2O) :: drot, F0
        type(ten_3D2Osym) :: de
        integer :: istep, nsteps
        real(dp) :: ang
        real(dp), parameter :: ANG_LARGE(5) = [60.0_dp, 90.0_dp, 180.0_dp, 179.0_dp, 45.0_dp]
        real(dp), parameter :: ANG_PURE(5) = [30.0_dp, 90.0_dp, 180.0_dp, 1.0e-6_dp, 120.0_dp]
        character(len=:), allocatable :: kk
        character(len=16) :: buf

        call new_finite(ik, fk)
        call new_finite(ik, fk2)
        select case (iseq)
        case (4); F0%vals = case_F(13)
        case default; F0%vals = case_F(1)
        end select
        call fk%update(F0)
        call fk2%update(F0)
        nsteps = 5
        do istep = 1, nsteps
            select case (iseq)
            case (1)
                ang = 1.0_dp
                call de%init(xx=1.0e-3_dp, yy=-3.0e-4_dp, zz=-3.0e-4_dp, xy=2.0e-4_dp, yz=0.0_dp, xz=-1.0e-4_dp)
            case (2)
                ang = ANG_LARGE(istep)
                call de%init(xx=0.02_dp, yy=-0.01_dp, zz=0.005_dp, xy=0.01_dp, yz=-0.004_dp, xz=0.003_dp)
            case (3)
                ang = ANG_PURE(istep)
                de%vals = 0.0_dp
            case default
                ang = 15.0_dp*real(istep, dp)
                call de%init(xx=-0.03_dp, yy=0.02_dp, zz=0.01_dp, xy=-0.02_dp, yz=0.01_dp, xz=0.015_dp)
            end select
            drot%vals = rotation([1.0_dp, real(istep, dp), -0.5_dp], ang)
            call fk%update(de, drot, 0.1_dp)
            call fk2%update(de, drot)
            write(buf, '(i0)') istep
            kk = kin_name(ik) // ' | incremental ' // name // ' step ' // trim(buf)
            call bits_equal(kk // ' without dt', state_of(fk2), state_of(fk))
            call record_finite(kk, fk, .true.)
        end do
    end subroutine run_sequence

end program test_kin_golden
