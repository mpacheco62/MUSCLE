! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module math_golden_derivatives_funcs
    !! Deterministic, non-trivial test functions of the golden test of `muscle_math_derivatives`
    !! (module procedures, not internal ones). Every component of their first (and, for the
    !! scalar functions of a symmetric tensor, second) derivatives is non-zero at generic points.
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2O, ten_3D2Osym
    implicit none
    private

    public :: fss_1, fss_2, f2o_1, f2o_2, fsym_1, fsym_2, gsym_1, gsym_2

    real(dp), parameter :: W6(6) = [0.9_dp, -0.7_dp, 0.5_dp, 1.3_dp, -1.1_dp, 0.8_dp]
        !! Weights of the linear forms of the symmetric tensor functions.

contains

    pure function fss_1(x) result(res)
        !! exp(0.7 x) sin(x) + x^3 - 2 x.
        real(dp), intent(in) :: x
            !! Argument.
        real(dp) :: res
            !! Value.
        res = exp(0.7_dp*x)*sin(x) + x**3 - 2.0_dp*x
    end function fss_1

    pure function fss_2(x) result(res)
        !! log(1 + x^2) + 1/(2 + 0.5 x) + cos(3 x).
        real(dp), intent(in) :: x
            !! Argument.
        real(dp) :: res
            !! Value.
        res = log(1.0_dp + x*x) + 1.0_dp/(2.0_dp + 0.5_dp*x) + cos(3.0_dp*x)
    end function fss_2

    pure function f2o_1(t) result(res)
        !! Weighted squares, a cubic term, a sine and an exponential of a general 3D tensor.
        type(ten_3D2O), intent(in) :: t
            !! Argument.
        real(dp) :: res
            !! Value.
        integer :: i
        res = 0.0_dp
        do i = 1, 9
            res = res + (1.0_dp + 0.1_dp*i)*t%vals(i)**2
        end do
        res = res + t%vals(1)*t%vals(5)*t%vals(9) + sin(t%vals(2)*t%vals(4) + t%vals(3)) &
            + exp(0.2_dp*(t%vals(3) + t%vals(6) + t%vals(7) - t%vals(8)))
    end function f2o_1

    pure function f2o_2(t) result(res)
        !! sqrt(1 + (w . t)^2) + 0.1 (sum t)^3 of a general 3D tensor.
        type(ten_3D2O), intent(in) :: t
            !! Argument.
        real(dp) :: res
            !! Value.
        real(dp) :: l
        integer :: i
        l = 0.0_dp
        do i = 1, 9
            l = l + (0.5_dp + 0.15_dp*i)*t%vals(i)
        end do
        res = sqrt(1.0_dp + l*l) + 0.1_dp*sum(t%vals)**3
    end function f2o_2

    pure function fsym_1(t) result(res)
        !! exp(0.3 w . t) + t11 t22 t33 + t12 t23 t13 + weighted squares of a symmetric tensor.
        type(ten_3D2Osym), intent(in) :: t
            !! Argument.
        real(dp) :: res
            !! Value.
        integer :: i
        res = exp(0.3_dp*sum(W6*t%vals)) + t%vals(1)*t%vals(2)*t%vals(3) + t%vals(4)*t%vals(5)*t%vals(6)
        do i = 1, 6
            res = res + 0.25_dp*i*t%vals(i)**2
        end do
    end function fsym_1

    pure function fsym_2(t) result(res)
        !! sqrt(1 + J2) + 0.1 I1^3 + sin(t12 + 2 t23 - t13) of a symmetric tensor.
        type(ten_3D2Osym), intent(in) :: t
            !! Argument.
        real(dp) :: res
            !! Value.
        real(dp) :: j2, i1
        i1 = t%vals(1) + t%vals(2) + t%vals(3)
        j2 = ((t%vals(1) - t%vals(2))**2 + (t%vals(2) - t%vals(3))**2 + (t%vals(3) - t%vals(1))**2)/6.0_dp &
            + t%vals(4)**2 + t%vals(5)**2 + t%vals(6)**2
        res = sqrt(1.0_dp + j2) + 0.1_dp*i1**3 + sin(t%vals(4) + 2.0_dp*t%vals(5) - t%vals(6))
    end function fsym_2

    pure function gsym_1(t) result(res)
        !! Symmetric tensor of components exp(0.2 W_i . t) + 0.1 i t_i^2 (W_ij = 1 + 0.1 (i - j)).
        type(ten_3D2Osym), intent(in) :: t
            !! Argument.
        type(ten_3D2Osym) :: res
            !! Value.
        integer :: i, j
        real(dp) :: l
        do i = 1, 6
            l = 0.0_dp
            do j = 1, 6
                l = l + (1.0_dp + 0.1_dp*(i - j))*t%vals(j)
            end do
            res%vals(i) = exp(0.2_dp*l) + 0.1_dp*i*t%vals(i)**2
        end do
    end function gsym_1

    pure function gsym_2(t) result(res)
        !! Symmetric tensor of components t_i |t|^2 + sin(t_(i mod 6 + 1)).
        type(ten_3D2Osym), intent(in) :: t
            !! Argument.
        type(ten_3D2Osym) :: res
            !! Value.
        integer :: i
        do i = 1, 6
            res%vals(i) = t%vals(i)*sum(t%vals**2) + sin(t%vals(mod(i, 6) + 1))
        end do
    end function gsym_2

end module math_golden_derivatives_funcs

program test_math_golden_derivatives
    !! Golden test of every public routine of `muscle_math_derivatives` (double precision):
    !! `derivative` (scalar-scalar, scalar-ten_3D2O, scalar-ten_3D2Osym, ten_3D2Osym-ten_3D2Osym),
    !! `derivative2O` (scalar-ten_3D2Osym) and `dsquare`, through the public generic names only,
    !! with the default step and an explicit `eps`. Results are compared bit for bit with
    !! `tests/reference/math/math_derivatives_dp.txt`, captured from the original code
    !! (`tests/reference/math/regenerate.sh`).
    !!
    !! Cases: arguments of magnitude ~1, ~0.05 and ~3 drawn with `park_miller`, and the zero
    !! argument (the default step falls back to its absolute minimum).
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2O, ten_3D2Osym, ten_3D4O2sym, ten_3D4O3sym
    use muscle_math_derivatives, only : derivative, derivative2O, dsquare
    use muscle_test_helpers, only : park_miller, rnd
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use math_golden_derivatives_funcs
    implicit none

    integer, parameter :: NCASES = 4
    real(dp), parameter :: MAGNITUDE(NCASES) = [1.0_dp, 0.05_dp, 3.0_dp, 0.0_dp]
    real(dp), parameter :: EPS1 = 1.0e-5_dp
        !! Explicit step of the first derivatives.
    real(dp), parameter :: EPS2 = 1.0e-3_dp
        !! Explicit step of the second derivative and of the tensor-tensor derivative.
    type(park_miller) :: rng
    type(ten_3D2O) :: a
    type(ten_3D2Osym) :: s
    type(ten_3D4O2sym) :: d42
    type(ten_3D4O3sym) :: d43
    real(dp) :: x
    character(len=16) :: c
    integer :: icase, i

    call golden_init('math/math_derivatives_dp.txt')
    rng = park_miller(20261031)
    do icase = 1, NCASES
        write(c, '(" | case ", i0)') icase

        x = MAGNITUDE(icase)*rnd(rng)
        call golden_check('input scalar' // trim(c), [x])
        call golden_check('derivative(fss_1, x)' // trim(c), [derivative(fss_1, x)])
        call golden_check('derivative(fss_1, x, eps)' // trim(c), [derivative(fss_1, x, EPS1)])
        call golden_check('derivative(fss_2, x)' // trim(c), [derivative(fss_2, x)])
        call golden_check('derivative(fss_2, x, eps)' // trim(c), [derivative(fss_2, x, eps=EPS1)])

        do i = 1, 9
            a%vals(i) = MAGNITUDE(icase)*rnd(rng)
        end do
        call golden_check('input ten_3D2O' // trim(c), a%vals)
        call emit_2o('derivative(f2o_1, ten_3D2O)' // trim(c), derivative(f2o_1, a))
        call emit_2o('derivative(f2o_1, ten_3D2O, eps)' // trim(c), derivative(f2o_1, a, EPS1))
        call emit_2o('derivative(f2o_2, ten_3D2O)' // trim(c), derivative(f2o_2, a))
        call emit_2o('derivative(f2o_2, ten_3D2O, eps)' // trim(c), derivative(f2o_2, a, eps=EPS1))

        do i = 1, 6
            s%vals(i) = MAGNITUDE(icase)*rnd(rng)
        end do
        call golden_check('input ten_3D2Osym' // trim(c), s%vals)
        call emit_2osym('derivative(fsym_1, ten_3D2Osym)' // trim(c), derivative(fsym_1, s))
        call emit_2osym('derivative(fsym_1, ten_3D2Osym, eps)' // trim(c), derivative(fsym_1, s, EPS1))
        call emit_2osym('derivative(fsym_2, ten_3D2Osym)' // trim(c), derivative(fsym_2, s))
        call emit_2osym('derivative(fsym_2, ten_3D2Osym, eps)' // trim(c), derivative(fsym_2, s, eps=EPS1))

        d43 = derivative2O(fsym_1, s)
        call golden_check('derivative2O(fsym_1, ten_3D2Osym)' // trim(c), d43%vals)
        d43 = derivative2O(fsym_1, s, EPS2)
        call golden_check('derivative2O(fsym_1, ten_3D2Osym, eps)' // trim(c), d43%vals)
        d43 = derivative2O(fsym_2, s)
        call golden_check('derivative2O(fsym_2, ten_3D2Osym)' // trim(c), d43%vals)
        d43 = derivative2O(fsym_2, s, eps=EPS2)
        call golden_check('derivative2O(fsym_2, ten_3D2Osym, eps)' // trim(c), d43%vals)

        d42 = derivative(gsym_1, s)
        call golden_check('derivative(gsym_1, ten_3D2Osym)' // trim(c), reshape(d42%vals, [36]))
        d42 = derivative(gsym_1, s, EPS2)
        call golden_check('derivative(gsym_1, ten_3D2Osym, eps)' // trim(c), reshape(d42%vals, [36]))
        d42 = derivative(gsym_2, s)
        call golden_check('derivative(gsym_2, ten_3D2Osym)' // trim(c), reshape(d42%vals, [36]))
        d42 = derivative(gsym_2, s, eps=EPS2)
        call golden_check('derivative(gsym_2, ten_3D2Osym, eps)' // trim(c), reshape(d42%vals, [36]))

        d43 = dsquare(s)
        call golden_check('dsquare(ten_3D2Osym)' // trim(c), d43%vals)
    end do
    call golden_finish()

contains

    subroutine emit_2o(key, t)
        !! Records a `ten_3D2O` result (exact dummy type: the generic must resolve to it).
        character(len=*), intent(in) :: key
            !! Key.
        type(ten_3D2O), intent(in) :: t
            !! Result.
        call golden_check(key, t%vals)
    end subroutine emit_2o

    subroutine emit_2osym(key, t)
        !! Records a `ten_3D2Osym` result (exact dummy type: the generic must resolve to it).
        character(len=*), intent(in) :: key
            !! Key.
        type(ten_3D2Osym), intent(in) :: t
            !! Result.
        call golden_check(key, t%vals)
    end subroutine emit_2osym

end program test_math_golden_derivatives
