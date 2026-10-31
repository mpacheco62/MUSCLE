! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_math_golden_spectral_derivs
    !! Golden test of the public routines of `muscle_math_spectral_derivs` (double precision):
    !! `dEigenvalues_dTensor`, `d2Eigenvalues_dTensor2`, `spectral_decomposition`,
    !! `spectral_gradient` and `spectral_hessian`, on tensors with distinct, double and triple
    !! eigenvalues (exact and rotated), both branches of the decomposition (most separated
    !! eigenvalue first or last), the multiple of the identity, and the repeated-eigenvalue limit
    !! of the Hessian. `sym_outer` and `cross` are private (exercised only through
    !! `spectral_hessian` and `spectral_decomposition`). Results are compared bit for bit with
    !! `tests/reference/math/math_spectral_derivs_dp.txt`, captured from the original code
    !! (`tests/reference/math/regenerate.sh`).
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym, ten_3D4O3sym
    use muscle_math_operations, only : eigenvals
    use muscle_math_spectral_derivs, only : dEigenvalues_dTensor, d2Eigenvalues_dTensor2, &
        spectral_decomposition, spectral_gradient, spectral_hessian
    use muscle_test_helpers, only : park_miller, rnd, rotation_ref
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    implicit none

    integer, parameter :: NRANDOM = 5
    real(dp), parameter :: MAGNITUDE(NRANDOM) = [1.0_dp, 0.05_dp, 3.0_dp, 250.0_dp, 1.0_dp]
    type(park_miller) :: rng
    type(ten_3D2Osym) :: s
    real(dp) :: v(3,3), lam(3), dphi(3), d2phi(3,3)
    character(len=16) :: c
    integer :: icase, i, j

    call golden_init('math/math_spectral_derivs_dp.txt')
    rng = park_miller(1031)

    do icase = 1, NRANDOM
        write(c, '(" | case ", i0)') icase
        do i = 1, 6
            s%vals(i) = MAGNITUDE(icase)*rnd(rng)
        end do
        call all_routines('random', c, s)
    end do
    call all_routines('diag distinct', ' | case 1', sym(3.0_dp, -1.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp))
    call all_routines('plane', ' | case 1', sym(1.2_dp, -0.4_dp, 0.7_dp, 0.9_dp, 0.0_dp, 0.0_dp))
    call all_routines('double exact', ' | case 1', sym(2.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, 0.0_dp, 0.0_dp))
    call all_routines('double exact', ' | case 2', sym(3.0_dp, 3.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp))
    call all_routines('double exact', ' | case 3', sym(-1.0_dp, 3.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp))
    call all_routines('triple exact', ' | case 1', sym(2.5_dp, 2.5_dp, 2.5_dp, 0.0_dp, 0.0_dp, 0.0_dp))
    call all_routines('triple exact', ' | case 2', sym(0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp))
    call all_routines('middle = mean', ' | case 1', sym(0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp))
    call all_routines('rotated, first most separated', ' | case 1', rotated([5.0_dp, 1.0_dp, 0.0_dp]))
    call all_routines('rotated, last most separated', ' | case 1', rotated([1.0_dp, 0.0_dp, -5.0_dp]))
    call all_routines('rotated double', ' | case 1', rotated([4.0_dp, 4.0_dp, -2.0_dp]))
    call all_routines('rotated double', ' | case 2', rotated([4.0_dp, -2.0_dp, -2.0_dp]))
    call all_routines('rotated nearly double', ' | case 1', rotated([4.0_dp, 4.0_dp + 1.0e-7_dp, -2.0_dp]))
    call all_routines('rotated triple', ' | case 1', rotated([-7.0_dp, -7.0_dp, -7.0_dp]))

    ! spectral_gradient / spectral_hessian on arbitrary data (independent of the decomposition)
    do icase = 1, 3
        write(c, '(" | case ", i0)') icase
        v = rotation_ref(rng)
        if (icase == 3) then
            do j = 1, 3
                do i = 1, 3
                    v(i, j) = rnd(rng)
                end do
            end do
        end if
        call draw_phi(dphi, d2phi)
        select case (icase)
        case (1)
            lam = [3.0_dp, 1.0_dp, -2.0_dp]
        case (2)
            lam = [2.0_dp, 2.0_dp*(1.0_dp + 1.0e-12_dp), -1.0_dp]
        case (3)
            lam = [1.5_dp, 1.5_dp, 1.5_dp]
        end select
        call golden_check('input V arbitrary' // trim(c), reshape(v, [9]))
        call golden_check('input lam, dphi, d2phi arbitrary' // trim(c), [lam, dphi, reshape(d2phi, [9])])
        s = spectral_gradient(v, dphi)
        call golden_check('spectral_gradient arbitrary' // trim(c), s%vals)
        call emit_hessian('spectral_hessian arbitrary' // trim(c), spectral_hessian(lam, v, dphi, d2phi))
    end do
    call golden_finish()

contains

    function sym(t11, t22, t33, t12, t23, t13) result(t)
        !! Symmetric tensor from its components.
        real(dp), intent(in) :: t11, t22, t33, t12, t23, t13
            !! Components.
        type(ten_3D2Osym) :: t
            !! Tensor.
        t%vals = [t11, t22, t33, t12, t23, t13]
    end function sym

    function rotated(l) result(t)
        !! Symmetric tensor Q diag(l) Q^T with a pseudo-random rotation Q of `rng`.
        real(dp), intent(in) :: l(3)
            !! Eigenvalues.
        type(ten_3D2Osym) :: t
            !! Tensor.
        real(dp) :: q(3,3), m(3,3)
        integer :: i, j
        q = rotation_ref(rng)
        do j = 1, 3
            do i = 1, 3
                m(i, j) = q(i,1)*l(1)*q(j,1) + q(i,2)*l(2)*q(j,2) + q(i,3)*l(3)*q(j,3)
            end do
        end do
        t%vals = [m(1,1), m(2,2), m(3,3), m(1,2), m(2,3), m(1,3)]
    end function rotated

    subroutine draw_phi(dphi, d2phi)
        !! Random first derivatives and symmetric second derivatives of phi(lam).
        real(dp), intent(out) :: dphi(3)
            !! phi_a.
        real(dp), intent(out) :: d2phi(3,3)
            !! phi_ab (symmetric).
        integer :: i, j
        dphi = [rnd(rng), rnd(rng), rnd(rng)]
        do j = 1, 3
            do i = 1, j
                d2phi(i, j) = rnd(rng)
                d2phi(j, i) = d2phi(i, j)
            end do
        end do
    end subroutine draw_phi

    subroutine all_routines(name, c, t)
        !! Records the input and the results of every public routine for the tensor `t`.
        character(len=*), intent(in) :: name
            !! Name of the case family.
        character(len=*), intent(in) :: c
            !! Case suffix.
        type(ten_3D2Osym), intent(in) :: t
            !! Tensor.
        real(dp) :: ev(3), l(3), vv(3,3), dp1(3), dp2(3,3)
        type(ten_3D2Osym) :: h(3), g
        type(ten_3D4O3sym) :: h2(3)

        call golden_check('input ' // name // trim(c), t%vals)
        ev = eigenvals(t)
        call golden_check('eigenvals ' // name // trim(c), ev)
        h = dEigenvalues_dTensor(t, ev)
        call golden_check('dEigenvalues_dTensor ' // name // trim(c), [h(1)%vals, h(2)%vals, h(3)%vals])
        h2 = d2Eigenvalues_dTensor2(t, ev, h)
        call golden_check('d2Eigenvalues_dTensor2 ' // name // trim(c), [h2(1)%vals, h2(2)%vals, h2(3)%vals])
        call spectral_decomposition(t, l, vv)
        call golden_check('spectral_decomposition lam ' // name // trim(c), l)
        call golden_check('spectral_decomposition V ' // name // trim(c), reshape(vv, [9]))
        call draw_phi(dp1, dp2)
        call golden_check('input dphi, d2phi ' // name // trim(c), [dp1, reshape(dp2, [9])])
        g = spectral_gradient(vv, dp1)
        call golden_check('spectral_gradient ' // name // trim(c), g%vals)
        call emit_hessian('spectral_hessian ' // name // trim(c), spectral_hessian(l, vv, dp1, dp2))
    end subroutine all_routines

    subroutine emit_hessian(key, x)
        !! Records a `ten_3D4O3sym` result.
        character(len=*), intent(in) :: key
            !! Key.
        type(ten_3D4O3sym), intent(in) :: x
            !! Result.
        call golden_check(key, x%vals)
    end subroutine emit_hessian

end program test_math_golden_spectral_derivs
