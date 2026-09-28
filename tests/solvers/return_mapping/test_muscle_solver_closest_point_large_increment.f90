! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_large_increment
    implicit none
    logical :: passed

    call test_closest_point_cpb06_large_increment_1(passed)
    if (.not. passed) STOP 1

    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_large_increment

subroutine test_closest_point_cpb06_large_increment_1(passed)
    ! Single strain increments of 100 to 1000 yield strains from the virgin state,
    ! CPB06 for Ti-6Al-4V (Tuninetti et al. 2013, preprint Complas 2013, Table 2,
    ! row Wp = 48.66). Elasticity and Swift hardening are synthetic verification
    ! inputs, not a calibrated material.
    ! Every converged state is checked with quantities that do not depend on the
    ! solver iterations (Voigt order xx, yy, zz, xy, yz, xz, tensorial shear):
    !   yield consistency   |stress_eq(stress) - sigma_y(strain_pf)| / sigma_y0
    !   elastic law         stress = C : (strain - strain_p)
    !   associated flow     strain_p = strain_pf * dstressEq_dstress(stress)
    !   incompressibility   tr(strain_p) = 0
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_cpb06, only : CPB06
    use muscle_hard_swift, only : Swift_hardening
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    implicit none

    logical, intent(out) :: passed

    real(real64), parameter :: YOUNG = 70000D0, POISSON = 0.33D0
    real(real64), parameter :: SW_K = 500D0, SW_N = 0.2D0, SW_E0 = 0.01D0
    ! Flow residual accepted by Closest_point%iter (absolute, Euclidean norm of vals)
    real(real64), parameter :: TOL_FLOW = 1D-5
    real(real64), parameter :: TOL_YIELD = 1D-6
    real(real64), parameter :: TOL_ELASTIC = 1D-10
    ! The gradient may be the finite-difference fallback of Base_yield_critera,
    ! whose trace is zero only to about 1e-10 relative
    real(real64), parameter :: TOL_TRACE = 1D-8

    type(CPB06) :: cpb
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain, strain_p, stress, flow, law, zero
    real(real64) :: sigma_y0, eps_y, strain_pf, d(6,6), lambdas(3)
    real(real64) :: err_yield, err_elastic, err_flow, err_trace
    integer :: i, j, status, iters

    passed = .False.

    call cpb%init(c11=1.000D0, c12=-2.428D0, c13=-2.920D0, &
                  c21=-2.428D0, c22=1.652D0, c23=-2.236D0, &
                  c31=-2.920D0, c32=-2.236D0, c33=1.003D0, &
                  c44=-3.996D0, c55=-3.996D0, c66=-3.996D0, &
                  k=-0.165D0, a=2.0D0)
    call elas%set_parameters(young=YOUNG, poisson=POISSON)
    sw = Swift_hardening(k=SW_K, n=SW_N, e0=SW_E0)
    call solver%init(elasticity=elas, hardening=sw, yield=cpb)

    sigma_y0 = sw%stress(0D0)
    eps_y = sigma_y0/YOUNG

    ! Strain directions: uniaxial RD, uniaxial TD, equibiaxial, plane strain RD,
    ! shear xy, general 3D. Each is scaled to max|component| = lambda*eps_y.
    d(:,1) = [ 1.0D0, -0.5D0, -0.5D0, 0.0D0,  0.0D0, 0.0D0]
    d(:,2) = [-0.5D0,  1.0D0, -0.5D0, 0.0D0,  0.0D0, 0.0D0]
    d(:,3) = [ 1.0D0,  1.0D0, -2.0D0, 0.0D0,  0.0D0, 0.0D0]
    d(:,4) = [ 1.0D0,  0.0D0, -1.0D0, 0.0D0,  0.0D0, 0.0D0]
    d(:,5) = [ 0.0D0,  0.0D0,  0.0D0, 1.0D0,  0.0D0, 0.0D0]
    d(:,6) = [ 0.8D0, -0.3D0, -0.5D0, 0.4D0, -0.2D0, 0.3D0]
    lambdas = [100D0, 300D0, 1000D0]

    call zero%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)

    do j = 1, 6
        do i = 1, 3
            call strain%init(lambdas(i)*eps_y*d(:,j)/maxval(abs(d(:,j))))
            call history%init(strain_p=zero, strain_pf=0D0)
            call solver%solve(strain=strain, history=history, status=status, iters=iters)

            passed = status == STATUS_CONVERGED
            if (.not. passed) print*, "Not converged, direction", j, " lambda", lambdas(i), new_line('A'), &
                                      "Status:", status, " Iters:", iters
            if (.not. passed) return

            stress    = history%state_np1%stress
            strain_p  = history%state_np1%strain_p
            strain_pf = history%state_np1%strain_pf

            err_yield = abs(cpb%stress_eq(stress) - sw%stress(strain_pf))/sigma_y0
            law = elas%stress(strain - strain_p)
            err_elastic = maxval(abs(stress%vals - law%vals))/maxval(abs(law%vals))
            flow = strain_p - strain_pf*cpb%dstressEq_dstress(stress)
            err_flow = sqrt(sum(flow%vals**2))
            err_trace = abs(strain_p%xx() + strain_p%yy() + strain_p%zz())/sqrt(sum(strain_p%vals**2))

            passed = err_yield < TOL_YIELD .and. err_elastic < TOL_ELASTIC .and. &
                     err_flow < TOL_FLOW .and. err_trace < TOL_TRACE
            if (.not. passed) print*, "Converged state is not admissible, direction", j, &
                                      " lambda", lambdas(i), new_line('A'), &
                                      "Yield:", err_yield, " Elastic:", err_elastic, new_line('A'), &
                                      "Flow:", err_flow, " Trace:", err_trace
            if (.not. passed) return
        end do
    end do
end subroutine test_closest_point_cpb06_large_increment_1
