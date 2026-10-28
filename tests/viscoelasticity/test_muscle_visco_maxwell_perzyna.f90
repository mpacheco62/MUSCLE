! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_visco_maxwell_perzyna
    !! Generalized Maxwell network in series with the Perzyna return mapping, through the step cycle
    !! of `muscle_visco_base`: coefficients, freeze, solver%init, solver%solve, update(eps - eps_p),
    !! commit. Material in every case: von Mises, modified Voce, Perzyna (eta = 1, N = 2).
    implicit none
    logical :: passed

    call test_maxwell_perzyna_shear_history(passed)
    if (.not. passed) STOP 1

    call test_maxwell_perzyna_coupled_tangent(passed)
    if (.not. passed) STOP 2

    call test_maxwell_perzyna_update_elastic_strain(passed)
    if (.not. passed) STOP 3

    call test_maxwell_perzyna_rejected_step_retry(passed)
    if (.not. passed) STOP 4

    call test_maxwell_perzyna_no_branches_linear(passed)
    if (.not. passed) STOP 5

    print*, "Passed!", passed
    STOP 0
end program test_muscle_visco_maxwell_perzyna

subroutine test_maxwell_perzyna_shear_history(passed)
    !! Pure shear eps_xy, seven steps: elastic loading, two viscoplastic steps, a hold with flow at
    !! constant strain, a hold with relaxation only, elastic unloading and reverse flow.
    ! Why: every trial state carries the history stress of the branches, and flow continues during
    ! the hold (step 4) and restarts in reverse (step 7). Expected values: the same discrete model
    ! reduced to sigma_xy (sigma_eq = sqrt(3)|sigma_xy|), branch integral by quadrature and Perzyna
    ! root by findroot, Python with mpmath at 50 digits.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_yield_vonmises, only : VonMises
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point, STATUS_ELASTIC_CASE, STATUS_CONVERGED
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: vhist
    type(Plastic_material_history) :: phist
    type(Visco_coefficients) :: coefs
    type(Closest_point) :: solver
    type(Perzyna_overstress) :: pz
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain
    real(real64), parameter :: STRAIN_XY(7) = [0.02D0, 0.04D0, 0.06D0, 0.06D0, 0.06D0, 0.03D0, 0D0]
    real(real64), parameter :: DT(7) = [0.1D0, 0.1D0, 0.1D0, 0.5D0, 2D0, 0.1D0, 0.1D0]
    real(real64), parameter :: SIGMA_XY(7) = [23.211035215529665D0, 38.907355993392154D0, &
        47.603106882016906D0, 37.548329220635627D0, 30.822913658327813D0, -4.1072977399239949D0, &
        -36.840072584919045D0]
    real(real64), parameter :: STRAIN_P_XY(7) = [0D0, 0.0052750377377581620D0, 0.015905248422661713D0, &
        0.017114805445713135D0, 0.017114805445713135D0, 0.017114805445713135D0, 0.017025434718236899D0]
    real(real64), parameter :: STRAIN_PF(7) = [0D0, 0.0060910889157602187D0, 0.018365798916703221D0, &
        0.019762475062421103D0, 0.019762475062421103D0, 0.019762475062421103D0, 0.019865671489553257D0]
    integer, parameter :: EXPECTED_STATUS(7) = [STATUS_ELASTIC_CASE, STATUS_CONVERGED, STATUS_CONVERGED, &
        STATUS_CONVERGED, STATUS_ELASTIC_CASE, STATUS_ELASTIC_CASE, STATUS_CONVERGED]
    real(real64) :: sigma_xy_got(7), strain_p_xy_got(7), strain_pf_got(7)
    integer :: status(7), i

    passed = .true.
    call model%set_parameters(bulk=2000D0, shear_inf=300D0, shear=[200D0, 100D0], tau=[0.5D0, 5D0])
    pz%hard_law = Voce_modified_hardening(Sy=50D0, k=100D0, q=30D0, n=20D0)
    pz%eta = 1D0
    pz%n = 2D0
    call model%init_history(vhist)

    do i = 1, 7
        call strain%init(xx=0D0, yy=0D0, zz=0D0, xy=STRAIN_XY(i), yz=0D0, xz=0D0)
        coefs = model%coefficients(DT(i))
        call solver%init(elasticity=model%freeze(coefs, vhist), viscoplastic=pz, yield=vm, tol=1D-13)
        call solver%solve(strain=strain, history=phist, dt=DT(i), status=status(i))
        call model%update(coefs, strain - phist%state_np1%strain_p, vhist)
        call vhist%commit()
        call phist%commit()
        sigma_xy_got(i) = phist%state_n%stress%xy()
        strain_p_xy_got(i) = phist%state_n%strain_p%xy()
        strain_pf_got(i) = phist%state_n%strain_pf
    end do

    call check_true(passed, "shear history: elastic / converged pattern", all(status == EXPECTED_STATUS))
    ! Solver residuals below 1D-13 (1.9e-8 with the default 1D-8); Measured: 2.6e-15 or less
    call check_approx(passed, "shear history: sigma_xy", sigma_xy_got, SIGMA_XY, tol=1D-12)
    call check_approx(passed, "shear history: plastic strain xy", strain_p_xy_got, STRAIN_P_XY, tol=1D-12)
    call check_approx(passed, "shear history: equivalent plastic strain", strain_pf_got, STRAIN_PF, tol=1D-12)
end subroutine test_maxwell_perzyna_shear_history

subroutine test_maxwell_perzyna_coupled_tangent(passed)
    !! Consistent tangent of a viscoplastic step with branch history against the numerical derivative
    !! of the step map eps -> sigma at the frozen elasticity, with normal and shear components.
    ! Why: the tangent goes through Xi = [C_alg^-1 + dlambda ddf]^-1 with G_alg (not G_0), at a stress
    ! shifted by the history stress, and the previous step left plastic strain and branch stresses.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_yield_vonmises, only : VonMises
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point, STATUS_CONVERGED
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: vhist
    type(Plastic_material_history) :: phist
    type(Visco_coefficients) :: coefs
    type(Closest_point) :: solver
    type(Perzyna_overstress) :: pz
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64), parameter :: DT = 0.2D0
    integer :: status

    passed = .true.
    call model%set_parameters(bulk=2000D0, shear_inf=300D0, shear=[200D0, 100D0], tau=[0.5D0, 5D0])
    pz%hard_law = Voce_modified_hardening(Sy=50D0, k=100D0, q=30D0, n=20D0)
    pz%eta = 1D0
    pz%n = 2D0
    call model%init_history(vhist)

    ! Previous plastic step: leaves eps_p, ep and branch stresses
    call strain%init(xx=0.05D0, yy=-0.02D0, zz=-0.01D0, xy=0.02D0, yz=-0.01D0, xz=0.015D0)
    coefs = model%coefficients(DT)
    call solver%init(elasticity=model%freeze(coefs, vhist), viscoplastic=pz, yield=vm)
    call solver%solve(strain=strain, history=phist, dt=DT, status=status)
    call model%update(coefs, strain - phist%state_np1%strain_p, vhist)
    call vhist%commit()
    call phist%commit()

    ! Step under test: a different direction, so n is not the one of the previous step
    call strain%init(xx=0.07D0, yy=-0.01D0, zz=-0.03D0, xy=0.035D0, yz=0.01D0, xz=0.005D0)
    call solver%init(elasticity=model%freeze(coefs, vhist), viscoplastic=pz, yield=vm)
    call solver%solve(strain=strain, history=phist, dt=DT, status=status)
    call check_true(passed, "coupled tangent: plastic step converged", &
                    status == STATUS_CONVERGED .and. phist%state_np1%strain_pf > phist%state_n%strain_pf)

    call solver%tangent(strain=strain, history=phist, dt=DT, tangent=tangent)
    call solver%tangent_numerical(strain=strain, history=phist, dt=DT, tangent=numerical_tangent)
    ! Central differences with h = 1D-4*|eps| over solves converged to 1D-8 (as in the Perzyna solver
    ! test); Measured: 3.1e-9
    call check_approx(passed, "coupled tangent: analytical vs numerical", tangent, numerical_tangent, tol=1D-7)
end subroutine test_maxwell_perzyna_coupled_tangent

subroutine test_maxwell_perzyna_update_elastic_strain(passed)
    !! After a viscoplastic step the network is updated with eps - eps_p: the returned stress equals
    !! K tr(eps_e) I + 2 G_inf dev(eps_e) + sum_k s_k^(n+1), with the updated branch stresses.
    ! Why: updating with eps instead of eps - eps_p leaves this identity off by 2 G_k gamma_k dev(eps_p);
    ! the state has a volumetric part and the three shears, which case 1 does not load.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_yield_vonmises, only : VonMises
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: vhist
    type(Plastic_material_history) :: phist
    type(Visco_coefficients) :: coefs
    type(Closest_point) :: solver
    type(Perzyna_overstress) :: pz
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_e, vol, branches
    real(real64), parameter :: K = 2000D0, G_INF = 300D0, DT = 0.2D0
    real(real64) :: tr

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=G_INF, shear=[200D0, 100D0], tau=[0.5D0, 5D0])
    pz%hard_law = Voce_modified_hardening(Sy=50D0, k=100D0, q=30D0, n=20D0)
    pz%eta = 1D0
    pz%n = 2D0
    call model%init_history(vhist)

    call strain%init(xx=0.06D0, yy=-0.01D0, zz=-0.02D0, xy=0.03D0, yz=-0.02D0, xz=0.01D0)
    coefs = model%coefficients(DT)
    call solver%init(elasticity=model%freeze(coefs, vhist), viscoplastic=pz, yield=vm)
    call solver%solve(strain=strain, history=phist, dt=DT)
    call check_true(passed, "update eps_e: plastic step", phist%state_np1%strain_pf > 0D0)

    strain_e = strain - phist%state_np1%strain_p
    call model%update(coefs, strain_e, vhist)
    branches%vals = vhist%state_np1%h(:, 1) + vhist%state_np1%h(:, 2)
    tr = strain_e%xx() + strain_e%yy() + strain_e%zz()
    call vol%init(xx=K*tr, yy=K*tr, zz=K*tr, xy=0D0, yz=0D0, xz=0D0)
    ! Same affine map summed in another order; Measured: 3.5e-17
    call check_approx(passed, "update eps_e: returned stress = network stress", phist%state_np1%stress, &
                      vol + 2D0*G_INF*(.dev. strain_e) + branches, tol=1D-12)
end subroutine test_maxwell_perzyna_update_elastic_strain

subroutine test_maxwell_perzyna_rejected_step_retry(passed)
    !! A step rejected by the global iteration is discarded with rollback in both histories and
    !! retried with a smaller strain: the result equals the retried step taken directly.
    ! Why: the rejected attempt writes state_np1 of both histories (plastic strain and branches);
    ! a commit instead of a rollback, or a rollback that leaves state_np1, carries it into the retry.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_yield_vonmises, only : VonMises
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: vhist, vhist_ref
    type(Plastic_material_history) :: phist, phist_ref
    type(Visco_coefficients) :: coefs
    type(Closest_point) :: solver
    type(Perzyna_overstress) :: pz
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain
    real(real64), parameter :: DT = 0.2D0

    passed = .true.
    call model%set_parameters(bulk=2000D0, shear_inf=300D0, shear=[200D0, 100D0], tau=[0.5D0, 5D0])
    pz%hard_law = Voce_modified_hardening(Sy=50D0, k=100D0, q=30D0, n=20D0)
    pz%eta = 1D0
    pz%n = 2D0
    call model%init_history(vhist)
    coefs = model%coefficients(DT)

    ! Accepted plastic step
    call strain%init(xx=0.05D0, yy=-0.02D0, zz=-0.01D0, xy=0.02D0, yz=-0.01D0, xz=0.015D0)
    call solver%init(elasticity=model%freeze(coefs, vhist), viscoplastic=pz, yield=vm)
    call solver%solve(strain=strain, history=phist, dt=DT)
    call model%update(coefs, strain - phist%state_np1%strain_p, vhist)
    call vhist%commit()
    call phist%commit()
    vhist_ref = vhist
    phist_ref = phist

    ! Rejected attempt: three times the strain, then rollback of both histories
    call solver%init(elasticity=model%freeze(coefs, vhist), viscoplastic=pz, yield=vm)
    call solver%solve(strain=3D0*strain, history=phist, dt=DT)
    call model%update(coefs, 3D0*strain - phist%state_np1%strain_p, vhist)
    call vhist%rollback()
    call phist%rollback()
    call check_approx(passed, "rollback: branch stresses", vhist%state_np1%h, vhist%state_n%h, tol=0D0)
    call check_approx(passed, "rollback: plastic strain", phist%state_np1%strain_p, phist%state_n%strain_p, &
                      tol=0D0)

    ! Retry with 1.5 times the strain, on the rolled back histories and on the untouched copies
    call solver%init(elasticity=model%freeze(coefs, vhist), viscoplastic=pz, yield=vm)
    call solver%solve(strain=1.5D0*strain, history=phist, dt=DT)
    call model%update(coefs, 1.5D0*strain - phist%state_np1%strain_p, vhist)
    call vhist%commit()
    call phist%commit()
    call solver%solve(strain=1.5D0*strain, history=phist_ref, dt=DT)
    call model%update(coefs, 1.5D0*strain - phist_ref%state_np1%strain_p, vhist_ref)
    call vhist_ref%commit()
    call phist_ref%commit()
    ! Same operations on the same data; Measured: 0
    call check_approx(passed, "retry: stress", phist%state_n%stress, phist_ref%state_n%stress, tol=0D0)
    call check_approx(passed, "retry: branch stresses", vhist%state_n%h, vhist_ref%state_n%h, tol=0D0)
end subroutine test_maxwell_perzyna_rejected_step_retry

subroutine test_maxwell_perzyna_no_branches_linear(passed)
    !! With no branches the cycle is Perzyna with Elasticity_linear (K, G_inf): two viscoplastic
    !! steps with normal and shear components, stress, plastic strain and tangent.
    ! Why: freeze must give K_alg = K, G_alg = G_inf and a zero history stress, and update/commit
    ! with zero branches must not disturb the next step.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_hard_vocemod, only : Voce_modified_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_yield_vonmises, only : VonMises
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: vhist
    type(Plastic_material_history) :: phist, phist_lin
    type(Elasticity_linear) :: lin
    type(Visco_coefficients) :: coefs
    type(Closest_point) :: solver, solver_lin
    type(Perzyna_overstress) :: pz
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain
    type(ten_3D4O2sym) :: tangent, tangent_lin
    real(real64), parameter :: DT = 0.2D0, NONE(0) = [real(real64) ::]
    real(real64) :: load
    integer :: i

    passed = .true.
    call model%set_parameters_normalized(young=1000D0, poisson=0.3D0, g=NONE, tau=NONE)
    call lin%set_parameters(young=1000D0, poisson=0.3D0)
    pz%hard_law = Voce_modified_hardening(Sy=50D0, k=100D0, q=30D0, n=20D0)
    pz%eta = 1D0
    pz%n = 2D0
    call model%init_history(vhist)
    call solver_lin%init(elasticity=lin, viscoplastic=pz, yield=vm)

    call strain%init(xx=0.1D0, yy=-0.04D0, zz=-0.02D0, xy=0.05D0, yz=-0.02D0, xz=0.03D0)
    do i = 1, 2
        load = real(i, real64)
        coefs = model%coefficients(DT)
        call solver%init(elasticity=model%freeze(coefs, vhist), viscoplastic=pz, yield=vm)
        call solver%solve(strain=load*strain, history=phist, dt=DT)
        call model%update(coefs, load*strain - phist%state_np1%strain_p, vhist)
        call vhist%commit()
        call phist%commit()
        call solver_lin%solve(strain=load*strain, history=phist_lin, dt=DT)
        call phist_lin%commit()
    end do
    call solver%tangent(strain=2D0*strain, history=phist, dt=DT, tangent=tangent)
    call solver_lin%tangent(strain=2D0*strain, history=phist_lin, dt=DT, tangent=tangent_lin)

    call check_true(passed, "no branches: plastic steps", phist%state_n%strain_pf > 0D0)
    ! Lame constants from (K, G) instead of (E, nu): round-off only; Measured: 2.6e-16 (tangent: 0)
    call check_approx(passed, "no branches: stress", phist%state_n%stress, phist_lin%state_n%stress, tol=1D-12)
    call check_approx(passed, "no branches: plastic strain", phist%state_n%strain_p, &
                      phist_lin%state_n%strain_p, tol=1D-12)
    call check_approx(passed, "no branches: tangent", tangent, tangent_lin, tol=1D-12)
end subroutine test_maxwell_perzyna_no_branches_linear
