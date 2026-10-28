program test_muscle_visco_maxwell
    implicit none
    logical :: passed

    call test_relaxation(passed)
    if (.not. passed) STOP 1

    call test_constant_rate(passed)
    if (.not. passed) STOP 2

    call test_time_step_limits(passed)
    if (.not. passed) STOP 3

    call test_freeze_update_consistency(passed)
    if (.not. passed) STOP 4

    call test_history(passed)
    if (.not. passed) STOP 5

    call test_closest_point_coupling(passed)
    if (.not. passed) STOP 6

    call test_maxwell_fea_wrapper_cycle(passed)
    if (.not. passed) STOP 7

    call test_maxwell_frozen_tangents_2d_and_dev(passed)
    if (.not. passed) STOP 8

    call test_maxwell_gamma_series_small_x(passed)
    if (.not. passed) STOP 9

    call test_maxwell_init_history_strain_e(passed)
    if (.not. passed) STOP 10

    call test_maxwell_normalized_sum_g_one(passed)
    if (.not. passed) STOP 11

    print*, "Passed!", passed
end program test_muscle_visco_maxwell

subroutine test_relaxation(passed)
    !! Step strain and hold: sigma(t) = K tr(eps) I + 2 G(t) dev(eps), exact for any dt
    !! (only the first step, which applies the strain, has to be short).
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(ten_3D2Osym) :: strain, stress, expected, vol
    real(real64), parameter :: K = 2000D0, DT1 = 1D-12
    real(real64), parameter :: DTS(5) = [0.05D0, 0.3D0, 0.7D0, 2D0, 5D0]
    real(real64) :: t, tr
    integer :: i
    character(len=40) :: label

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=100D0, shear=[200D0, 300D0], tau=[0.1D0, 1D0])
    call model%init_history(hist)

    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    tr = strain%xx() + strain%yy() + strain%zz()
    call vol%init(xx=K*tr, yy=K*tr, zz=K*tr, xy=0D0, yz=0D0, xz=0D0)

    ! Apply the strain in a very short step
    t = DT1
    elas = model%freeze_dt(DT1, hist)
    stress = elas%stress(strain)
    call model%update_dt(DT1, strain, hist)
    call hist%commit()
    expected = vol + 2D0*model%shear_instantaneous()*(.dev. strain)
    call check_approx(passed, "relaxation: instantaneous response", stress, expected, tol=1D-9)

    ! Hold the strain with coarse steps
    do i = 1, size(DTS)
        t = t + DTS(i)
        elas = model%freeze_dt(DTS(i), hist)
        stress = elas%stress(strain)
        call model%update_dt(DTS(i), strain, hist)
        call hist%commit()
        expected = vol + 2D0*model%relaxation_shear(t)*(.dev. strain)
        write(label, "(A, ES9.2)") "relaxation: stress at t =", t
        call check_approx(passed, trim(label), stress, expected, tol=1D-9)
    end do
end subroutine test_relaxation

subroutine test_constant_rate(passed)
    !! Constant strain rate D: the update is exact for any dt,
    !! sigma(t) = K tr(D) t I + 2 dev(D) [G_inf t + sum G_k tau_k (1 - exp(-t/tau_k))].
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: rate, strain, stress, expected, vol
    real(real64), parameter :: K = 2000D0, G_INF = 100D0
    real(real64), parameter :: G(2) = [200D0, 300D0], TAU(2) = [0.1D0, 1D0]
    real(real64), parameter :: DTS(5) = [0.03D0, 0.2D0, 0.5D0, 1.5D0, 3D0]
    real(real64) :: t, tr, gt
    integer :: i
    character(len=40) :: label

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=G_INF, shear=G, tau=TAU)
    call model%init_history(hist)

    call rate%init(xx=2D-3, yy=-1D-3, zz=5D-4, xy=1D-3, yz=3D-4, xz=-2D-4)
    tr = rate%xx() + rate%yy() + rate%zz()

    t = 0D0
    do i = 1, size(DTS)
        t = t + DTS(i)
        strain = t*rate
        coefs = model%coefficients(DTS(i))
        elas = model%freeze(coefs, hist)
        stress = elas%stress(strain)
        call model%update(coefs, strain, hist)
        call hist%commit()

        gt = G_INF*t + sum(G*TAU*(1D0 - exp(-t/TAU)))
        call vol%init(xx=K*tr*t, yy=K*tr*t, zz=K*tr*t, xy=0D0, yz=0D0, xz=0D0)
        expected = vol + 2D0*gt*(.dev. rate)
        write(label, "(A, ES9.2)") "constant rate: stress at t =", t
        call check_approx(passed, trim(label), stress, expected, tol=1D-10)
    end do
end subroutine test_constant_rate

subroutine test_time_step_limits(passed)
    !! dt -> 0 gives the instantaneous (linear elastic) response, dt -> inf the long-term one.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Elasticity_linear) :: lin
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain
    type(ten_3D4O3sym) :: tan_frozen, tan_lin
    real(real64), parameter :: E = 1000D0, NU = 0.3D0
    real(real64) :: g0

    passed = .true.
    call model%set_parameters_normalized(young=E, poisson=NU, g=[0.2D0, 0.3D0], tau=[0.1D0, 1D0])
    call model%init_history(hist)
    call lin%set_parameters(young=E, poisson=NU)
    g0 = 0.5D0*E/(1D0 + NU)

    call check_approx(passed, "normalized: G_0", model%shear_instantaneous(), g0, tol=1D-12)
    call check_approx(passed, "normalized: G(inf)", model%relaxation_shear(1D12), 0.5D0*g0, tol=1D-12)

    coefs = model%coefficients(0D0)
    call check_approx(passed, "dt = 0: G_alg = G_0", coefs%g_alg, g0, tol=1D-12)
    call check_approx(passed, "dt = 0: K_alg", coefs%k_alg, E/(3D0*(1D0 - 2D0*NU)), tol=1D-12)

    ! dt = 0 without history must be exactly the linear elastic law
    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    elas = model%freeze(coefs, hist)
    call check_approx(passed, "dt = 0: stress = linear elastic", elas%stress(strain), &
                      lin%stress(strain), tol=1D-12)
    tan_frozen = elas%dstress_dstrain(strain)
    tan_lin = lin%dstress_dstrain(strain)
    call check_approx(passed, "dt = 0: tangent = linear elastic", tan_frozen, tan_lin, tol=1D-12)

    coefs = model%coefficients(1D12)
    call check_approx(passed, "dt -> inf: G_alg = G_inf", coefs%g_alg, 0.5D0*g0, tol=1D-9)
end subroutine test_time_step_limits

subroutine test_freeze_update_consistency(passed)
    !! The frozen stress at eps_e_n+1 equals K tr(eps) I + 2 G_inf dev(eps) + sum s_k^(n+1)
    !! after the update. Also: analytical vs numerical tangent, deviatoric part and plane strain.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_math_derivatives
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain, branches, expected, vol
    type(ten_2D2Osym) :: strain_2d, stress_2d, expected_2d
    type(ten_3D4O2sym) :: tan_num, tan_an
    real(real64), parameter :: K = 2000D0, G_INF = 100D0
    real(real64) :: tr
    integer :: i

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=G_INF, shear=[200D0, 300D0, 50D0], &
                              tau=[0.1D0, 1D0, 10D0])
    call model%init_history(hist)

    ! Build a non-trivial (plane strain) history
    call strain%init(xx=1D-3, yy=-4D-4, zz=0D0, xy=6D-4, yz=0D0, xz=0D0)
    call model%update_dt(0.2D0, strain, hist)
    call hist%commit()
    call strain%init(xx=4D-4, yy=2D-4, zz=0D0, xy=-3D-4, yz=0D0, xz=0D0)
    call model%update_dt(0.5D0, strain, hist)
    call hist%commit()

    call strain%init(xx=7D-4, yy=-1D-4, zz=0D0, xy=2D-4, yz=0D0, xz=0D0)
    coefs = model%coefficients(0.3D0)
    elas = model%freeze(coefs, hist)
    call model%update(coefs, strain, hist)

    branches%vals = 0D0
    do i = 1, model%n_branches()
        branches%vals = branches%vals + hist%state_np1%h(:, i)
    end do
    tr = strain%xx() + strain%yy() + strain%zz()
    call vol%init(xx=K*tr, yy=K*tr, zz=K*tr, xy=0D0, yz=0D0, xz=0D0)
    expected = vol + 2D0*G_INF*(.dev. strain) + branches
    call check_approx(passed, "freeze(eps_n+1) = K tr I + 2 G_inf e + sum s_k^(n+1)", &
                      elas%stress(strain), expected, tol=1D-12)

    call check_approx(passed, "stress_dev = dev(stress)", elas%stress_dev(strain), &
                      .dev. elas%stress(strain), tol=1D-12)

    tan_an = elas%dstress_dstrain(strain)
    tan_num = derivative(wrapper, strain)
    call check_approx(passed, "frozen tangent: analytical vs numerical", tan_an, tan_num, tol=1D-7)

    call strain_2d%init(xx=strain%xx(), yy=strain%yy(), zz=strain%zz(), xy=strain%xy())
    stress_2d = elas%stress(strain_2d)
    call expected_2d%init(xx=expected%xx(), yy=expected%yy(), zz=expected%zz(), xy=expected%xy())
    call check_approx(passed, "plane strain stress = 3D components", stress_2d, expected_2d, tol=1D-12)

    contains
    pure function wrapper(x) result(res)
        type(ten_3D2Osym), intent(in) :: x
        type(ten_3D2Osym) :: res
        res = elas%stress(x)
    end function wrapper
end subroutine test_freeze_update_consistency

subroutine test_history(passed)
    !! commit / rollback and FEA pack / unpack round trip.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist, hist2
    type(ten_3D2Osym) :: strain
    real(real64) :: hsv(18), h_n(6, 2)

    passed = .true.
    call model%set_parameters(bulk=2000D0, shear_inf=100D0, shear=[200D0, 300D0], tau=[0.1D0, 1D0])
    call model%init_history(hist)
    call check_true(passed, "init: n_branches", hist%state_n%n_branches == 2)

    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    call model%update_dt(0.1D0, strain, hist)
    call hist%commit()
    h_n = hist%state_n%h(:, 1:2)

    ! A failed global iteration must not touch state_n
    call model%update_dt(0.1D0, 2D0*strain, hist)
    call hist%rollback()
    call check_approx(passed, "rollback: history tensors", hist%state_np1%h(:, 1:2), h_n, tol=0D0)
    call check_approx(passed, "rollback: elastic strain", hist%state_np1%strain_e, strain, tol=0D0)

    call hist%pack_to_fea(hsv)
    call model%init_history(hist2)
    call hist2%unpack_from_fea(hsv)
    call check_approx(passed, "pack/unpack: history tensors", hist2%state_n%h(:, 1:2), h_n, tol=0D0)
    call check_approx(passed, "pack/unpack: elastic strain", hist2%state_n%strain_e, strain, tol=0D0)
end subroutine test_history

subroutine test_closest_point_coupling(passed)
    !! The frozen elasticity (with a non-zero history stress) plugged into the rate-independent
    !! closest point: the return mapping converges, f = 0, its consistent tangent matches the
    !! numerical one, and the viscoelastic update at eps - eps_p reproduces the returned stress.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: vhist
    type(Plastic_material_history) :: phist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(Swift_hardening) :: sw
    type(ten_3D2Osym) :: strain, strain_e, strain_p0, stress_hist
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: f
    integer :: status, iters

    passed = .true.
    call model%set_parameters_normalized(young=1000D0, poisson=0.3D0, g=[0.3D0, 0.2D0], &
                                         tau=[0.1D0, 1D0])
    call model%init_history(vhist)

    ! Previous step leaves a non-zero viscous history
    call strain%init(xx=2D-2, yy=-5D-3, zz=-5D-3, xy=4D-3, yz=0D0, xz=1D-3)
    call model%update_dt(0.05D0, strain, vhist)
    call vhist%commit()

    coefs = model%coefficients(0.05D0)
    elas = model%freeze(coefs, vhist)
    stress_hist = elas%history_stress()
    call check_true(passed, "coupling: history stress is non-zero", stress_hist%norm() > 0D0)

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    call strain_p0%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call phist%init(strain_p=strain_p0, strain_pf=0D0)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)

    call strain%init(xx=0.15D0, yy=-0.05D0, zz=-0.05D0, xy=1D-2, yz=0D0, xz=0D0)
    call solver%solve(strain=strain, history=phist, status=status, iters=iters)
    call check_true(passed, "coupling: return mapping converged", status == STATUS_CONVERGED)
    call check_true(passed, "coupling: plastic step", phist%state_np1%strain_pf > 0D0)

    f = vm%stress_eq(phist%state_np1%stress) - sw%stress(phist%state_np1%strain_pf)
    call check_approx(passed, "coupling: f = 0", f, 0D0, tol=1D-5)

    call solver%tangent(strain=strain, history=phist, tangent=tangent)
    call solver%tangent_numerical(strain=strain, history=phist, tangent=numerical_tangent)
    call check_approx(passed, "coupling: consistent tangent", tangent, numerical_tangent, tol=1D-6)

    ! Viscoelastic update with the converged elastic strain
    strain_e = strain - phist%state_np1%strain_p
    call model%update(coefs, strain_e, vhist)
    call check_approx(passed, "coupling: frozen stress = returned stress", elas%stress(strain_e), &
                      phist%state_np1%stress, tol=1D-12)
    call check_true(passed, "coupling: update stores eps_e", &
                    vhist%state_np1%strain_e%is_approx(strain_e, tol=0D0))
end subroutine test_closest_point_coupling

subroutine test_maxwell_fea_wrapper_cycle(passed)
    !! Cycle of an FEA wrapper (UMAT/USERMAT style): in every call the state is a new local
    !! object rebuilt from the flat history array -- init_history, unpack_from_fea, freeze,
    !! update, pack_to_fea -- and the array is the only thing that survives between calls.
    !! Reference: the same three steps with a single state kept in memory.
    !! Hard case: after the update the trial state differs from the converged one, so packing
    !! state_n instead of state_np1 makes the next step start from a stale history.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist_mem, hist, hist_new
    type(Elasticity_frozen) :: elas, elas_new
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: eps(3), strain, stress_mem, stress_fea
    real(real64), parameter :: K = 2000D0, G_INF = 100D0
    real(real64), parameter :: DTS(3) = [0.05D0, 0.3D0, 2D0]
    real(real64) :: hsv(18), hsv_out(18)
    integer :: step
    character(len=40) :: label

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=G_INF, shear=[200D0, 300D0], tau=[0.1D0, 1D0])

    ! Three steps with normals and shears, and a history that is not small
    call eps(1)%init(xx=1D-3, yy=-4D-4, zz=2D-4, xy=6D-4, yz=-2D-4, xz=1D-4)
    call eps(2)%init(xx=2D-3, yy=-6D-4, zz=1D-4, xy=3D-4, yz=1D-4, xz=-5D-4)
    call eps(3)%init(xx=1D-3, yy=-1D-4, zz=5D-4, xy=-8D-4, yz=2D-4, xz=3D-4)

    hsv = 0D0                        ! undeformed point with no viscous history
    call model%init_history(hist_mem)

    do step = 1, size(DTS)
        strain = eps(step)
        coefs = model%coefficients(DTS(step))

        ! Reference: one state kept in memory
        elas = model%freeze(coefs, hist_mem)
        stress_mem = elas%stress(strain)
        call model%update(coefs, strain, hist_mem)
        call hist_mem%commit()

        ! Wrapper: new local state on every call, all the information in the FEA array
        call model%init_history(hist)
        call hist%unpack_from_fea(hsv)
        elas = model%freeze(coefs, hist)
        stress_fea = elas%stress(strain)
        call model%update(coefs, strain, hist)
        call hist%pack_to_fea(hsv)   ! trial state back to the array, then global convergence
        call hist%commit()

        write(label, "(A, I0)") "FEA cycle: stress, step ", step
        ! Same operations on the same numbers: the wrapper must reproduce the state in memory
        ! exactly. Measured: 0 (bit to bit) in the three steps.
        call check_approx(passed, trim(label), stress_fea, stress_mem, tol=0D0)
    end do

    ! Measured: 0 (bit to bit), the memory cycle and the wrapper end in the same state.
    call check_approx(passed, "FEA cycle: history tensors of state_n", hist%state_n%h(:, 1:2), &
                      hist_mem%state_n%h(:, 1:2), tol=0D0)
    call check_approx(passed, "FEA cycle: elastic strain of state_n", hist%state_n%strain_e, &
                      hist_mem%state_n%strain_e, tol=0D0)
    call check_approx(passed, "FEA cycle: strain_e of the last pack", hsv(1:6), &
                      hist_mem%state_n%strain_e%vals, tol=0D0)
    call check_approx(passed, "FEA cycle: branch 1 of the last pack", hsv(7:12), &
                      hist_mem%state_n%h(:, 1), tol=0D0)
    call check_approx(passed, "FEA cycle: branch 2 of the last pack", hsv(13:18), &
                      hist_mem%state_n%h(:, 2), tol=0D0)

    ! A step that fails before any local solve writes the array back unchanged: the candidate
    ! state starts equal to the converged one, so packing right after unpacking is the identity.
    ! Measured: 0 (bit to bit).
    call model%init_history(hist)
    call hist%unpack_from_fea(hsv)
    hsv_out = -1D0                   ! every position has to be written by the pack
    call hist%pack_to_fea(hsv_out)
    call check_approx(passed, "FEA cycle: pack right after unpack returns the same array", &
                      hsv_out, hsv, tol=0D0)

    ! A new material point reuses the same state object: init_history must leave it as a brand
    ! new one (zero history), so the frozen stress equals the one of an untouched state.
    ! Measured: 0 (bit to bit).
    coefs = model%coefficients(0.3D0)
    call model%init_history(hist)          ! still holds the history of the last step
    call model%init_history(hist_new)
    elas = model%freeze(coefs, hist)
    elas_new = model%freeze(coefs, hist_new)
    call check_approx(passed, "FEA cycle: init_history clears the history", &
                      elas%stress(strain), elas_new%stress(strain), tol=0D0)
end subroutine test_maxwell_fea_wrapper_cycle

subroutine test_maxwell_frozen_tangents_2d_and_dev(passed)
    !! bulk(), shear(), the plane-strain tangent and the deviatoric tangent of
    !! Elasticity_frozen, against Elasticity_linear with the same pair (K_alg, G_alg), in a state
    !! with normals and shears coupled and a non-zero history stress.
    !! Hard case: the frozen tangent is not the plain elastic one, and the 2-D and deviatoric
    !! forms are not covered by the 3-D tangent check.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Elasticity_linear) :: lin
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain
    type(ten_2D2Osym) :: strain_2d
    real(real64), parameter :: K = 2000D0, G_INF = 100D0
    real(real64) :: e_alg, nu_alg, sum_kg

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=G_INF, shear=[200D0, 300D0], tau=[0.1D0, 1D0])
    call model%init_history(hist)

    ! A previous step leaves a non-zero history stress, so the frozen moduli are (K, G_alg)
    call strain%init(xx=1D-3, yy=-4D-4, zz=0D0, xy=6D-4, yz=0D0, xz=0D0)
    call model%update_dt(0.2D0, strain, hist)
    call hist%commit()
    call strain%init(xx=7D-4, yy=-1D-4, zz=2D-4, xy=2D-4, yz=-3D-4, xz=5D-4)
    coefs = model%coefficients(0.3D0)
    elas = model%freeze(coefs, hist)

    ! Measured: 0 (the frozen object stores the coefficients it was given).
    call check_approx(passed, "frozen bulk() = K_alg", elas%bulk(), coefs%k_alg, tol=0D0)
    call check_approx(passed, "frozen shear() = G_alg", elas%shear(), coefs%g_alg, tol=0D0)

    ! Reference: the linear elasticity of the library with the same pair of moduli
    sum_kg = 3D0*coefs%k_alg + coefs%g_alg
    e_alg  = 9D0*coefs%k_alg*coefs%g_alg/sum_kg
    nu_alg = (3D0*coefs%k_alg - 2D0*coefs%g_alg)/(2D0*sum_kg)
    call lin%set_parameters(young=e_alg, poisson=nu_alg)

    ! The reference is checked before being used: C : strain returns its own stress
    ! Measured: 0 and 1.3e-16 or less.
    call strain_2d%init(xx=strain%xx(), yy=strain%yy(), zz=strain%zz(), xy=strain%xy())
    call check_approx(passed, "reference: 2-D tangent .ddot. strain = 2-D stress", &
                      lin%dstress_dstrain_2D(strain_2d) .ddot. strain_2d, &
                      lin%stress_2D(strain_2d), tol=1D-14)
    call check_approx(passed, "reference: deviatoric tangent .ddot. strain = deviatoric stress", &
                      lin%dstress_dstrain_dev_3D(strain) .ddot. strain, &
                      lin%stress_dev_3D(strain), tol=1D-14)

    ! Measured: 5.3e-16 and 2.0e-16 relative (K_alg, G_alg -> E_alg, nu_alg round trip).
    call check_approx(passed, "frozen 2-D tangent = linear elastic 2-D tangent", &
                      elas%dstress_dstrain_2D(strain_2d), lin%dstress_dstrain_2D(strain_2d), &
                      tol=1D-14)
    call check_approx(passed, "frozen deviatoric tangent = linear elastic deviatoric tangent", &
                      elas%dstress_dstrain_dev_3D(strain), lin%dstress_dstrain_dev_3D(strain), &
                      tol=1D-14)
end subroutine test_maxwell_frozen_tangents_2d_and_dev

subroutine test_maxwell_gamma_series_small_x(passed)
    !! One step of constant deviatoric strain rate with x = dt/tau = 1e-6, where the model uses
    !! the three-term series of gamma, against the closed form
    !! sigma(t) = 2 dev(D) G_1 tau (1 - exp(-x)) of a single Maxwell element (G_inf = 0).
    !! A second probe at x = 9.9e-4, just below the series threshold, is the worst case for a
    !! wider threshold: the truncated series loses x^3/24 there.
    !! Hard case: no case of this file reaches the series branch (the smallest x of
    !! test_constant_rate is 0.03), so a wrong series is invisible until x is small; the sign of
    !! the linear term and the position of the threshold show up at x ~ 1e-6 and x ~ 1e-3.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: rate, strain, stress, expected
    real(real64), parameter :: K = 2000D0, G_1 = 100D0, TAU = 1D4
    real(real64), parameter :: DTS(2) = [1D-2, 9.9D0]     ! x = 1e-6 and x = 9.9e-4
    real(real64) :: x
    integer :: i
    character(len=40) :: label

    passed = .true.
    ! Single Maxwell element and a purely deviatoric rate: the whole stress comes from the
    ! branch, so the relative error of gamma is the relative error of the stress.
    call model%set_parameters(bulk=K, shear_inf=0D0, shear=[G_1], tau=[TAU])
    call rate%init(xx=2D-4, yy=-1D-4, zz=-1D-4, xy=1D-4, yz=3D-5, xz=-2D-5)

    do i = 1, size(DTS)
        x = DTS(i)/TAU
        strain = DTS(i)*rate
        call model%init_history(hist)       ! virgin state: the step starts from zero strain
        coefs = model%coefficients(DTS(i))
        elas = model%freeze(coefs, hist)
        stress = elas%stress(strain)
        ! 1 - exp(-x) = 2 exp(-x/2) sinh(x/2): no cancellation for small x
        expected = 2D0*G_1*TAU*(2D0*exp(-x/2D0)*sinh(x/2D0))*(.dev. rate)
        write(label, "(A, ES9.2)") "gamma series: stress at x =", x
        ! Covers the series branch (x = 1e-6) and the exact branch next to the threshold
        ! (x = 9.9e-4). Measured: 1.1e-19 and 4.3e-15 (2.7e-16 and 1.1e-14 relative).
        call check_approx(passed, trim(label), stress, expected, tol=1D-12)
    end do
end subroutine test_maxwell_gamma_series_small_x

subroutine test_maxwell_init_history_strain_e(passed)
    !! init_history with the initial elastic strain: the state stores it, the frozen stress at
    !! that strain is the fully relaxed one (the branches carry no stress, so only G_inf
    !! responds) and a step that keeps the strain does not create branch history.
    !! Hard case: with the initial strain dropped the same sequence looks like a step from zero
    !! strain and returns the instantaneous response instead of the relaxed one.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain, expected, vol
    real(real64), parameter :: K = 2000D0, G_INF = 100D0
    real(real64) :: tr

    passed = .true.
    call model%set_parameters(bulk=K, shear_inf=G_INF, shear=[200D0, 300D0], tau=[0.1D0, 1D0])
    call strain%init(xx=1D-3, yy=-4D-4, zz=2D-4, xy=6D-4, yz=-2D-4, xz=1D-4)
    call model%init_history(hist, strain_e=strain)   ! elastic strain, no branch stress

    ! Measured: 0 (bit to bit).
    call check_approx(passed, "init_history: initial elastic strain stored", &
                      hist%state_n%strain_e, strain, tol=0D0)

    tr = strain%xx() + strain%yy() + strain%zz()
    call vol%init(xx=K*tr, yy=K*tr, zz=K*tr, xy=0D0, yz=0D0, xz=0D0)
    expected = vol + 2D0*G_INF*(.dev. strain)   ! relaxed state: only K and G_inf respond

    coefs = model%coefficients(0.5D0)
    elas = model%freeze(coefs, hist)
    ! Measured: 2.2e-16 (1.3e-16 relative); G_alg and sum G_k gamma_k cancel to G_inf.
    call check_approx(passed, "init_history: stress at the initial strain is the relaxed one", &
                      elas%stress(strain), expected, tol=1D-12)

    ! The increment of the next step is measured from the initial elastic strain
    call model%update(coefs, strain, hist)
    call hist%commit()
    elas = model%freeze(coefs, hist)
    call check_approx(passed, "init_history: a step at constant strain adds no branch stress", &
                      elas%stress(strain), expected, tol=1D-12)
end subroutine test_maxwell_init_history_strain_e

subroutine test_maxwell_normalized_sum_g_one(passed)
    !! Normalized input with sum(g) = 1, that is, a Maxwell fluid: G_inf = 0, G_0 from E and nu,
    !! and the deviatoric stress relaxes to zero, leaving only K tr(eps) I.
    !! Hard case: sum([0.2, 0.4, 0.3, 0.1]) = 1.0000000000000002 in double precision, so the
    !! rounding of the sum is the whole problem, and a negative G_inf would show up here.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils
    use muscle_visco_state, only : Viscoelastic_history
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_visco_base, only : Visco_coefficients
    use muscle_visco_maxwell, only : Generalized_maxwell
    implicit none
    logical, intent(out) :: passed

    type(Generalized_maxwell) :: model
    type(Viscoelastic_history) :: hist
    type(Elasticity_frozen) :: elas
    type(Visco_coefficients) :: coefs
    type(ten_3D2Osym) :: strain, stress, vol
    real(real64), parameter :: E = 1000D0, NU = 0.3D0
    real(real64), parameter :: K = E/(3D0*(1D0 - 2D0*NU)), G_0 = 0.5D0*E/(1D0 + NU)
    real(real64) :: tr

    passed = .true.
    call model%set_parameters_normalized(young=E, poisson=NU, g=[0.2D0, 0.4D0, 0.3D0, 0.1D0], &
                                         tau=[0.1D0, 1D0, 10D0, 100D0])

    ! Measured: 0 and 0 (bit to bit).
    call check_approx(passed, "sum(g) = 1: G_0", model%shear_instantaneous(), G_0, tol=0D0)
    call check_approx(passed, "sum(g) = 1: G(1e12) = 0", model%relaxation_shear(1D12), 0D0, tol=0D0)

    ! A strain applied instantaneously and then held: only the bulk part is left
    call model%init_history(hist)
    call strain%init(xx=1D-3, yy=-2D-4, zz=3D-4, xy=5D-4, yz=-1D-4, xz=2D-4)
    call model%update_dt(1D-9, strain, hist)
    call hist%commit()
    coefs = model%coefficients(1D12)          ! much longer than every tau
    elas = model%freeze(coefs, hist)
    stress = elas%stress(strain)
    tr = strain%xx() + strain%yy() + strain%zz()
    call vol%init(xx=K*tr, yy=K*tr, zz=K*tr, xy=0D0, yz=0D0, xz=0D0)
    ! Measured: 1.1e-16 (1.2e-16 relative).
    call check_approx(passed, "sum(g) = 1: stress after a long hold = K tr(eps) I", &
                      stress, vol, tol=1D-9)
end subroutine test_maxwell_normalized_sum_g_one
