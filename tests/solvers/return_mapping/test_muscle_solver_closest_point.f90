program test_muscle_solver_closest_point
    implicit none
    logical :: passed

    call test_closest_point_vonmises_uniaxial_tensile(passed)
    if (.not. passed) STOP 1

    call test_closest_point_vonmises_zero_strain(passed)
    if (.not. passed) STOP 2

    call test_closest_point_vonmises_elastic_strain(passed)
    if (.not. passed) STOP 3

    ! call test_closest_point_druckerPrager_uniaxial_tensile(passed)
    ! if (.not. passed) STOP 4
 
    call test_closest_point_packed_shear(passed)
    if (.not. passed) STOP 5

    call test_closest_point_keeps_state_n(passed)
    if (.not. passed) STOP 6

    call test_closest_point_hill48_tension_shear(passed)
    if (.not. passed) STOP 7

    print*, "Passed!", passed
end program test_muscle_solver_closest_point

subroutine test_closest_point_vonmises_uniaxial_tensile(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS=1e-5
    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    integer :: status, iters
    type(ten_3D4O2sym) :: tangent, numerical_tangent

    passed = .False.
    strain_pf = 0D0
    error = .False.
    

    call expected_stress%init(xx=90D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    expected_strain_effective = 0.3485784401D0
    call expected_strain_plastic%init(xx=expected_strain_effective,        &
                                      yy=-0.5D0*expected_strain_effective, &
                                      zz=-0.5D0*expected_strain_effective, &
                                      xy=0D0, yz=0D0, xz=0D0)
                                      
    call strain%init(xx=0.43857844D0, yy=-0.20128922D0, zz=-0.20128922D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call history%init(strain_p=strain_p, strain_pf=strain_pf)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    
    stress    = history%state_np1%stress
    strain_p  = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf

    passed = stress%is_approx(expected_stress, tol=EPS)
    if (.not. passed) print*, "Error: Stress is no equal", new_line('A'),          &
                              "Expected:", expected_stress, new_line('A'),  &
                              "Actual Value:", stress, new_line('A'),       &
                              "Difference", stress - expected_stress
    if (.not. passed) return
    
    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p%is_approx(expected_strain_plastic, tol=EPS)
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return
    
    passed = iters .le. 6
    if(.not. passed) print*, "Iterations are greater than expected", new_line('A'), &
                             "Iters:", iters, " Maximum:", 6 , new_line('A'), &
                             "Status code: ", status
    if(.not. passed) return


    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)

    passed = tangent%is_approx(numerical_tangent, tol=1.0D-7)
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return


    return
end subroutine


subroutine test_closest_point_vonmises_zero_strain(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    integer :: status, iters


    passed = .False.
    strain_pf = 0D0
    error = .False.

    expected_strain_effective = 0.0D0
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, &
                                      xy=0D0, yz=0D0, xz=0D0)
                                      
    call strain%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call history%init(strain_p=strain_p, strain_pf=strain_pf)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    
    stress    = history%state_np1%stress
    strain_p  = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf

    passed = stress .approx. expected_stress
    if (.not. passed) print*, "Stress is no equal", new_line('A'),          &
                              "Expected:", expected_stress, new_line('A'),  &
                              "Actual Value:", stress, new_line('A'),       &
                              "Difference", stress - expected_stress
    if (.not. passed) return
    
    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p .approx. expected_strain_plastic
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return


    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)

    passed = tangent .approx. numerical_tangent
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return
    

    return
end subroutine



subroutine test_closest_point_vonmises_elastic_strain(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    integer :: status, iters


    passed = .False.
    strain_pf = 0D0
    error = .False.

    expected_strain_effective = 0.0D0
    call expected_stress%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, &
                                      xy=0D0, yz=0D0, xz=0D0)
                                      
    call strain%init(xx=0.001D0, yy=-0.0003D0, zz=-0.0003D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call history%init(strain_p=strain_p, strain_pf=strain_pf)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)
    
    stress    = history%state_np1%stress
    strain_p  = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf

    passed = stress .approx. expected_stress
    if (.not. passed) print*, "Stress is no equal", new_line('A'),          &
                              "Expected:", expected_stress, new_line('A'),  &
                              "Actual Value:", stress, new_line('A'),       &
                              "Difference", stress - expected_stress
    if (.not. passed) return
    
    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p .approx. expected_strain_plastic
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return


    call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
    call solver%tangent(strain=strain, history=history, tangent=tangent)

    passed = tangent .approx. numerical_tangent
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return

    return
end subroutine

subroutine test_closest_point_druckerPrager_uniaxial_tensile(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_bilinear
    use muscle_yield_druckerprager
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point
    implicit none

    real(real64), parameter :: EPS=1e-5
    logical, intent(out) :: passed
    type(Plastic_material_history) :: history
    type(Closest_point) :: solver
    type(DruckerPrager) :: dp
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Bilinear_hardening) :: bilinear
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    real(real64) :: beta, K_ratio

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    integer :: status, iters
    real(real64) :: strain_xx, strain_yy

    passed = .False.
    strain_pf = 0D0

    call elas%set_parameters(young=210D3, poisson=0.25D0)
    bilinear = Bilinear_hardening(y0=300D0, k=1000D0)
    beta = 16.0D0
    K_ratio = 0.85D0

    ! ---------------------------------------------------
    ! --------  Tensile type  ---------------------------
    ! ---------------------------------------------------

    strain_xx = 0.09531D0
    strain_yy = -0.03666D0
    call dp%init(beta_deg=beta, K=K_ratio, hardening_mode=DP_HARDENING_TENSION)


    call expected_stress%init(xx=393.4D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0.09344D0,  &
                                      yy=-0.03619D0, &
                                      zz=-0.03619D0, &
                                      xy=0D0, yz=0D0, xz=0D0)
    expected_strain_effective = 0.09344D0

    call strain%init(xx=strain_xx, yy=strain_yy, zz=strain_yy, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)

    call history%init(strain_p=strain_p, strain_pf=strain_pf)
    call solver%init(elasticity=elas, hardening=bilinear, yield=dp)
    call solver%solve(strain=strain, history=history, status=status, iters=iters)

    stress    = history%state_np1%stress
    strain_p  = history%state_np1%strain_p
    strain_pf = history%state_np1%strain_pf

    ! passed = stress%is_approx(expected_stress, tol=EPS)
    ! if (.not. passed) print*, "Error: Stress is no equal", iters, new_line('A'),          &
    !                           "Expected:", expected_stress, new_line('A'),  &
    !                           "Actual Value:", stress, new_line('A'),       &
    !                           "Difference", stress - expected_stress
    ! if (.not. passed) return

    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p%is_approx(expected_strain_plastic, tol=EPS)
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return

    passed = iters .le. 6
    if(.not. passed) print*, "Iterations are greater than expected", new_line('A'), &
                             "Iters:", iters, " Maximum:", 6 , new_line('A'), &
                             "Status code: ", status
    if(.not. passed) return

    return
end subroutine test_closest_point_druckerPrager_uniaxial_tensile

! FE boundary: the FE program reads the result through pack_to_fea and commits only after global
! equilibrium. Both cases use von Mises with Swift hardening and a strain with all six components
! nonzero, written as a multiple of the yield strain eps_y along a fixed direction.

subroutine test_closest_point_packed_shear(passed)
    ! The FE program reads the plastic strain from slots 7-12. Its shears must be the tensorial
    ! components (eps_xy, not gamma_xy = 2*eps_xy): with them the packed stress is the elastic
    ! law of the total minus the plastic strain.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_hard_swift
    use muscle_yield_vonmises
    use muscle_elasticity_linear
    use muscle_plastic_history
    use muscle_solver_closest_point
    implicit none
    logical, intent(out) :: passed

    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: hardening
    type(VonMises) :: vm
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain, packed_stress, packed_strain_p
    real(real64) :: hsv(13), eps_y
    integer :: status
    real(real64), parameter :: TOL = 1.0D-10
    ! Newton tolerance of the solver
    real(real64), parameter :: TOL_NW = 1.0D-5

    call elas%set_parameters(young=70000.0D0, poisson=0.33D0)
    hardening = Swift_hardening(k=500.0D0, n=0.2D0, e0=0.01D0)
    call solver%init(elasticity=elas, hardening=hardening, yield=vm)
    eps_y = hardening%stress(0.0D0)/70000.0D0

    ! Five times the yield strain: plastic, with large plastic shears
    call strain%init(xx=4.0D0*eps_y, yy=-1.5D0*eps_y, zz=-2.5D0*eps_y, &
                     xy=2.0D0*eps_y, yz=-1.0D0*eps_y, xz=1.5D0*eps_y)
    call solver%solve(strain=strain, history=history, status=status)
    if (status /= STATUS_CONVERGED) then
        print*, "FAIL: Packed shear case did not converge, status:", status
        passed = .false.
        return
    end if

    hsv = 0.0D0
    call history%pack_to_fea(hsv)
    call packed_stress%init(hsv(1:6))
    call packed_strain_p%init(hsv(7:12))

    passed = packed_stress%is_approx(elas%stress(strain - packed_strain_p), tol=TOL)
    if (.not. passed) then
        print*, "FAIL: Packed stress is not the elastic law of the packed (tensorial) plastic strain"
        return
    end if

    ! One step from the virgin state: eps_p = dgamma*n with the von Mises normal n of unit
    ! equivalent norm, so slot 13 is sqrt(2/3 eps_p:eps_p)
    passed = abs(hsv(13) - sqrt(2.0D0/3.0D0*(packed_strain_p .ddot. packed_strain_p))) < TOL_NW
    if (.not. passed) then
        print*, "FAIL: Slot 13 is not the equivalent plastic strain of slots 7-12"
        return
    end if
end subroutine test_closest_point_packed_shear


subroutine test_closest_point_keeps_state_n(passed)
    ! The FE program commits only after global equilibrium. Before that, neither a failed attempt
    ! (rolled back and retried) nor a converged one may write the committed state_n.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_hard_swift
    use muscle_yield_vonmises
    use muscle_elasticity_linear
    use muscle_plastic_history
    use muscle_solver_closest_point
    implicit none
    logical, intent(out) :: passed

    type(Elasticity_linear) :: elas
    type(Swift_hardening) :: hardening
    type(VonMises) :: vm
    type(Closest_point) :: solver, one_iteration
    type(Plastic_material_history) :: history
    type(Plastic_state) :: state_before
    type(ten_3D2Osym) :: strain_1, strain_2
    real(real64) :: eps_y
    integer :: status
    real(real64), parameter :: TOL = 1.0D-12

    call elas%set_parameters(young=70000.0D0, poisson=0.33D0)
    hardening = Swift_hardening(k=500.0D0, n=0.2D0, e0=0.01D0)
    call solver%init(elasticity=elas, hardening=hardening, yield=vm)
    ! One Newton iteration is not enough for step 2: forces a failed attempt
    call one_iteration%init(elasticity=elas, hardening=hardening, yield=vm, iter_nw=1)
    eps_y = hardening%stress(0.0D0)/70000.0D0

    ! Step 1, converged and committed, so state_n is a plastic state and not the virgin one
    call strain_1%init(xx=1.6D0*eps_y, yy=-0.6D0*eps_y, zz=-1.0D0*eps_y, &
                       xy=0.8D0*eps_y, yz=-0.4D0*eps_y, xz=0.6D0*eps_y)
    call solver%solve(strain=strain_1, history=history, status=status)
    if (status /= STATUS_CONVERGED) then
        print*, "FAIL: Step 1 did not converge, status:", status
        passed = .false.
        return
    end if
    call history%commit()
    state_before = history%state_n

    ! Step 2: failed attempt with a single Newton iteration, rollback, converged retry
    call strain_2%init(xx=3.2D0*eps_y, yy=-1.2D0*eps_y, zz=-2.0D0*eps_y, &
                       xy=1.6D0*eps_y, yz=-0.8D0*eps_y, xz=1.2D0*eps_y)
    call one_iteration%solve(strain=strain_2, history=history, status=status)
    if (status /= STATUS_NONCONVERGED) then
        print*, "FAIL: Forced attempt was expected not to converge, status:", status
        passed = .false.
        return
    end if
    call history%rollback()
    call solver%solve(strain=strain_2, history=history, status=status)
    if (status /= STATUS_CONVERGED) then
        print*, "FAIL: Retry of step 2 did not converge, status:", status
        passed = .false.
        return
    end if

    passed = history%state_n%stress%is_approx(state_before%stress, tol=TOL) .and. &
             history%state_n%strain_p%is_approx(state_before%strain_p, tol=TOL) .and. &
             abs(history%state_n%strain_pf - state_before%strain_pf) <= TOL*state_before%strain_pf
    if (.not. passed) then
        print*, "FAIL: An attempt before commit changed state_n"
        return
    end if
end subroutine test_closest_point_keeps_state_n

subroutine test_closest_point_hill48_tension_shear(passed)
    ! Hill48 with the AA2090-T3 r-values of test_Hill48_gradient (Yoon and Barlat,
    ! 2006), one Closest_point step from zero plastic strain.
    ! Expected stress and plastic strain: associated flow of Hill (1948, Eq. 5)
    ! at equivalent plastic strain 0.05, Swift hardening, evaluated with mpmath
    ! (50 digits). Not taken from this solver.
    ! Elasticity and Swift match test_closest_point_vonmises_uniaxial_tensile.
    ! RD exercises the normal block. 45 deg also turns on N and the xy shear.
    ! Pure shear in xy, yz and xz checks N, L and M one at a time.
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point, only : Closest_point, STATUS_CONVERGED
    implicit none

    real(real64), parameter :: EPS = 1.0D-5
    real(real64), parameter :: EPS_PF = 1.0D-7
    real(real64), parameter :: EPS_R = 1.0D-6
    real(real64), parameter :: EPS_TAN = 1.0D-7
    logical, intent(out) :: passed
    type(Hill48) :: h48
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: r0, r45, r90

    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    ! The r-values do not give L and M: they take distinct values, also distinct from N,
    ! so that each pure shear below is resisted by its own coefficient.
    h48 = Hill48(f=r0/(r90*(1.0D0 + r0)), g=1.0D0/(1.0D0 + r0), h=r0/(1.0D0 + r0), &
                 l=1.25D0, m=1.75D0, &
                 n=(r0 + r90)*(1.0D0 + 2.0D0*r45)/(2.0D0*r90*(1.0D0 + r0)))
    call elas%set_parameters(young=1000.0D0, poisson=0.3D0)
    sw = Swift_hardening(k=100.0D0, n=0.1D0, e0=1.0D-4)

    ! 45 deg first: N and the tensorial xy shear are active. Width strain uses eps_p,xy.
    call check_direction(passed, h48, sw, elas, 45.0D0, r45, &
        3.1455217937819118D+1, 3.1455217937819118D+1, 0.0D0, 3.1455217937819118D+1, &
        3.9531449486028268D-2, 2.7368871465329165D-2, -4.1736146601102139D-2, &
        8.8375980913959970D-2, &
        1.7512796929554885D-2, 5.3502189088557826D-3, -2.2863015838410668D-2, &
        4.7484197594795116D-2)
    if (.not. passed) return

    ! RD: normal anisotropy only. r0 = H/G is eps_p,yy / eps_p,zz.
    call check_direction(passed, h48, sw, elas, 0.0D0, r0, &
        7.4128254276130228D+1, 0.0D0, 0.0D0, 0.0D0, &
        1.2412825427613023D-1, -3.0967324817713192D-2, -6.3509627747964946D-2, 0.0D0, &
        5.0000000000000000D-2, -8.7288485348741230D-3, -4.1271151465125877D-2, 0.0D0)
    if (.not. passed) return

    ! Pure shear (Voigt slot 4, 5 or 6): tau = sigma_y/sqrt(2 C) and eps_p = 0.05 sqrt(C/2),
    ! with C = N, L or M. A swapped coefficient or an engineering factor 2 changes tau.
    call check_shear(passed, h48, sw, elas, 4, &
        3.5037546580065136D+1, 9.8440823203099628D-2, 5.2892012649014951D-2)
    if (.not. passed) return

    call check_shear(passed, h48, sw, elas, 5, &
        4.6882824496937552D+1, 1.0047614259812356D-1, 3.9528470752104742D-2)
    if (.not. passed) return

    call check_shear(passed, h48, sw, elas, 6, &
        3.9623218597277097D+1, 9.8280901511134493D-2, 4.6770717334674267D-2)

contains

    ! Resolve and check one prescribed uniaxial Hill48 state (Hill, 1948, Eq. 5).
    ! Called by test_closest_point_hill48_tension_shear.
    subroutine check_direction(ok, h48, sw, elas, theta_deg, r_in, &
                               sxx, syy, szz, sxy, exx, eyy, ezz, exy, &
                               pxx, pyy, pzz, pxy)
        implicit none
        logical, intent(out) :: ok
        type(Hill48), intent(in) :: h48
        type(Swift_hardening), intent(in) :: sw
        type(Elasticity_linear), intent(in) :: elas
        real(real64), intent(in) :: theta_deg, r_in
        real(real64), intent(in) :: sxx, syy, szz, sxy
        real(real64), intent(in) :: exx, eyy, ezz, exy
        real(real64), intent(in) :: pxx, pyy, pzz, pxy
        type(ten_3D2Osym) :: strain, strain_p
        type(ten_3D2Osym) :: expected_stress, expected_plastic
        type(ten_3D2Osym) :: width_projection ! t (x) t, transverse to the specimen axis
        real(real64) :: r_theta, th, sn, cs, width

        ! Uniaxial tension along the specimen axis (xy = sigma/2 at 45 deg).
        ! At 45 deg the xy component of the flow is the one that carries N.
        call expected_stress%init(xx=sxx, yy=syy, zz=szz, xy=sxy, yz=0.0D0, xz=0.0D0)
        call expected_plastic%init(xx=pxx, yy=pyy, zz=pzz, xy=pxy, yz=0.0D0, xz=0.0D0)
        call strain%init(xx=exx, yy=eyy, zz=ezz, xy=exy, yz=0.0D0, xz=0.0D0)
        call check_step(ok, h48, sw, elas, theta_deg, strain, expected_stress, &
                        expected_plastic, strain_p)
        if (.not. ok) return

        ! Lankford r = width/thickness must equal the r used to build F, G, H, N.
        th = theta_deg*acos(-1.0D0)/180.0D0
        sn = sin(th)
        cs = cos(th)
        call width_projection%init(xx=sn*sn, yy=cs*cs, zz=0.0D0, &
                                   xy=-sn*cs, yz=0.0D0, xz=0.0D0)
        width = width_projection .ddot. strain_p
        r_theta = width/strain_p%zz()
        ok = abs(r_theta - r_in) < EPS_R
        if (.not. ok) print *, "Hill48 r", theta_deg, r_theta, r_in
    end subroutine check_direction

    ! Resolve and check one prescribed pure shear Hill48 state (Hill, 1948, Eq. 5): only the
    ! tensorial shear in Voigt slot i is nonzero, in the strain, the stress and the flow.
    ! Called by test_closest_point_hill48_tension_shear.
    subroutine check_shear(ok, h48, sw, elas, i, tau, e, p)
        implicit none
        logical, intent(out) :: ok
        type(Hill48), intent(in) :: h48
        type(Swift_hardening), intent(in) :: sw
        type(Elasticity_linear), intent(in) :: elas
        integer, intent(in) :: i          ! Voigt slot: 4 = xy, 5 = yz, 6 = xz
        real(real64), intent(in) :: tau   ! Expected shear stress
        real(real64), intent(in) :: e     ! Prescribed tensorial shear strain
        real(real64), intent(in) :: p     ! Expected tensorial plastic shear strain
        type(ten_3D2Osym) :: strain, strain_p, expected_stress, expected_plastic
        real(real64) :: slot(6)           ! Unit vector on slot i

        slot = 0.0D0
        slot(i) = 1.0D0
        call strain%init(e*slot)
        call expected_stress%init(tau*slot)
        call expected_plastic%init(p*slot)
        call check_step(ok, h48, sw, elas, real(i, real64), strain, expected_stress, &
                        expected_plastic, strain_p)
    end subroutine check_shear

    ! One Closest_point step from zero plastic strain: status, stress, strain_pf, strain_p,
    ! iterations and consistent tangent. tag only labels the messages (angle or Voigt slot).
    ! Called by check_direction and check_shear.
    subroutine check_step(ok, h48, sw, elas, tag, strain, expected_stress, expected_plastic, &
                          strain_p)
        implicit none
        logical, intent(out) :: ok
        type(Hill48), intent(in) :: h48
        type(Swift_hardening), intent(in) :: sw
        type(Elasticity_linear), intent(in) :: elas
        real(real64), intent(in) :: tag
        type(ten_3D2Osym), intent(in) :: strain, expected_stress, expected_plastic
        type(ten_3D2Osym), intent(out) :: strain_p
        type(Plastic_material_history) :: history
        type(Closest_point) :: solver
        type(ten_3D2Osym) :: strain_p0, stress
        type(ten_3D4O2sym) :: tangent, numerical_tangent
        real(real64) :: strain_pf
        integer :: status, iters

        ok = .false.
        call strain_p0%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, yz=0.0D0, xz=0.0D0)
        call history%init(strain_p=strain_p0, strain_pf=0.0D0)
        call solver%init(elasticity=elas, hardening=sw, yield=h48)
        call solver%solve(strain=strain, history=history, status=status, iters=iters)

        ! This prescribed step is plastic and must report convergence.
        ok = status == STATUS_CONVERGED
        if (.not. ok) print *, "Hill48 status", tag, status
        if (.not. ok) return

        stress = history%state_np1%stress
        strain_p = history%state_np1%strain_p
        strain_pf = history%state_np1%strain_pf

        ok = stress%is_approx(expected_stress, tol=EPS)
        if (.not. ok) print *, "Hill48 stress", tag, stress%vals
        if (.not. ok) return

        ! Plastic multiplier equal to the prescribed equivalent plastic strain.
        ok = abs(strain_pf - 0.05D0) < EPS_PF
        if (.not. ok) print *, "Hill48 strain_pf", tag, strain_pf
        if (.not. ok) return

        ! Associated flow.
        ok = strain_p%is_approx(expected_plastic, tol=EPS)
        if (.not. ok) print *, "Hill48 strain_p", tag, strain_p%vals
        if (.not. ok) return

        ! Same iteration bound as the von Mises uniaxial case.
        ok = iters .le. 6
        if (.not. ok) print *, "Hill48 iterations", tag, iters, "status", status
        if (.not. ok) return

        ! Consistent tangent against tangent_numerical.
        call solver%tangent_numerical(strain=strain, history=history, tangent=numerical_tangent)
        call solver%tangent(strain=strain, history=history, tangent=tangent)
        ok = tangent%is_approx(numerical_tangent, tol=EPS_TAN)
        if (.not. ok) print *, "Hill48 tangent", tag
    end subroutine check_step

end subroutine test_closest_point_hill48_tension_shear
