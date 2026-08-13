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