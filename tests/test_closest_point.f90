program test_closest_point
    implicit none
    logical :: passed

    call test_closest_point_vonmises_uniaxial_tensile(passed)
    if (.not. passed) STOP 1

    call test_closest_point_vonmises_zero_strain(passed)
    if (.not. passed) STOP 2

    call test_closest_point_vonmises_elastic_strain(passed)
    if (.not. passed) STOP 3
 
    print*, "Passed!", passed
end program test_closest_point

subroutine test_closest_point_vonmises_uniaxial_tensile(passed)
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening, only : Swift_hardening
    use mod_vonMises, only : VonMises
    use mod_elasticity_linear, only : Elasticity_linear
    ! use mod_closest_point, only : closest_point2
    use mod_closest_point
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed
    type(Closest_point_data) :: data
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

    call data%init(strain_pf=strain_pf, strain_p=strain_p)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress,        & 
                  strain_pf=strain_pf,  &
                  strain_p=strain_p,    &
                  status=status,        &
                  iters=iters           &
                  )

    passed = stress .isequal. expected_stress
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

    passed = strain_p .isequal. expected_strain_plastic
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


    call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
    call solver%tangent(strain=strain, data=data, tangent=tangent)

    passed = tangent .isequal. numerical_tangent
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return


    return
end subroutine


subroutine test_closest_point_vonmises_zero_strain(passed)
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening, only : Swift_hardening
    use mod_vonMises, only : VonMises
    use mod_elasticity_linear, only : Elasticity_linear
    use mod_closest_point, only : Closest_point, Closest_point_data
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(Closest_point_data) :: data
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

    call data%init(strain_pf=strain_pf, strain_p=strain_p)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress,        & 
                  strain_pf=strain_pf,  &
                  strain_p=strain_p     &
                  )

    passed = stress .isequal. expected_stress
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

    passed = strain_p .isequal. expected_strain_plastic
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return


    call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
    call solver%tangent(strain=strain, data=data, tangent=tangent)

    passed = tangent .isequal. numerical_tangent
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return
    

    return
end subroutine



subroutine test_closest_point_vonmises_elastic_strain(passed)
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening, only : Swift_hardening
    use mod_vonMises, only : VonMises
    use mod_elasticity_linear, only : Elasticity_linear
    use mod_closest_point, only : Closest_point, Closest_point_data
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(Closest_point_data) :: data
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

    call data%init(strain_pf=strain_pf, strain_p=strain_p)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress,        & 
                  strain_pf=strain_pf,  &
                  strain_p=strain_p     &
                  )

    passed = stress .isequal. expected_stress
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

    passed = strain_p .isequal. expected_strain_plastic
    if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return


    call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
    call solver%tangent(strain=strain, data=data, tangent=tangent)

    passed = tangent .isequal. numerical_tangent
    if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return

    return
end subroutine