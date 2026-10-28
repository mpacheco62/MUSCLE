program test_muscle_solver_closest_point_time
    implicit none
    logical :: passed

    call test_closest_point_time_vonmises_uniaxial_tensile(passed)
    if (.not. passed) STOP 1

    call test_closest_point_time_non_converged(passed)
    if (.not. passed) STOP 2

    ! call test_closest_point_time_vonmises_zero_strain(passed)
    ! if (.not. passed) STOP 3

    ! call test_closest_point_time_vonmises_elastic_strain(passed)
    ! if (.not. passed) STOP 4
 
    call test_closest_point_time_reinit_other_types(passed)
    if (.not. passed) STOP 5
    call test_closest_point_time_vonmises_mrk_uniaxial_tensile(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
end program test_muscle_solver_closest_point_time

subroutine test_closest_point_time_vonmises_uniaxial_tensile(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_jc, only : JC_viscoplastic
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    ! use muscle_solver_closest_point, only : closest_point2
    use muscle_solver_closest_point_time
    implicit none

    real(real64), parameter :: EPS=1e-8
    real(real64), parameter :: TOL_TEN=1D-6  ! input strain has 8 digits: |stress_yy| ~ 1e-7*|stress_xx|
    logical, intent(out) :: passed
    type(Closest_point_data) :: data
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress

    type(JC_viscoplastic) :: jc
    type(Swift_hardening), target :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64) :: expected_strain_effective
    integer :: status, iters
    type(ten_3D4O2sym) :: tangent, numerical_tangent

    real(real64) :: xy

    ! passed = .False.
    ! strain_pf = 0D0
    ! error = .False.
    

    ! call expected_stress%init(xx=90D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    ! expected_strain_effective = 0.3485784401D0
    ! call expected_strain_plastic%init(xx=expected_strain_effective,        &
    !                                   yy=-0.5D0*expected_strain_effective, &
    !                                   zz=-0.5D0*expected_strain_effective, &
    !                                   xy=0D0, yz=0D0, xz=0D0)
                                      
    ! call strain%init(xx=0.43857844D0, yy=-0.20128922D0, zz=-0.20128922D0, xy=0D0, yz=0D0, xz=0D0)
    ! call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    ! call elas%set_parameters(young=1000D0, poisson=0.3D0)

    ! sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    ! jc = JC_viscoplastic(C=0D0, epdmax=1D0)
    ! jc%hard_law = sw  ! Error in gfortran-12 if not assigned here

    ! call data%init(strain_pf=strain_pf, strain_p=strain_p, dt=1.0D-3)
    ! call solver%init(elasticity=elas, hardening=jc, yield=vm)
    ! call solver%solve(strain=strain, data=data)
    ! call data%get(stress=stress,        & 
    !               strain_pf=strain_pf,  &
    !               strain_p=strain_p,    &
    !               status=status,        &
    !               iters=iters           &
    !               )

    ! passed = stress .approx. expected_stress
    ! if (.not. passed) print*, "Case 1 Stress is no equal", new_line('A'),          &
    !                           "Expected:", expected_stress, new_line('A'),  &
    !                           "Actual Value:", stress, new_line('A'),       &
    !                           "Difference", stress - expected_stress
    ! if (.not. passed) return
    
    ! passed = abs(strain_pf - expected_strain_effective) < EPS
    ! if (.not. passed) print*, "Case 2 Effective plastic Strain is not equal", new_line('A'), &
    !                           "Expected:", expected_strain_effective, new_line('A'),  &
    !                           "Actual Value:", strain_pf, new_line('A'),              &
    !                           "Difference", strain_pf - expected_strain_effective
    ! if (.not. passed) return

    ! passed = strain_p .approx. expected_strain_plastic
    ! if (.not. passed) print*, "Case 3 Plastic Strain is not equal", new_line('A'), &
    !                           "Expected:", expected_strain_plastic, new_line('A'), &
    !                           "Actual Value:", strain_p, new_line('A'),            & 
    !                           "Difference:", strain_p - expected_strain_plastic
    ! if (.not. passed) return
    
    ! passed = iters .le. 6
    ! if(.not. passed) print*, "Case 4 Iterations are greater than expected", new_line('A'), &
    !                          "Iters:", iters, " Maximum:", 6 , new_line('A'), &
    !                          "Status code: ", status
    ! if(.not. passed) return


    ! call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
    ! call solver%tangent(strain=strain, data=data, tangent=tangent)

    ! passed = tangent .approx. numerical_tangent
    ! if (.not. passed) print*, "Case 5 Tangent is no equal", new_line('A'),          &
    !                           "Analitical:", tangent, new_line('A'),  &
    !                           "Numerical:", numerical_tangent, new_line('A'),       &
    !                           "Difference", tangent - numerical_tangent
    ! if (.not. passed) return


    ! The test now implement rate-dependant constants
    

    data = Closest_point_data()

    passed = .False.
    strain_pf = 0D0
    error = .False.
    

    ! Uniaxial tension with JC: E*(exx - ep) = Swift(ep)*(1 + C*log(ep/dt/epdmax))
    call expected_stress%init(xx=100.6464834098D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    expected_strain_effective = 0.3379319566D0
    call expected_strain_plastic%init(xx=expected_strain_effective,        &
                                      yy=-0.5D0*expected_strain_effective, &
                                      zz=-0.5D0*expected_strain_effective, &
                                      xy=0D0, yz=0D0, xz=0D0)
    
    xy = -0.19915993D0
    call strain%init(xx=0.43857844D0, yy=xy, zz=xy, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    jc = JC_viscoplastic(C=1D-1, epdmax=1.0D-3)
    jc%hard_law = sw  ! Error in gfortran-12 if not assigned here

    ! data = Closest_point_data(strain_pf=strain_pf, strain_p=strain_p, dt=1.0D+2)
    call data%init(strain_pf=strain_pf, strain_p=strain_p, dt=1.0D+2)
    call solver%init(elasticity=elas, hardening=jc, yield=vm)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress,        & 
                  strain_pf=strain_pf,  &
                  strain_p=strain_p,    &
                  status=status,        &
                  iters=iters           &
                  )

    passed = stress%is_approx(expected_stress, tol=TOL_TEN)
    print*, "strain_pf", strain_pf, sw%stress(strain_pf), jc%flow_stress(strain_pf, strain_pf/1D+2)
    if (.not. passed) print*, "Case 6 Stress is no equal", new_line('A'),          &
                              "Expected:", expected_stress, new_line('A'),  &
                              "Actual Value:", stress, new_line('A'),       &
                              "Difference", stress - expected_stress
    if (.not. passed) return
    
    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) print*, "Case 7 Effective plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_effective, new_line('A'),  &
                              "Actual Value:", strain_pf, new_line('A'),              &
                              "Difference", strain_pf - expected_strain_effective
    if (.not. passed) return

    passed = strain_p%is_approx(expected_strain_plastic, tol=TOL_TEN)
    if (.not. passed) print*, "Case 8 Plastic Strain is not equal", new_line('A'), &
                              "Expected:", expected_strain_plastic, new_line('A'), &
                              "Actual Value:", strain_p, new_line('A'),            & 
                              "Difference:", strain_p - expected_strain_plastic
    if (.not. passed) return
    
    passed = iters .le. 6
    if(.not. passed) print*, "Case 9 Iterations are greater than expected", new_line('A'), &
                             "Iters:", iters, " Maximum:", 6 , new_line('A'), &
                             "Status code: ", status
    if(.not. passed) return


    call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
    call solver%tangent(strain=strain, data=data, tangent=tangent)

    passed = tangent%is_approx(numerical_tangent, tol=1.0D-7)
    if (.not. passed) print*, "Case 10 Tangent is no equal", new_line('A'),          &
                              "Analitical:", tangent, new_line('A'),  &
                              "Numerical:", numerical_tangent, new_line('A'),       &
                              "Difference", tangent - numerical_tangent
    if (.not. passed) return
end subroutine

subroutine test_closest_point_time_non_converged(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_jc, only : JC_viscoplastic
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    ! use muscle_solver_closest_point, only : closest_point2
    use muscle_solver_closest_point_time
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed
    type(Closest_point_data) :: data
    type(Closest_point) :: solver
    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p

    type(JC_viscoplastic) :: jc
    type(Swift_hardening), target :: sw
    type(Elasticity_linear) :: elas

    integer :: status, iters


    passed = .False.
    

    call strain%init(xx=0.43857844D0, yy=-0.19915993D0, zz=-0.19915993D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    jc = JC_viscoplastic(C=1D-1, epdmax=1.0D-3)
    jc%hard_law = sw  ! Error in gfortran-12 if not assigned here

    ! data = Closest_point_data(strain_pf=strain_pf, strain_p=strain_p, dt=1.0D+2)
    call data%init(strain_pf=0D0, strain_p=strain_p, dt=1.0D+2)
    ! Same state as the uniaxial case, which converges in 5 iterations
    call solver%init(elasticity=elas, hardening=jc, yield=vm, iter_nw=4)
    call solver%solve(strain=strain, data=data)
    call data%get(status=status,        &
                  iters=iters           &
                  )

    passed = (status .eq. STATUS_NONCONVERGED)
    if (.not. passed) print*, "Case 1 Converged when non converged must occur", new_line('A'),          &
                              "iters:", iters, new_line('A'),  &
                              "status:", status
    if (.not. passed) return
end subroutine

! subroutine test_closest_point_time_vonmises_zero_strain(passed)
!     use muscle_tensors
!     use, intrinsic :: iso_fortran_env, only : real64
!     use muscle_hard_swift, only : Swift_hardening
!     use muscle_yield_vonmises, only : VonMises
!     use muscle_elasticity_linear, only : Elasticity_linear
!     use muscle_solver_closest_point, only : Closest_point, Closest_point_data
!     implicit none

!     real(real64), parameter :: EPS=1e-8
!     logical, intent(out) :: passed

!     type(Closest_point_data) :: data
!     type(Closest_point) :: solver
!     type(VonMises) :: vm
!     type(ten_3D2Osym) :: strain, strain_p, stress
!     type(Swift_hardening) :: sw
!     type(Elasticity_linear) :: elas
!     real(real64) :: strain_pf
!     logical :: error

!     type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
!     real(real64) :: expected_strain_effective
!     type(ten_3D4O2sym) :: tangent, numerical_tangent


!     passed = .False.
!     strain_pf = 0D0
!     error = .False.

!     expected_strain_effective = 0.0D0
!     call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
!     call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, &
!                                       xy=0D0, yz=0D0, xz=0D0)
                                      
!     call strain%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=0D0, yz=0D0, xz=0D0)
!     call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
!     call elas%set_parameters(young=1000D0, poisson=0.3D0)
!     sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

!     call data%init(strain_pf=strain_pf, strain_p=strain_p)
!     call solver%init(elasticity=elas, hardening=sw, yield=vm)
!     call solver%solve(strain=strain, data=data)
!     call data%get(stress=stress,        & 
!                   strain_pf=strain_pf,  &
!                   strain_p=strain_p     &
!                   )

!     passed = stress .approx. expected_stress
!     if (.not. passed) print*, "Stress is no equal", new_line('A'),          &
!                               "Expected:", expected_stress, new_line('A'),  &
!                               "Actual Value:", stress, new_line('A'),       &
!                               "Difference", stress - expected_stress
!     if (.not. passed) return
    
!     passed = abs(strain_pf - expected_strain_effective) < EPS
!     if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
!                               "Expected:", expected_strain_effective, new_line('A'),  &
!                               "Actual Value:", strain_pf, new_line('A'),              &
!                               "Difference", strain_pf - expected_strain_effective
!     if (.not. passed) return

!     passed = strain_p .approx. expected_strain_plastic
!     if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
!                               "Expected:", expected_strain_plastic, new_line('A'), &
!                               "Actual Value:", strain_p, new_line('A'),            & 
!                               "Difference:", strain_p - expected_strain_plastic
!     if (.not. passed) return


!     call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
!     call solver%tangent(strain=strain, data=data, tangent=tangent)

!     passed = tangent .approx. numerical_tangent
!     if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
!                               "Analitical:", tangent, new_line('A'),  &
!                               "Numerical:", numerical_tangent, new_line('A'),       &
!                               "Difference", tangent - numerical_tangent
!     if (.not. passed) return
    

!     return
! end subroutine



! subroutine test_closest_point_time_vonmises_elastic_strain(passed)
!     use muscle_tensors
!     use, intrinsic :: iso_fortran_env, only : real64
!     use muscle_hard_swift, only : Swift_hardening
!     use muscle_yield_vonmises, only : VonMises
!     use muscle_elasticity_linear, only : Elasticity_linear
!     use muscle_solver_closest_point, only : Closest_point, Closest_point_data
!     implicit none

!     real(real64), parameter :: EPS=1e-8
!     logical, intent(out) :: passed

!     type(Closest_point_data) :: data
!     type(Closest_point) :: solver
!     type(VonMises) :: vm
!     type(ten_3D2Osym) :: strain, strain_p, stress
!     type(Swift_hardening) :: sw
!     type(Elasticity_linear) :: elas
!     real(real64) :: strain_pf
!     logical :: error

!     type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
!     real(real64) :: expected_strain_effective
!     type(ten_3D4O2sym) :: tangent, numerical_tangent


!     passed = .False.
!     strain_pf = 0D0
!     error = .False.

!     expected_strain_effective = 0.0D0
!     call expected_stress%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
!     call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, &
!                                       xy=0D0, yz=0D0, xz=0D0)
                                      
!     call strain%init(xx=0.001D0, yy=-0.0003D0, zz=-0.0003D0, xy=0D0, yz=0D0, xz=0D0)
!     call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
!     call elas%set_parameters(young=1000D0, poisson=0.3D0)
!     sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

!     call data%init(strain_pf=strain_pf, strain_p=strain_p)
!     call solver%init(elasticity=elas, hardening=sw, yield=vm)
!     call solver%solve(strain=strain, data=data)
!     call data%get(stress=stress,        & 
!                   strain_pf=strain_pf,  &
!                   strain_p=strain_p     &
!                   )

!     passed = stress .approx. expected_stress
!     if (.not. passed) print*, "Stress is no equal", new_line('A'),          &
!                               "Expected:", expected_stress, new_line('A'),  &
!                               "Actual Value:", stress, new_line('A'),       &
!                               "Difference", stress - expected_stress
!     if (.not. passed) return
    
!     passed = abs(strain_pf - expected_strain_effective) < EPS
!     if (.not. passed) print*, "Effective plastic Strain is not equal", new_line('A'), &
!                               "Expected:", expected_strain_effective, new_line('A'),  &
!                               "Actual Value:", strain_pf, new_line('A'),              &
!                               "Difference", strain_pf - expected_strain_effective
!     if (.not. passed) return

!     passed = strain_p .approx. expected_strain_plastic
!     if (.not. passed) print*, "Plastic Strain is not equal", new_line('A'), &
!                               "Expected:", expected_strain_plastic, new_line('A'), &
!                               "Actual Value:", strain_p, new_line('A'),            & 
!                               "Difference:", strain_p - expected_strain_plastic
!     if (.not. passed) return


!     call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
!     call solver%tangent(strain=strain, data=data, tangent=tangent)

!     passed = tangent .approx. numerical_tangent
!     if (.not. passed) print*, "Tangent is no equal", new_line('A'),          &
!                               "Analitical:", tangent, new_line('A'),  &
!                               "Numerical:", numerical_tangent, new_line('A'),       &
!                               "Difference", tangent - numerical_tangent
!     if (.not. passed) return

!     return
! end subroutine
subroutine test_closest_point_time_reinit_other_types(passed)
    ! A second init with Elasticity_linear, a larger extension of the JC law and Hill48 after
    ! Elasticity_frozen, JC and VonMises must give the same step as a solver initialized with the
    ! second set only. Without reallocating each component, the new data overflow the old storage.
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_vp_jc, only : JC_viscoplastic
    use muscle_hard_swift, only : Swift_hardening
    use muscle_yield_vonmises, only : VonMises
    use muscle_yield_hill48, only : Hill48
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_visco_frozen, only : Elasticity_frozen
    use muscle_solver_closest_point_time
    use muscle_test_utils, only : check_approx, check_true
    implicit none

    type, extends(JC_viscoplastic) :: Larger_jc
        real(real64) :: extra(32) = 0D0
    end type Larger_jc

    logical, intent(out) :: passed
    type(Closest_point_data) :: data, data_ref
    type(Closest_point) :: solver, fresh
    type(VonMises) :: vm
    type(Hill48) :: hill
    type(ten_3D2Osym) :: strain, strain_p, stress, stress_ref, strain_p_ref
    type(JC_viscoplastic) :: jc
    type(Larger_jc) :: larger
    type(Elasticity_frozen) :: frozen
    type(Elasticity_linear) :: lin
    type(ten_3D4O2sym) :: tangent, tangent_ref
    real(real64) :: strain_pf, strain_pf_ref
    integer :: status, status_ref

    passed = .true.
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call frozen%set(k_alg=1000D0/1.2D0, g_alg=1000D0/2.6D0, stress_hist=strain_p)
    call lin%set_parameters(young=1000D0, poisson=0.3D0)
    hill = Hill48(f=0.25D0, g=0.35D0, h=0.65D0, l=1.3D0, m=1.7D0, n=1.45D0)
    jc = JC_viscoplastic(C=1D-1, epdmax=1.0D-3)
    jc%hard_law = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    larger%C = 1D-1
    larger%epdmax = 1.0D-3
    larger%hard_law = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call solver%init(elasticity=frozen, hardening=jc, yield=vm)
    call solver%init(elasticity=lin, hardening=larger, yield=hill)
    call fresh%init(elasticity=lin, hardening=larger, yield=hill)
    call strain%init(xx=0.08D0, yy=-0.024D0, zz=-0.01D0, xy=0.012D0, yz=-0.007D0, xz=0.004D0)
    call data%init(strain_pf=0D0, strain_p=strain_p, dt=1.0D+2)
    call data_ref%init(strain_pf=0D0, strain_p=strain_p, dt=1.0D+2)
    call solver%solve(strain=strain, data=data)
    call fresh%solve(strain=strain, data=data_ref)
    call data%get(stress=stress, strain_pf=strain_pf, strain_p=strain_p, status=status)
    call data_ref%get(stress=stress_ref, strain_pf=strain_pf_ref, strain_p=strain_p_ref, status=status_ref)
    call check_true(passed, "reinit: converged", status == STATUS_CONVERGED)
    call check_true(passed, "fresh: converged", status_ref == STATUS_CONVERGED)
    call check_true(passed, "reinit: plastic step", strain_pf > 0D0)

    ! Same operations on the same data; Measured: 0
    call solver%tangent(strain=strain, data=data, tangent=tangent)
    call fresh%tangent(strain=strain, data=data_ref, tangent=tangent_ref)
    call check_approx(passed, "reinit: stress", stress, stress_ref, tol=0D0)
    call check_approx(passed, "reinit: plastic strain", strain_p, strain_p_ref, tol=0D0)
    call check_approx(passed, "reinit: equivalent plastic strain", strain_pf, strain_pf_ref, tol=0D0)
    call check_approx(passed, "reinit: tangent", tangent, tangent_ref, tol=0D0)
end subroutine


subroutine test_closest_point_time_vonmises_mrk_uniaxial_tensile(passed)
    !! Uniaxial MRK tension from a previous plastic state: the flow stress uses the total
    !! equivalent plastic strain, the rate uses only the step increment (ep - ep_n)/dt, and the
    !! consistent tangent covers both derivatives.
    ! Von Mises with associated flow and MRK as the flow stress. Step from ep_n = 0.03 to 0.05 in
    ! dt = 80 s, so the rate 2.5D-4 1/s sits inside the rate domain where nu1 and nu2 both act. A
    ! rate built with the total ep instead of the increment would raise the stress by 1.06 MPa
    ! (3.9D-3), which this case would catch. Expected values evaluated independently with mpmath
    ! (60 digits), not with this solver.
    ! Called by test_muscle_solver_closest_point_time, STOP 6.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_test_utils, only : check_approx, check_true
    use muscle_vp_mrk, only : MRK_viscoplastic
    use muscle_yield_vonmises, only : VonMises
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_solver_closest_point_time
    implicit none
    logical, intent(out) :: passed
    ! Steel, so the elastic part is small next to the plastic part and the stress is not fitted.
    ! Measured: 1.1D-8 or less on stress and plastic strain.
    real(real64), parameter :: TOL_STRESS = 1D-6
    ! Absolute equivalent plastic strain error, scalar scale = 1. Measured: 1.8D-11.
    real(real64), parameter :: TOL_PF = 1D-8
    ! Both derivatives of MRK enter the tangent; the numerical one is a central difference of the
    ! same solve, so the tolerance only has to absorb its truncation error. Measured: 1.5D-9.
    real(real64), parameter :: TOL_TANGENT = 1D-7
    ! Newton bound for a smooth point inside the rate domain: the JC case of this file uses 6.
    ! Measured: 3 iterations.
    integer, parameter :: MAX_ITERS = 6
    type(Closest_point) :: solver
    type(Closest_point_data) :: data
    type(VonMises) :: vm
    type(MRK_viscoplastic) :: mrk
    type(Elasticity_linear) :: elas
    type(ten_3D2Osym) :: strain, stress, strain_p, expected_stress, expected_strain_p
    type(ten_3D4O2sym) :: tangent, numerical_tangent
    real(real64) :: strain_pf
    integer :: status, iters

    passed = .true.
    ! Prior state: ep_n = 0.03 with the transverse plastic strains at -ep_n/2.
    call expected_stress%init((/ 2.73922882128596387D+02, 0D0, 0D0, 0D0, 0D0, 0D0 /))
    call expected_strain_p%init((/ 5.00000000000000000D-02, -2.50000000000000000D-02, &
                                  -2.50000000000000000D-02, 0D0, 0D0, 0D0 /))
    call strain%init((/ 5.13043946768028399D-02, -2.53913184030408520D-02, &
                       -2.53913184030408520D-02, 0D0, 0D0, 0D0 /))
    call strain_p%init((/ 0.03D0, -0.015D0, -0.015D0, 0D0, 0D0, 0D0 /))
    mrk = MRK_viscoplastic(B01=100.020701D0, B02=100D0, epdmax=0.01D0, &
                           nu1=0.01439459D0, nu2=0.14313556D0, n0=3.95295384D0, D2=5.6375D-12, &
                           epdmin=1D-5, sig_u=253.772550D0, chi1=9.1503D-4, chi2=0.01942062D0)
    call elas%set_parameters(young=210000D0, poisson=0.3D0)
    call solver%init(elasticity=elas, hardening=mrk, yield=vm)
    call data%init(strain_pf=0.03D0, strain_p=strain_p, dt=80D0)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress, strain_p=strain_p, strain_pf=strain_pf, status=status, iters=iters)
    print *, 'von Mises/MRK uniaxial status:', status, ' iterations:', iters
    ! A plastic step must finish converged, not elastic, not exhausted.
    call check_true(passed, "Converged status", status == STATUS_CONVERGED)
    if (status /= STATUS_CONVERGED) return
    call check_approx(passed, "Stress", stress, expected_stress, tol=TOL_STRESS)
    ! The scalar adds the step increment to the previous 0.03.
    call check_approx(passed, "Equivalent plastic strain", strain_pf, 0.05D0, tol=TOL_PF)
    ! Associated Von Mises flow gives transverse plastic strains of -ep/2.
    call check_approx(passed, "Plastic strain", strain_p, expected_strain_p, tol=TOL_STRESS)
    ! Bound the Newton work for this smooth point inside the rate domain.
    call check_true(passed, "Iteration bound", iters <= MAX_ITERS)
    ! The consistent tangent must include both the strain and the rate derivative of MRK.
    call solver%tangent(strain=strain, data=data, tangent=tangent)
    call solver%tangent_numerical(strain=strain, data=data, tangent=numerical_tangent)
    call check_approx(passed, "Consistent tangent", tangent, numerical_tangent, tol=TOL_TANGENT)
end subroutine test_closest_point_time_vonmises_mrk_uniaxial_tensile
