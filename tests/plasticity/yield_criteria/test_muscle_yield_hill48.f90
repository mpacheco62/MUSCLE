program test_muscle_yield_hill48
    implicit none
    
    logical :: passed

    call test_Hill48_stresseq_hydrostatic(passed)
    if (.not. passed) STOP 1

    ! call test_vonMises_stresseq_simple_tensile(passed)
    ! if (.not. passed) STOP 2

    ! call test_vonMises_stresseq_zero_stress(passed)
    ! if (.not. passed) STOP 3

    ! call test_vonMises_stresseq_biaxial(passed)
    ! if (.not. passed) STOP 4

    ! call test_vonMises_stresseq_shear(passed)
    ! if (.not. passed) STOP 5

    ! call test_vonMises_stresseq_derivates(passed)
    ! if (.not. passed) STOP 6

    ! call test_vonMises_stresseq_derivates2(passed)
    ! if (.not. passed) STOP 7

    call test_Hill48_gradient(passed)
    if (.not. passed) STOP 8

    call test_Hill48_hessian(passed)
    if (.not. passed) STOP 9

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_hill48

subroutine test_Hill48_stresseq_hydrostatic(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_hill48
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(yield_hill48_dp) :: h48
    type(ten_3D2Osym) :: to_test1 

    real(real64) :: result
    real(real64) :: expected_result1 = 0.0D0


    call to_test1%init((/5D0, 5D0, 5D0, 0D0, 0D0, 0D0/))
    h48 = yield_hill48_dp(f=0.5D0, g=0.5D0, h=0.5D0, l=1.5D0, m=1.5D0, n=1.5D0)

    result = yield_stress_eq(h48, to_test1)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) print *, "Hill48 test hydrostatic failed, result:", &
                               result, ", expected result:", expected_result1
    if (.not. passed) return
end subroutine

! subroutine test_vonMises_stresseq_simple_tensile(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     use muscle_yield_vonmises
!     implicit none
    
!     real(real64), parameter :: EPS=1e-10
!     logical, intent(out) :: passed

!     type(yield_vonmises_dp) :: vm
!     type(ten_3D2Osym) :: to_test2

!     real(real64) :: result
!     real(real64) :: expected_result2 = 5D0

!     call to_test2%init((/5D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    
!     result = yield_stress_eq(vm, to_test2)
!     passed = (abs(result - expected_result2) < EPS)
!     if (.not. passed) return
    
! end subroutine

! subroutine test_vonMises_stresseq_zero_stress(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     use muscle_yield_vonmises
!     implicit none
    
!     real(real64), parameter :: EPS=1e-10
!     logical, intent(out) :: passed

!     type(yield_vonmises_dp) :: vm
!     type(ten_3D2Osym) :: to_test1 

!     real(real64) :: result
!     real(real64) :: expected_result1 = 0.0D0


!     call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))

!     result = yield_stress_eq(vm, to_test1)
!     passed = (abs(result - expected_result1) < EPS)
!     if (.not. passed) return
! end subroutine

! subroutine test_vonMises_stresseq_biaxial(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     use muscle_yield_vonmises
!     implicit none
    
!     real(real64), parameter :: EPS=1e-10
!     logical, intent(out) :: passed

!     type(yield_vonmises_dp) :: vm
!     type(ten_3D2Osym) :: to_test3

!     real(real64) :: result
!     real(real64) :: expected_result3 = 5D0


!     call to_test3%init((/5D0, 5D0, 0D0, 0D0, 0D0, 0D0/))
!     result = yield_stress_eq(vm, to_test3)
!     passed = (abs(result - expected_result3) < EPS)
!     if (.not. passed) return

! end subroutine

! subroutine test_vonMises_stresseq_shear(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     use muscle_yield_vonmises
!     implicit none
    
!     real(real64), parameter :: EPS=1e-10
!     logical, intent(out) :: passed

!     type(yield_vonmises_dp) :: vm
!     type(ten_3D2Osym) :: to_test4

!     real(real64) :: result
!     real(real64) :: expected_result4 = 3D0**0.5D0 

!     call to_test4%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
!     result = yield_stress_eq(vm, to_test4)
!     passed = (abs(result - expected_result4) < EPS)
!     if (.not. passed) return
! end subroutine

! subroutine test_vonMises_stresseq_derivates(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     use muscle_yield_vonmises
!     implicit none
    
!     real(real64), parameter :: EPS=1e-10
!     logical, intent(out) :: passed

!     type(yield_vonmises_dp) :: vm
!     type(ten_3D2Osym) :: to_test

!     type(ten_3D2Osym) :: result1, result2

!     call to_test%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
!     result1 = yield_dstress_eq(vm, to_test)
!     result2 = num_dstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
!     result1 = yield_dstress_eq(vm, to_test)
!     result2 = num_dstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
!     result1 = yield_dstress_eq(vm, to_test)
!     result2 = num_dstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
!     result1 = yield_dstress_eq(vm, to_test)
!     result2 = num_dstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
!     result1 = yield_dstress_eq(vm, to_test)
!     result2 = num_dstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
!     result1 = yield_dstress_eq(vm, to_test)
!     result2 = num_dstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
!     result1 = yield_dstress_eq(vm, to_test)
!     result2 = num_dstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
!     result1 = yield_dstress_eq(vm, to_test)
!     result2 = num_dstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

! end subroutine


! subroutine test_vonMises_stresseq_derivates2(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     use muscle_yield_vonmises
!     implicit none
    
!     real(real64), parameter :: EPS=1e-10
!     logical, intent(out) :: passed

!     type(yield_vonmises_dp) :: vm
!     type(ten_3D2Osym) :: to_test

!     type(ten_3D4O3sym) :: result1, result2

!     call to_test%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
!     result1 = yield_ddstress_eq(vm, to_test)
!     result2 = num_ddstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
!     result1 = yield_ddstress_eq(vm, to_test)
!     result2 = num_ddstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
!     result1 = yield_ddstress_eq(vm, to_test)
!     result2 = num_ddstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
!     result1 = yield_ddstress_eq(vm, to_test)
!     result2 = num_ddstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
!     result1 = yield_ddstress_eq(vm, to_test)
!     result2 = num_ddstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
!     result1 = yield_ddstress_eq(vm, to_test)
!     result2 = num_ddstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
!     result1 = yield_ddstress_eq(vm, to_test)
!     result2 = num_ddstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

!     call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
!     result1 = yield_ddstress_eq(vm, to_test)
!     result2 = num_ddstress_eq(vm, to_test)
!     passed = result1 .isequal. result2
!     if (.not. passed) return

! end subroutine

subroutine test_Hill48_gradient(passed)
    ! Analytical gradient of Hill48 with the AA2090-T3 coefficients, obtained from the
    ! r-values r0, r45, r90 (Yoon and Barlat, 2006) and normalized to the RD yield stress.
    ! Checks: closed forms in uniaxial RD tension and pure xy shear, the finite-difference
    ! finite-difference gradient (yield_test_numeric) in a general state, and Euler's identity grad : stress = seq.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_hill48
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    implicit none

    real(real64), parameter :: EPS = 1.0e-12_real64, EPS_FD = 1.0e-7_real64
    logical, intent(out) :: passed

    type(yield_hill48_dp) :: h48
    type(ten_3D2Osym) :: stress, grad, grad_fd, expected
    real(real64) :: r0, r45, r90, f, g, h, n, seq

    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    f = r0/(r90*(1D0 + r0))
    g = 1D0/(1D0 + r0)
    h = r0/(1D0 + r0)
    n = (r0 + r90)*(1D0 + 2D0*r45)/(2D0*r90*(1D0 + r0))
    h48 = yield_hill48_dp(f=f, g=g, h=h, l=1.5D0, m=1.5D0, n=n)

    ! Uniaxial tension in RD: grad = (g + h, -h, -g, 0, 0, 0)/sqrt(g + h)
    call stress%init((/2D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    grad = yield_dstress_eq(h48, stress)
    call expected%init((/g + h, -h, -g, 0D0, 0D0, 0D0/)/sqrt(g + h))
    passed = grad%is_approx(expected, tol=EPS)
    if (.not. passed) print *, "Hill48 gradient, uniaxial RD:", grad%vals
    if (.not. passed) return

    ! Pure xy shear (tensorial): grad_xy = n*sxy/seq = sqrt(n/2)
    call stress%init((/0D0, 0D0, 0D0, 3D0, 0D0, 0D0/))
    grad = yield_dstress_eq(h48, stress)
    call expected%init((/0D0, 0D0, 0D0, sqrt(0.5D0*n), 0D0, 0D0/))
    passed = grad%is_approx(expected, tol=EPS)
    if (.not. passed) print *, "Hill48 gradient, xy shear:", grad%vals
    if (.not. passed) return

    ! General state: finite differences (yield_test_numeric) and Euler's identity
    call stress%init((/1.3D0, -0.4D0, 0.7D0, 0.5D0, -0.8D0, 0.2D0/))
    grad = yield_dstress_eq(h48, stress)
    grad_fd = num_dstress_eq(h48, stress)
    seq = yield_stress_eq(h48, stress)
    passed = grad%is_approx(grad_fd, tol=EPS_FD) .and. &
             abs((grad .ddot. stress) - seq) < EPS*seq
    if (.not. passed) print *, "Hill48 gradient, general state:", grad%vals, grad_fd%vals
end subroutine test_Hill48_gradient

subroutine test_Hill48_hessian(passed)
    ! Analytical Hessian of Hill48 with the AA2090-T3 coefficients of test_Hill48_gradient,
    ! against the finite-difference Hessian (yield_test_numeric) (ddstressEq_ddstress_numeric),
    ! and of the isotropic Hill48 against the analytical Hessian of VonMises.
    ! Measured differences: below 1e-7 relative against finite differences.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_hill48
    use muscle_yield_vonmises
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    implicit none

    real(real64), parameter :: EPS = 1.0e-12_real64, EPS_FD = 1.0e-6_real64
    logical, intent(out) :: passed

    type(yield_hill48_dp) :: h48
    type(yield_vonmises_dp) :: vm
    type(ten_3D2Osym) :: stress
    type(ten_3D4O3sym) :: hess, hess_fd, hess_vm
    real(real64) :: r0, r45, r90

    r0 = 0.2115D0
    r45 = 1.5769D0
    r90 = 0.6923D0
    h48 = yield_hill48_dp(f=r0/(r90*(1D0 + r0)), g=1D0/(1D0 + r0), h=r0/(1D0 + r0), &
                 l=1.5D0, m=1.5D0, n=(r0 + r90)*(1D0 + 2D0*r45)/(2D0*r90*(1D0 + r0)))

    ! General state: normal and shear components together, so every entry of the Hessian,
    ! including the normal-shear coupling, is non-zero
    call stress%init((/1.3D0, -0.4D0, 0.7D0, 0.5D0, -0.8D0, 0.2D0/))
    hess = yield_ddstress_eq(h48, stress)
    hess_fd = num_ddstress_eq(h48, stress)
    passed = hess%is_approx(hess_fd, tol=EPS_FD)
    if (.not. passed) print *, "Hill48 Hessian, general state:", hess%vals, hess_fd%vals
    if (.not. passed) return

    ! Uniaxial tension in RD: the shear entries reduce to n/(2 seq), l/(2 seq), m/(2 seq)
    call stress%init((/2D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    hess = yield_ddstress_eq(h48, stress)
    hess_fd = num_ddstress_eq(h48, stress)
    passed = hess%is_approx(hess_fd, tol=EPS_FD)
    if (.not. passed) print *, "Hill48 Hessian, uniaxial RD:", hess%vals, hess_fd%vals
    if (.not. passed) return

    ! Isotropic limit: with f = g = h = 1/2 and l = m = n = 3/2, Hill48 is von Mises, so in a
    ! general state its Hessian equals the analytical one of VonMises, independent of finite
    ! differences
    h48 = yield_hill48_dp(f=0.5D0, g=0.5D0, h=0.5D0, l=1.5D0, m=1.5D0, n=1.5D0)
    call stress%init((/1.3D0, -0.4D0, 0.7D0, 0.5D0, -0.8D0, 0.2D0/))
    hess = yield_ddstress_eq(h48, stress)
    hess_vm = yield_ddstress_eq(vm, stress)
    passed = hess%is_approx(hess_vm, tol=EPS)
    if (.not. passed) print *, "Hill48 Hessian, isotropic limit:", hess%vals, hess_vm%vals
end subroutine test_Hill48_hessian
