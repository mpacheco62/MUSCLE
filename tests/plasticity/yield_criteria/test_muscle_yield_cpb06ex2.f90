! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_cpb06ex2
    implicit none
    logical :: passed
    call test_CPB06ex2_reduces_to_CPB06(passed)
    if (.not. passed) STOP 1
    call test_CPB06ex2_isotropic_vonmises(passed)
    if (.not. passed) STOP 2
    call test_CPB06ex2_paper_tables(passed)
    if (.not. passed) STOP 3
    call test_CPB06ex2_tension_compression_asymmetry(passed)
    if (.not. passed) STOP 4
    call test_CPB06ex2_pressure_insensitive(passed)
    if (.not. passed) STOP 5
    call test_CPB06ex2_gradient_vs_numeric(passed)
    if (.not. passed) STOP 6
    call test_CPB06ex2_gradient_analytical(passed)
    if (.not. passed) STOP 7
    call test_CPB06ex2_hessian_vs_numeric(passed)
    if (.not. passed) STOP 8
    call test_CPB06ex2_hessian_analytical(passed)
    if (.not. passed) STOP 9
    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_cpb06ex2

subroutine test_CPB06ex2_reduces_to_CPB06(passed)
    ! With Cp = C and kp = k (the primed map), CPB06ex2 is CPB06: both sums are equal and
    ! normalized the same way.
    ! General nonsymmetric normal block and three distinct shear factors expose map errors.
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 1)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    logical, intent(out) :: passed
    ! Same eigenvalues and powers in both models: round-off only. Measured: 1.6e-16
    real(real64), parameter :: TOL = 1D-12
    type(yield_cpb06ex2_dp) :: yld
    type(yield_cpb06_dp) :: cpb
    type(ten_3D2Osym) :: s
    passed = .true.
    call yield_cpb06_init_dp(cpb, 1D0,-0.2D0,0.3D0,0.15D0,0.8D0,-0.4D0,-0.1D0,0.2D0,1.2D0, &
                  1.3D0,0.7D0,1.1D0,0.25D0,4D0, status=status)
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06_init failed"
        return
    end if
    call yield_cpb06ex2_init_dp(yld, 1D0,-0.2D0,0.3D0,0.15D0,0.8D0,-0.4D0,-0.1D0,0.2D0,1.2D0, &
                  1.3D0,0.7D0,1.1D0,0.25D0, &
                  1D0,-0.2D0,0.3D0,0.15D0,0.8D0,-0.4D0,-0.1D0,0.2D0,1.2D0, &
                  1.3D0,0.7D0,1.1D0,0.25D0,4D0, status=status)
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06ex2_init failed"
        return
    end if
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, 'CPB06 reduction coupled', yield_stress_eq(yld, s), yield_stress_eq(cpb, s), &
                      tol=TOL)
    s = -s
    call check_approx(passed, 'CPB06 reduction compression', yield_stress_eq(yld, s), yield_stress_eq(cpb, s), &
                      tol=TOL)
end subroutine test_CPB06ex2_reduces_to_CPB06

subroutine test_CPB06ex2_isotropic_vonmises(passed)
    ! Isotropic a=2 checks all tensorial shear weights against the closed-form VonMises.
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 2)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_vonmises
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    logical, intent(out) :: passed
    ! Closed form against eigenvalues: round-off only. Measured: 1.7e-16
    real(real64), parameter :: TOL = 1D-12
    type(yield_cpb06ex2_dp) :: yld
    type(yield_vonmises_dp) :: vm
    type(ten_3D2Osym) :: s
    passed = .true.
    call yield_cpb06ex2_init_dp(yld, 1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0,1D0,1D0,1D0,0D0, &
                  1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0,1D0,1D0,1D0,0D0,2D0, status=status)
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06ex2_init failed"
        return
    end if
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, 'VonMises coupled', yield_stress_eq(yld, s), yield_stress_eq(vm, s), tol=TOL)
end subroutine test_CPB06ex2_isotropic_vonmises

subroutine CPB06ex2_test_materials(yld)
    ! Plunkett et al. (2008), Tables 2, 8, 1: steel, 6111-T4 and AZ31B, respectively.
    ! The shear order yz,xz,xy of the paper is permuted to library xy,yz,xz.
    ! AZ31B lacks out-of-plane shear data: unit factors there are synthetic test inputs.
    ! Called by test_CPB06ex2_paper_tables, test_CPB06ex2_tension_compression_asymmetry,
    ! test_CPB06ex2_pressure_insensitive and the gradient and Hessian tests
    use muscle_yield_cpb06ex2
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    type(yield_cpb06ex2_dp), intent(out) :: yld(3)
    call yield_cpb06ex2_init_dp(yld(1), 1.283D0,0.208D0,-0.421D0,0.208D0,0.535D0,-0.791D0,-0.421D0,-0.791D0,0.746D0, &
                     1.156D0,1.418D0,1.329D0,0.030D0, &
                     1.283D0,0.380D0,-0.001D0,0.380D0,2.145D0,0.792D0,-0.001D0,0.792D0,0.218D0, &
                     1.721D0,1.418D0,1.329D0,-0.027D0,5D0, status=status)
    if (status /= STATUS_OK) error stop 'yield_cpb06ex2_init failed'
    call yield_cpb06ex2_init_dp(yld(2), 0.250D0,-1.399D0,-0.971D0,-1.399D0,-0.375D0,-2.254D0, &
                     -0.971D0,-2.254D0,-0.831D0, &
                     1.316D0,1.403D0,1.101D0,0D0, &
                     0.250D0,-0.840D0,-1.330D0,-0.840D0,0.780D0,-0.565D0, &
                     -1.330D0,-0.565D0,0.479D0, &
                     1.413D0,1.094D0,1.448D0,0D0,12D0, status=status)
    if (status /= STATUS_OK) error stop 'yield_cpb06ex2_init failed'
    call yield_cpb06ex2_init_dp(yld(3), 0.771D0,-0.914D0,-1.312D0,-0.914D0,0.654D0,-1.209D0, &
                     -1.312D0,-1.209D0,-0.293D0, &
                     1.617D0,1D0,1D0,0.295D0, &
                     0.771D0,0.337D0,-0.764D0,0.337D0,0.646D0,-0.401D0,-0.764D0,-0.401D0,-0.022D0, &
                     0.499D0,1D0,1D0,0.900D0,4D0, status=status)
    if (status /= STATUS_OK) error stop 'yield_cpb06ex2_init failed'
end subroutine CPB06ex2_test_materials

subroutine test_CPB06ex2_paper_tables(passed)
    ! Rounded steel parameters and out-of-plane 6111-T4 loadings detect the shear permutation.
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 3)
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(yield_cpb06ex2_dp) :: yld(3)
    type(ten_3D2Osym) :: s
    passed = .true.
    call CPB06ex2_test_materials(yld)
    ! Table 2 has three-digit parameters: 2e-3 allows their rounding, not solver error.
    ! Measured: 9.7e-4 or less (Y_C)
    call s%init((/1D0,0D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 X_T', 1D0/yield_stress_eq(yld(1), s), 1D0, tol=2D-3)
    call s%init((/0D0,1D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 Y_T', 1D0/yield_stress_eq(yld(1), s), 0.914D0, tol=2D-3)
    call s%init((/-1D0,0D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 X_C', 1D0/yield_stress_eq(yld(1), s), 0.963D0, tol=2D-3)
    call s%init((/0D0,-1D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 Y_C', 1D0/yield_stress_eq(yld(1), s), 0.960D0, tol=2D-3)
    call s%init((/0D0,0D0,-1D0,0D0,0D0,0D0/))
    call check_approx(passed, 'steel Table 3 Z_C', 1D0/yield_stress_eq(yld(1), s), 0.958D0, tol=2D-3)
    ! Table 9 has two decimals: half a last digit plus rounding of Table 8. Measured: 2.6e-4 or less
    call s%init((/0D0,0.5D0,0.5D0,0D0,0.5D0,0D0/))
    call check_approx(passed, '6111 Table 9 YZ 45', 1D0/yield_stress_eq(yld(2), s), 1.07D0, tol=6D-3)
    call s%init((/0.5D0,0D0,0.5D0,0D0,0D0,0.5D0/))
    call check_approx(passed, '6111 Table 9 ZX 45', 1D0/yield_stress_eq(yld(2), s), 1.04D0, tol=6D-3)
    call s%init((/0D0,0D0,0D0,0D0,1D0,0D0/))
    call check_approx(passed, '6111 Table 9 YZ shear', 1D0/yield_stress_eq(yld(2), s), 0.67D0, tol=6D-3)
    call s%init((/0D0,0D0,0D0,0D0,0D0,1D0/))
    call check_approx(passed, '6111 Table 9 ZX shear', 1D0/yield_stress_eq(yld(2), s), 0.65D0, tol=6D-3)
    ! Independent 90-digit trace-of-powers reference, all normal and shear components coupled.
    ! Measured: 1.7e-16
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, '6111 independent scalar', yield_stress_eq(yld(2), s), &
                      170.33623760110963846252D0, tol=1D-12)
end subroutine test_CPB06ex2_paper_tables

subroutine test_CPB06ex2_tension_compression_asymmetry(passed)
    ! kp = 0.9 (the primed k) makes AZ31B compression much weaker than tension in RD (Table 1).
    ! The RD normalization of init is exact: f = 1 in unit tension along RD (measured: 0).
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 4)
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(yield_cpb06ex2_dp) :: yld(3)
    type(ten_3D2Osym) :: s
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call s%init((/1D0,0D0,0D0,0D0,0D0,0D0/))
    call check_approx(passed, 'AZ31 RD normalization', yield_stress_eq(yld(3), s), 1D0, tol=1D-12)
    call check_true(passed, 'AZ31 tension compression asymmetry', &
                    yield_stress_eq(yld(3), -s) > 1.5D0*yield_stress_eq(yld(3), s))
end subroutine test_CPB06ex2_tension_compression_asymmetry

subroutine test_CPB06ex2_pressure_insensitive(passed)
    ! Both maps act on the deviator (Eq. 9): a hydrostatic stress gives zero, and adding a
    ! pressure to a coupled state leaves the equivalent stress unchanged. AZ31B, kp = 0.9.
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 5)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! linear_transform subtracts the zz stress before the products, so the hydrostatic state
    ! maps to an exact zero and the shifted state differs only by the rounding of the
    ! differences. Measured: exact (0) in both checks
    real(real64), parameter :: TOL = 1D-12
    type(yield_cpb06ex2_dp) :: yld(3)
    type(ten_3D2Osym) :: s, p
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call p%init((/250D0,250D0,250D0,0D0,0D0,0D0/))
    call check_approx(passed, 'AZ31 hydrostatic stress = 0', yield_stress_eq(yld(3), p), 0D0, tol=TOL)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, 'AZ31 coupled state plus pressure 250', yield_stress_eq(yld(3), s + p), &
                      yield_stress_eq(yld(3), s), tol=TOL)
end subroutine test_CPB06ex2_pressure_insensitive

subroutine test_CPB06ex2_gradient_vs_numeric(passed)
    ! Gradient against the finite-difference gradient (yield_test_numeric), for the three
    ! materials in a state with normal and shear components together, and for 6111-T4 in a
    ! rotated state where two principal values of Sigma differ by about 1e-8 relatively
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 6)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! Central differences (yield_test_numeric): truncation and the rounding of eigenvals.
    ! Measured: 1.2e-10 or less
    real(real64), parameter :: TOL = 1D-8
    type(yield_cpb06ex2_dp) :: yld(3)
    type(ten_3D2Osym) :: s
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    ! Steel: k and kp of opposite sign, a = 5
    call check_approx(passed, 'gradient vs numeric: steel', yield_dstress_eq(yld(1), s), &
                      num_dstress_eq(yld(1), s), tol=TOL)
    ! 6111-T4: a = 12, the largest powers of the three
    call check_approx(passed, 'gradient vs numeric: 6111-T4', yield_dstress_eq(yld(2), s), &
                      num_dstress_eq(yld(2), s), tol=TOL)
    ! AZ31B: kp = 0.9, the second transformation almost one-sided
    call check_approx(passed, 'gradient vs numeric: AZ31B', yield_dstress_eq(yld(3), s), &
                      num_dstress_eq(yld(3), s), tol=TOL)
    call s%init((/197.95561360136386D0,-58.04950942270237D0,-139.90610417866154D0, &
                  -50.24571382486021D0,11.771438358965758D0,-51.20999705196427D0/))
    call check_approx(passed, 'gradient vs numeric: 6111-T4, close pair', &
                      yield_dstress_eq(yld(2), s), num_dstress_eq(yld(2), s), tol=TOL)
end subroutine test_CPB06ex2_gradient_vs_numeric

subroutine test_CPB06ex2_gradient_analytical(passed)
    ! Gradient of 6111-T4 against Euler's identity for the homogeneous f (grad : stress = f)
    ! and against an independent 90-digit reference computed with traces of powers of
    ! Sigma and Sigmap (no eigenvalues), in a state with all components
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 7)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! Exact references: round-off only. Measured: 2.8e-16 or less
    real(real64), parameter :: TOL = 1D-12
    type(yield_cpb06ex2_dp) :: yld(3)
    type(ten_3D2Osym) :: s, expected
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call check_approx(passed, 'gradient Euler identity: 6111-T4', &
                      yield_dstress_eq(yld(2), s) .ddot. s, yield_stress_eq(yld(2), s), tol=TOL)
    call expected%init((/0.63491427307714846058430D0, -0.8014554330786033866660D0, &
                         0.16654116000145492608174D0, 0.29454772703368041095476D0, &
                         -0.2025437658446382337026D0, 0.28693805299310551585208D0/))
    call check_approx(passed, 'gradient vs reference: 6111-T4', yield_dstress_eq(yld(2), s), &
                      expected, tol=TOL)
end subroutine test_CPB06ex2_gradient_analytical

subroutine test_CPB06ex2_hessian_vs_numeric(passed)
    ! Hessian against the finite-difference Hessian (yield_test_numeric), for the three
    ! materials in a state with normal and shear components together (it fills the
    ! normal-shear entries and the mixed block between Sigma and Sigmap)
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 8)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    ! Central second differences (yield_test_numeric): rounding of eps*f/h**2 and a truncation
    ! O(h**2) that grows with a. Measured: 2.3e-7 or less
    real(real64), parameter :: TOL = 1D-6
    type(yield_cpb06ex2_dp) :: yld(3)
    type(ten_3D2Osym) :: s
    passed = .true.
    call CPB06ex2_test_materials(yld)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    ! Steel: k and kp of opposite sign, a = 5
    call check_approx(passed, 'Hessian vs numeric: steel', yield_ddstress_eq(yld(1), s), &
                      num_ddstress_eq(yld(1), s), tol=TOL)
    ! 6111-T4: a = 12, the largest powers of the three
    call check_approx(passed, 'Hessian vs numeric: 6111-T4', yield_ddstress_eq(yld(2), s), &
                      num_ddstress_eq(yld(2), s), tol=TOL)
    ! AZ31B: kp = 0.9, the second transformation almost one-sided
    call check_approx(passed, 'Hessian vs numeric: AZ31B', yield_ddstress_eq(yld(3), s), &
                      num_ddstress_eq(yld(3), s), tol=TOL)
end subroutine test_CPB06ex2_hessian_vs_numeric

subroutine test_CPB06ex2_hessian_analytical(passed)
    ! Hessian against analytical values: the 21 components for 6111-T4 against an independent
    ! 90-digit reference, von Mises in the isotropic limit at a zero principal value, and the
    ! mean of the one-sided curvatures at a zero principal value with k, kp /= 0
    ! Called by program test_muscle_yield_cpb06ex2 (STOP 9)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06ex2
    use muscle_yield_vonmises
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    logical, intent(out) :: passed
    ! Exact references: round-off only. Measured: 8.8e-16 or less
    real(real64), parameter :: TOL = 1D-12
    type(yield_cpb06ex2_dp) :: yld(3), iso, az
    type(yield_vonmises_dp) :: vm
    type(ten_3D2Osym) :: s
    type(ten_3D4O3sym) :: expected
    passed = .true.
    call CPB06ex2_test_materials(yld)
    ! 6111-T4, all components: reference computed with traces of powers of Sigma and Sigmap
    ! in 90-digit arithmetic (no eigenvalues)
    call s%init((/120.5D0,-35.25D0,60.75D0,42.5D0,-28D0,33.25D0/))
    call expected%init((/3.448675381240494907D-03, 2.429944194238824531D-03, &
                         5.686182425050850114D-03, 3.890903566445265936D-03, &
                         4.758349797967605446D-03, 3.843197135993091807D-03, &
                         -9.621857521423459413D-05, -2.333725619024589869D-03, &
                         6.862937982881940494D-04, -6.260620910707611047D-04, &
                         1.638792459368608933D-03, -3.352456806026260245D-03, &
                         1.406435033139943414D-03, -3.900271823623490557D-04, &
                         -1.589894276519464777D-03, -2.092728831428137572D-03, &
                         -1.203204133117417298D-03, -1.562479102753076962D-03, &
                         1.593231315479766300D-03, 7.834179435264228015D-04, &
                         7.790611592266542687D-04/))
    call check_approx(passed, 'Hessian vs reference: 6111-T4', yield_ddstress_eq(yld(2), s), &
                      expected, tol=TOL)
    ! Isotropic limit with a = 2, pure shear xy: one principal value is exactly zero, where
    ! the curvature has a finite limit
    call yield_cpb06ex2_init_dp(iso, 1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0,1D0,1D0,1D0,0D0, &
                  1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0,1D0,1D0,1D0,0D0,2D0, status=status)
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06ex2_init failed"
        return
    end if
    call s%init((/0D0,0D0,0D0,50D0,0D0,0D0/))
    call check_approx(passed, 'Hessian vs von Mises: isotropic, zero principal value', &
                      yield_ddstress_eq(iso, s), yield_ddstress_eq(vm, s), tol=TOL)
    ! AZ31B maps with a = 2 under pure shear xy: with k, kp /= 0 the zero principal value has
    ! no classical Hessian, and curvature takes the mean of the one-sided values. Reference:
    ! symmetric second differences of Eq. (8) in 90-digit arithmetic, which give that mean
    ! independently of the spectral formula
    call yield_cpb06ex2_init_dp(az, 0.771D0,-0.914D0,-1.312D0,-0.914D0,0.654D0,-1.209D0, &
                 -1.312D0,-1.209D0,-0.293D0, &
                 1.617D0,1D0,1D0,0.295D0, &
                 0.771D0,0.337D0,-0.764D0,0.337D0,0.646D0,-0.401D0,-0.764D0,-0.401D0,-0.022D0, &
                 0.499D0,1D0,1D0,0.900D0,2D0, status=status)
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06ex2_init failed"
        return
    end if
    call expected%init((/1.607067269481341886D-02, 1.104520186736500050D-02, &
                         1.654252201683665822D-02, 0D0, &
                         8.379554812243580838D-03, 8.379554812243580838D-03, &
                         -5.286676272670880566D-03, -5.758525594694119934D-03, &
                         0D0, 0D0, &
                         -6.913000751205860530D-03, -1.078399642214253829D-02, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call check_approx(passed, 'Hessian vs reference: AZ31B maps, a = 2, zero principal value', &
                      yield_ddstress_eq(az, s), expected, tol=TOL)
end subroutine test_CPB06ex2_hessian_analytical
