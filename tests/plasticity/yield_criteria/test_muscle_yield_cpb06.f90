! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
program test_muscle_yield_cpb06
    implicit none
    
    logical :: passed

    call test_CPB06_stresseq_hydrostatic_1(passed)
    if (.not. passed) STOP 1

    call test_CPB06_stresseq_simple_tensile_2(passed)
    if (.not. passed) STOP 2

    call test_CPB06_stresseq_biaxial_3(passed)
    if (.not. passed) STOP 3

    call test_CPB06_stresseq_shear_4(passed)
    if (.not. passed) STOP 4

    call test_CPB06_gradient_numeric_7(passed)
    if (.not. passed) STOP 7

    call test_CPB06_gradient_analytical_8(passed)
    if (.not. passed) STOP 8

    call test_CPB06_hessian_numeric_9(passed)
    if (.not. passed) STOP 9

    call test_CPB06_hessian_analytical_10(passed)
    if (.not. passed) STOP 10

    call test_CPB06_stresseq_derivates_5(passed)
    if (.not. passed) STOP 5

    call test_CPB06_stresseq_derivates2_6(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_cpb06

subroutine test_CPB06_stresseq_hydrostatic_1(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    
    ! Analytical value of the isotropic limit (C = I, k = 0, a = 2: von Mises). Measured: exact;
    ! the relative tolerance 1e-12 (~5000 eps) leaves margin for the eigenvalue round-off
    real(real64), parameter :: TOL = 1D-12
    logical, intent(out) :: passed

    type(yield_cpb06_dp) :: cpb
    type(ten_3D2Osym) :: to_test1 

    real(real64) :: result
    real(real64) :: expected_result1 = 0.0D0

    passed = .true.
    call yield_cpb06_init_dp(cpb, c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0, status=status &
                  )
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06_init failed"
        return
    end if


    call to_test1%init((/5D0, 5D0, 5D0, 0D0, 0D0, 0D0/))

    result = yield_stress_eq(cpb, to_test1)
    call check_approx(passed, "stress_eq, hydrostatic (5,5,5,0,0,0) = 0", result, expected_result1, tol=TOL)
    if (.not. passed) return
end subroutine

subroutine test_CPB06_stresseq_simple_tensile_2(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    
    ! Analytical value of the isotropic limit (C = I, k = 0, a = 2: von Mises). Measured: exact;
    ! the relative tolerance 1e-12 (~5000 eps) leaves margin for the eigenvalue round-off
    real(real64), parameter :: TOL = 1D-12
    logical, intent(out) :: passed

    type(yield_cpb06_dp) :: cpb
    type(ten_3D2Osym) :: to_test2

    real(real64) :: result
    real(real64) :: expected_result2 = 5D0

    passed = .true.
    call yield_cpb06_init_dp(cpb, c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0, status=status &
                  )
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06_init failed"
        return
    end if

    call to_test2%init((/5D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    
    result = yield_stress_eq(cpb, to_test2)
    call check_approx(passed, "stress_eq, uniaxial (5,0,0,0,0,0) = 5", result, expected_result2, tol=TOL)
    if (.not. passed) return
    
end subroutine

subroutine test_CPB06_stresseq_biaxial_3(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    
    ! Analytical value of the isotropic limit (C = I, k = 0, a = 2: von Mises). Measured: exact;
    ! the relative tolerance 1e-12 (~5000 eps) leaves margin for the eigenvalue round-off
    real(real64), parameter :: TOL = 1.0e-12_real64
    logical, intent(out) :: passed

    type(yield_cpb06_dp) :: cpb
    type(ten_3D2Osym) :: to_test3

    real(real64) :: result
    real(real64) :: expected_result3 = 5D0

    passed = .true.
    call yield_cpb06_init_dp(cpb, c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0, status=status &
                  )
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06_init failed"
        return
    end if


    call to_test3%init((/5D0, 5D0, 0D0, 0D0, 0D0, 0D0/))
    result = yield_stress_eq(cpb, to_test3)
    call check_approx(passed, "stress_eq, equibiaxial (5,5,0,0,0,0) = 5", result, expected_result3, tol=TOL)
    if (.not. passed) return

end subroutine

subroutine test_CPB06_stresseq_shear_4(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    
    ! Analytical value of the isotropic limit (C = I, k = 0, a = 2: von Mises). Measured: exact;
    ! the relative tolerance 1e-12 (~5000 eps) leaves margin for the eigenvalue round-off
    real(real64), parameter :: TOL = 1.0e-12_real64
    logical, intent(out) :: passed

    type(yield_cpb06_dp) :: cpb
    type(ten_3D2Osym) :: to_test4

    real(real64) :: result
    real(real64) :: expected_result4 = 3D0**0.5D0 

    passed = .true.
    call yield_cpb06_init_dp(cpb, c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0, status=status &
                  )
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06_init failed"
        return
    end if

    call to_test4%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result = yield_stress_eq(cpb, to_test4)
    call check_approx(passed, "stress_eq, shear (0,0,0,1,0,0) = sqrt(3)", result, expected_result4, tol=TOL)
    if (.not. passed) return
end subroutine

subroutine test_CPB06_stresseq_derivates_5(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    
    ! Gradient against central differences (yield_test_numeric) (step h = 1e-6 ||sigma||, h/2 on
    ! the shears): the rounding error eps*f/(h/2) is ~4e-10 and the truncation O(h**2) ~1e-12.
    ! Measured: 1.3e-10 or less
    real(real64), parameter :: TOL = 1.0e-9_real64
    logical, intent(out) :: passed

    type(yield_cpb06_dp) :: cpb
    type(ten_3D2Osym) :: to_test

    type(ten_3D2Osym) :: result1, result2

    passed = .true.
    call yield_cpb06_init_dp(cpb, c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0, status=status &
                  )
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06_init failed"
        return
    end if

    call to_test%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_dstress_eq(cpb, to_test)
    result2 = num_dstress_eq(cpb, to_test)
    call check_approx(passed, "dstressEq_dstress vs numeric, stress = (1,0,0,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_dstress_eq(cpb, to_test)
    result2 = num_dstress_eq(cpb, to_test)
    call check_approx(passed, "dstressEq_dstress vs numeric, stress = (0,1,0,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    result1 = yield_dstress_eq(cpb, to_test)
    result2 = num_dstress_eq(cpb, to_test)
    call check_approx(passed, "dstressEq_dstress vs numeric, stress = (0,0,1,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result1 = yield_dstress_eq(cpb, to_test)
    result2 = num_dstress_eq(cpb, to_test)
    call check_approx(passed, "dstressEq_dstress vs numeric, stress = (0,0,0,1,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
    result1 = yield_dstress_eq(cpb, to_test)
    result2 = num_dstress_eq(cpb, to_test)
    call check_approx(passed, "dstressEq_dstress vs numeric, stress = (0,0,0,0,1,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
    result1 = yield_dstress_eq(cpb, to_test)
    result2 = num_dstress_eq(cpb, to_test)
    call check_approx(passed, "dstressEq_dstress vs numeric, stress = (0,0,0,0,0,1)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_dstress_eq(cpb, to_test)
    result2 = num_dstress_eq(cpb, to_test)
    call check_approx(passed, "dstressEq_dstress vs numeric, stress = (1,1,0,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_dstress_eq(cpb, to_test)
    result2 = num_dstress_eq(cpb, to_test)
    call check_approx(passed, "dstressEq_dstress vs numeric, stress = (1,1,0,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

end subroutine


subroutine test_CPB06_stresseq_derivates2_6(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    
    ! Hessian against central second differences (yield_test_numeric) (step h = 1e-4 ||sigma||,
    ! h/2 on the shears): the rounding error eps*f/h**2 is ~2e-8 and the truncation O(h**2)
    ! ~1e-8. Measured: 3.5e-8 or less
    real(real64), parameter :: TOL = 1.0e-7_real64
    logical, intent(out) :: passed

    type(yield_cpb06_dp) :: cpb
    type(ten_3D2Osym) :: to_test

    type(ten_3D4O3sym) :: result1, result2

    passed = .true.
    call yield_cpb06_init_dp(cpb, c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0, status=status &
                  )
    if (status /= STATUS_OK) then
        passed = .false.
        print *, "yield_cpb06_init failed"
        return
    end if

    call to_test%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_ddstress_eq(cpb, to_test)
    result2 = num_ddstress_eq(cpb, to_test)
    call check_approx(passed, "ddstressEq_ddstress vs numeric, stress = (1,0,0,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_ddstress_eq(cpb, to_test)
    result2 = num_ddstress_eq(cpb, to_test)
    call check_approx(passed, "ddstressEq_ddstress vs numeric, stress = (0,1,0,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    result1 = yield_ddstress_eq(cpb, to_test)
    result2 = num_ddstress_eq(cpb, to_test)
    call check_approx(passed, "ddstressEq_ddstress vs numeric, stress = (0,0,1,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result1 = yield_ddstress_eq(cpb, to_test)
    result2 = num_ddstress_eq(cpb, to_test)
    call check_approx(passed, "ddstressEq_ddstress vs numeric, stress = (0,0,0,1,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
    result1 = yield_ddstress_eq(cpb, to_test)
    result2 = num_ddstress_eq(cpb, to_test)
    call check_approx(passed, "ddstressEq_ddstress vs numeric, stress = (0,0,0,0,1,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
    result1 = yield_ddstress_eq(cpb, to_test)
    result2 = num_ddstress_eq(cpb, to_test)
    call check_approx(passed, "ddstressEq_ddstress vs numeric, stress = (0,0,0,0,0,1)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_ddstress_eq(cpb, to_test)
    result2 = num_ddstress_eq(cpb, to_test)
    call check_approx(passed, "ddstressEq_ddstress vs numeric, stress = (1,1,0,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_ddstress_eq(cpb, to_test)
    result2 = num_ddstress_eq(cpb, to_test)
    call check_approx(passed, "ddstressEq_ddstress vs numeric, stress = (1,1,0,0,0,0)", result1, result2, tol=TOL)
    if (.not. passed) return

end subroutine

subroutine CPB06_test_materials(cpb)
    ! Materials of the tests of the analytical derivatives
    ! Called by test_CPB06_gradient_numeric_7 and test_CPB06_gradient_analytical_8
    use, intrinsic :: iso_fortran_env
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_status, only : STATUS_OK
    implicit none
    integer :: status
    type(yield_cpb06_dp), intent(out) :: cpb(6)

    ! 1: isotropic limit (C = I, k = 0, a = 2), equal to von Mises
    call yield_cpb06_init_dp(cpb(1), c11=1D0, c12=0D0, c13=0D0, &
                     c21=0D0, c22=1D0, c23=0D0, &
                     c31=0D0, c32=0D0, c33=1D0, &
                     c44=1D0, c55=1D0, c66=1D0, &
                     k=0D0, a=2D0, status=status &
                     )
    if (status /= STATUS_OK) error stop 'yield_cpb06_init failed'
    ! 2: Ti-6Al-4V (Tuninetti et al., 2013, Table 2, row Wp = 1.857)
    call yield_cpb06_init_dp(cpb(2), c11=1.000D0, c12=-2.373D0, c13=-2.364D0, &
                     c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                     c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                     c44=3.607D0, c55=3.607D0, c66=3.607D0,   &
                     k=-0.136D0, a=2D0, status=status         &
                     )
    if (status /= STATUS_OK) error stop 'yield_cpb06_init failed'
    ! 3: material 2 with a non-symmetric C1 (c21 and c32 changed) and three different shear
    !    coefficients: the derivatives need C1^T and the C2 entry of each shear
    call yield_cpb06_init_dp(cpb(3), c11=1.000D0, c12=-2.373D0, c13=-2.364D0, &
                     c21=-1.500D0, c22=-1.838D0, c23=1.196D0, &
                     c31=-2.364D0, c32=0.800D0, c33=-2.444D0, &
                     c44=3.607D0, c55=2.900D0, c66=4.100D0,   &
                     k=-0.136D0, a=2D0, status=status         &
                     )
    if (status /= STATUS_OK) error stop 'yield_cpb06_init failed'
    ! 4: material 2 with a = 8: the powers psi**(a-1) and psi**(a-2) are not trivial
    call yield_cpb06_init_dp(cpb(4), c11=1.000D0, c12=-2.373D0, c13=-2.364D0, &
                     c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                     c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                     c44=3.607D0, c55=3.607D0, c66=3.607D0,   &
                     k=-0.136D0, a=8D0, status=status         &
                     )
    if (status /= STATUS_OK) error stop 'yield_cpb06_init failed'
    ! 5: isotropic with k = 1: psi = |lam| - lam vanishes for the positive principal values
    call yield_cpb06_init_dp(cpb(5), c11=1D0, c12=0D0, c13=0D0, &
                     c21=0D0, c22=1D0, c23=0D0, &
                     c31=0D0, c32=0D0, c33=1D0, &
                     c44=1D0, c55=1D0, c66=1D0, &
                     k=1D0, a=2D0, status=status &
                     )
    if (status /= STATUS_OK) error stop 'yield_cpb06_init failed'
    ! 6: isotropic with k = -1: psi = |lam| + lam vanishes for the negative principal values
    call yield_cpb06_init_dp(cpb(6), c11=1D0, c12=0D0, c13=0D0, &
                     c21=0D0, c22=1D0, c23=0D0, &
                     c31=0D0, c32=0D0, c33=1D0, &
                     c44=1D0, c55=1D0, c66=1D0, &
                     k=-1D0, a=2D0, status=status &
                     )
    if (status /= STATUS_OK) error stop 'yield_cpb06_init failed'
end subroutine CPB06_test_materials

subroutine test_CPB06_gradient_numeric_7(passed)
    ! Gradient against the finite-difference gradient (yield_test_numeric). Measured differences:
    ! 1.3e-10 or less in the general state and 1.2e-7 in pure shear, where the finite differences
    ! straddle the zero principal value
    ! Called by program test_muscle_yield_cpb06 (STOP 7)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0e-6_real64
    type(yield_cpb06_dp) :: cpb(6)
    type(ten_3D2Osym) :: to_test
    type(ten_3D2Osym) :: result1, result2

    passed = .true.
    call CPB06_test_materials(cpb)

    ! General state: every normal and shear component is non-zero, so the shears and C1^T
    ! act on all the components of the result
    call to_test%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))
    result1 = yield_dstress_eq(cpb(2), to_test)
    result2 = num_dstress_eq(cpb(2), to_test)
    call check_approx(passed, "CPB06 gradient vs numeric: Ti-6Al-4V, general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Same state, non-symmetric C1
    result1 = yield_dstress_eq(cpb(3), to_test)
    result2 = num_dstress_eq(cpb(3), to_test)
    call check_approx(passed, "CPB06 gradient vs numeric: non-symmetric C1, general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Same state, a = 8
    result1 = yield_dstress_eq(cpb(4), to_test)
    result2 = num_dstress_eq(cpb(4), to_test)
    call check_approx(passed, "CPB06 gradient vs numeric: a = 8, general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Same state, k = 1: psi = 0 for the positive principal value
    result1 = yield_dstress_eq(cpb(5), to_test)
    result2 = num_dstress_eq(cpb(5), to_test)
    call check_approx(passed, "CPB06 gradient vs numeric: k = 1, general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Same state, k = -1: psi = 0 for the negative principal values
    result1 = yield_dstress_eq(cpb(6), to_test)
    result2 = num_dstress_eq(cpb(6), to_test)
    call check_approx(passed, "CPB06 gradient vs numeric: k = -1, general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Pure shear xy: Sigma has a zero principal value, where h(lam) = (|lam| - k lam)**a
    ! changes branch
    call to_test%init((/0D0, 0D0, 0D0, 200D0, 0D0, 0D0/))
    result1 = yield_dstress_eq(cpb(2), to_test)
    result2 = num_dstress_eq(cpb(2), to_test)
    call check_approx(passed, "CPB06 gradient vs numeric: Ti-6Al-4V, pure shear", result1, result2, tol=TOL)
    if (.not. passed) return
end subroutine test_CPB06_gradient_numeric_7

subroutine test_CPB06_gradient_analytical_8(passed)
    ! Gradient against analytical values: von Mises in the isotropic limit, a high-precision
    ! reference and the exact identities df/dsigma : sigma = f and tr(df/dsigma) = 0
    ! Called by program test_muscle_yield_cpb06 (STOP 8)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    use muscle_yield_vonmises
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0e-12_real64
    type(yield_cpb06_dp) :: cpb(6)
    type(yield_vonmises_dp) :: vm
    type(ten_3D2Osym) :: to_test
    type(ten_3D2Osym) :: result1, result2
    real(real64) :: seq

    passed = .true.
    call CPB06_test_materials(cpb)

    ! Isotropic limit, uniaxial tension along x: Sigma has two equal principal values
    call to_test%init((/300D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_dstress_eq(cpb(1), to_test)
    result2 = yield_dstress_eq(vm, to_test)
    call check_approx(passed, "CPB06 isotropic gradient vs von Mises: uniaxial", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Isotropic limit, general state (normal and shear components together)
    call to_test%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))
    result1 = yield_dstress_eq(cpb(1), to_test)
    result2 = yield_dstress_eq(vm, to_test)
    call check_approx(passed, "CPB06 isotropic gradient vs von Mises: general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Ti-6Al-4V, general state: reference computed with 60-digit arithmetic (mpmath) from an
    ! independent implementation of Eqs. 8, 9 and 12 of Cazacu et al. (2006)
    result1 = yield_dstress_eq(cpb(2), to_test)
    call result2%init((/0.7936671101859936737886589D0, -0.5398406780715101404861498D0, &
                        -0.2538264321144835333025091D0, 0.4019790387089366846503347D0, &
                        -0.117991814026224128346693D0, 0.2714960025436936656127475D0/))
    call check_approx(passed, "CPB06 gradient vs reference: Ti-6Al-4V, general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Non-symmetric C1, general state: homogeneity of degree one and pressure insensitivity
    ! (the gradient is dimensionless, so its trace is compared with TOL directly)
    result1 = yield_dstress_eq(cpb(3), to_test)
    seq = yield_stress_eq(cpb(3), to_test)
    call check_approx(passed, "CPB06 gradient identity df/dsigma : sigma = f: non-symmetric C1", &
                      result1 .ddot. to_test, seq, tol=TOL)
    if (.not. passed) return
    call check_approx(passed, "CPB06 gradient identity tr(df/dsigma) = 0: non-symmetric C1", &
                      result1 .ddot. iden_2O(), 0D0, tol=TOL)
    if (.not. passed) return
end subroutine test_CPB06_gradient_analytical_8

subroutine test_CPB06_hessian_numeric_9(passed)
    ! Hessian against the finite-difference Hessian (yield_test_numeric), in a state with normal
    ! and shear components together (it fills the normal-shear entries)
    ! Called by program test_muscle_yield_cpb06 (STOP 9)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    ! Central second differences (yield_test_numeric) (step h = 1e-4 ||sigma||): rounding of
    ! eps*f/h**2 ~2e-8 relative, and a truncation O(h**2) that grows with a.
    ! Measured: 1.5e-7 or less (a = 8)
    real(real64), parameter :: TOL = 1.0e-6_real64
    type(yield_cpb06_dp) :: cpb(6)
    type(ten_3D2Osym) :: to_test
    type(ten_3D4O3sym) :: result1, result2

    passed = .true.
    call CPB06_test_materials(cpb)
    call to_test%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))

    ! Ti-6Al-4V
    result1 = yield_ddstress_eq(cpb(2), to_test)
    result2 = num_ddstress_eq(cpb(2), to_test)
    call check_approx(passed, "CPB06 Hessian vs numeric: Ti-6Al-4V", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Non-symmetric C1 and different shear coefficients
    result1 = yield_ddstress_eq(cpb(3), to_test)
    result2 = num_ddstress_eq(cpb(3), to_test)
    call check_approx(passed, "CPB06 Hessian vs numeric: non-symmetric C1", result1, result2, tol=TOL)
    if (.not. passed) return

    ! a = 8
    result1 = yield_ddstress_eq(cpb(4), to_test)
    result2 = num_ddstress_eq(cpb(4), to_test)
    call check_approx(passed, "CPB06 Hessian vs numeric: a = 8", result1, result2, tol=TOL)
    if (.not. passed) return

    ! k = 1: the positive principal value has psi = 0, so its curvature is zero
    result1 = yield_ddstress_eq(cpb(5), to_test)
    result2 = num_ddstress_eq(cpb(5), to_test)
    call check_approx(passed, "CPB06 Hessian vs numeric: k = 1", result1, result2, tol=TOL)
    if (.not. passed) return

    ! k = -1: the negative principal values have psi = 0, so their curvature is zero
    result1 = yield_ddstress_eq(cpb(6), to_test)
    result2 = num_ddstress_eq(cpb(6), to_test)
    call check_approx(passed, "CPB06 Hessian vs numeric: k = -1", result1, result2, tol=TOL)
    if (.not. passed) return
end subroutine test_CPB06_hessian_numeric_9

subroutine test_CPB06_hessian_analytical_10(passed)
    ! Hessian against analytical values: von Mises in the isotropic limit and a high-precision
    ! reference
    ! Called by program test_muscle_yield_cpb06 (STOP 10)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_vonmises
    use muscle_yield_criteria
    use yield_test_numeric, only : num_dstress_eq, num_ddstress_eq
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed

    ! Exact references: round-off only. Measured: 5.2e-16; the relative tolerance 1e-12
    ! (~5000 eps) leaves margin
    real(real64), parameter :: TOL = 1.0e-12_real64
    type(yield_cpb06_dp) :: cpb(6)
    type(yield_vonmises_dp) :: vm
    type(ten_3D2Osym) :: to_test
    type(ten_3D4O3sym) :: result1, result2

    passed = .true.
    call CPB06_test_materials(cpb)

    ! Isotropic limit, uniaxial tension along x: two equal principal values, where the
    ! quotient (df/dlam_a - df/dlam_b)/(lam_a - lam_b) takes its limit
    call to_test%init((/300D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = yield_ddstress_eq(cpb(1), to_test)
    result2 = yield_ddstress_eq(vm, to_test)
    call check_approx(passed, "CPB06 isotropic Hessian vs von Mises: uniaxial", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Isotropic limit, general state (normal and shear components together)
    call to_test%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))
    result1 = yield_ddstress_eq(cpb(1), to_test)
    result2 = yield_ddstress_eq(vm, to_test)
    call check_approx(passed, "CPB06 isotropic Hessian vs von Mises: general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Ti-6Al-4V, general state: reference computed with 60-digit arithmetic (mpmath) from an
    ! independent implementation of Eqs. 8, 9 and 12 of Cazacu et al. (2006)
    result1 = yield_ddstress_eq(cpb(2), to_test)
    call result2%init(xxxx=1.5321853636983782D-3, yyyy=1.8448464640656092D-3, &
                      zzzz=2.9111428406260130D-3, xyxy=2.7049699608016333D-3, &
                      yzyz=2.3227131577796801D-3, xzxz=2.9246521198296845D-3, &
                      xxyy=-2.3294449356898719D-4, yyzz=-1.6119019704966220D-3, &
                      zzxy=7.2682495498235987D-4, xyyz=3.0688325999276139D-4, &
                      yzxz=3.7960839816057824D-4, xxzz=-1.2992408701293910D-3, &
                      yyxy=6.2994434641387331D-4, zzyz=-1.1186522622978231D-4, &
                      xyxz=-6.1041972737013786D-4, xxxy=-1.3567693013962332D-3, &
                      yyyz=-2.3196787439809610D-4, zzxz=2.9753317956722756D-5, &
                      xxyz=3.4383310062787841D-4, yyxz=7.3900025148785610D-4, &
                      xxxz=-7.6875356944457885D-4)
    call check_approx(passed, "CPB06 Hessian vs reference: Ti-6Al-4V, general state", result1, result2, tol=TOL)
    if (.not. passed) return

    ! Ti-6Al-4V, pure shear xy: Sigma has a zero principal value, where (|lam| - k lam)**2
    ! has no second derivative and curvature takes the mean of the one-sided values.
    ! Reference: symmetric second differences of f in 65-digit arithmetic (mpmath), which
    ! give that mean, independently of the spectral formula
    call to_test%init((/0D0, 0D0, 0D0, 200D0, 0D0, 0D0/))
    result1 = yield_ddstress_eq(cpb(2), to_test)
    call result2%init((/2.7283349833090195D-03, 2.5146288022842053D-03, 2.5156698490754739D-03, &
                        0.0D0, 2.1558789278938954D-03, 2.1558789278938954D-03, &
                        -1.3636469682588755D-03, -1.1509818340253298D-03, 0.0D0, 0.0D0, &
                        5.7574999645274949D-04, -1.3646880150501440D-03, 0.0D0, 0.0D0, &
                        0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/))
    call check_approx(passed, "CPB06 Hessian vs reference: Ti-6Al-4V, pure shear", result1, result2, tol=TOL)
    if (.not. passed) return
end subroutine test_CPB06_hessian_analytical_10
