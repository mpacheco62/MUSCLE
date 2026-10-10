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
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: to_test1 

    real(real64) :: result
    real(real64) :: expected_result1 = 0.0D0

    call cpb%init(c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0               &
                  )


    call to_test1%init((/5D0, 5D0, 5D0, 0D0, 0D0, 0D0/))

    result = cpb%stress_eq(to_test1)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) print*, "CPB06 hydrostatic not pass, result =", result, "expected =", expected_result1
    if (.not. passed) return
end subroutine

subroutine test_CPB06_stresseq_simple_tensile_2(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: to_test2

    real(real64) :: result
    real(real64) :: expected_result2 = 5D0

    call cpb%init(c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0               &
                  )

    call to_test2%init((/5D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    
    result = cpb%stress_eq(to_test2)
    passed = (abs(result - expected_result2) < EPS)
    if (.not. passed) return
    
end subroutine

subroutine test_CPB06_stresseq_biaxial_3(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: to_test3

    real(real64) :: result
    real(real64) :: expected_result3 = 5D0

    call cpb%init(c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0               &
                  )


    call to_test3%init((/5D0, 5D0, 0D0, 0D0, 0D0, 0D0/))
    result = cpb%stress_eq(to_test3)
    passed = (abs(result - expected_result3) < EPS)
    if (.not. passed) return

end subroutine

subroutine test_CPB06_stresseq_shear_4(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: to_test4

    real(real64) :: result
    real(real64) :: expected_result4 = 3D0**0.5D0 

    call cpb%init(c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0               &
                  )

    call to_test4%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result = cpb%stress_eq(to_test4)
    passed = (abs(result - expected_result4) < EPS)
    if (.not. passed) print*, "CPB06 shear not pass, result =", result, "expected =", expected_result4
    if (.not. passed) return
end subroutine

subroutine test_CPB06_stresseq_derivates_5(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: to_test

    type(ten_3D2Osym) :: result1, result2

    call cpb%init(c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0               &
                  )

    call to_test%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb%dstressEq_dstress(to_test)
    result2 = cpb%dstressEq_dstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb%dstressEq_dstress(to_test)
    result2 = cpb%dstressEq_dstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    result1 = cpb%dstressEq_dstress(to_test)
    result2 = cpb%dstressEq_dstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result1 = cpb%dstressEq_dstress(to_test)
    result2 = cpb%dstressEq_dstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
    result1 = cpb%dstressEq_dstress(to_test)
    result2 = cpb%dstressEq_dstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
    result1 = cpb%dstressEq_dstress(to_test)
    result2 = cpb%dstressEq_dstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb%dstressEq_dstress(to_test)
    result2 = cpb%dstressEq_dstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb%dstressEq_dstress(to_test)
    result2 = cpb%dstressEq_dstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

end subroutine


subroutine test_CPB06_stresseq_derivates2_6(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(CPB06) :: cpb
    type(ten_3D2Osym) :: to_test

    type(ten_3D4O3sym) :: result1, result2

    call cpb%init(c11=1D0, c12=0D0, c13=0D0, &
                  c21=0D0, c22=1D0, c23=0D0, &
                  c31=0D0, c32=0D0, c33=1D0, &
                  c44=1D0, c55=1D0, c66=1D0, &
                  k=0D0, a=2D0               &
                  )

    call to_test%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb%ddstressEq_ddstress(to_test)
    result2 = cpb%ddstressEq_ddstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb%ddstressEq_ddstress(to_test)
    result2 = cpb%ddstressEq_ddstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    result1 = cpb%ddstressEq_ddstress(to_test)
    result2 = cpb%ddstressEq_ddstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result1 = cpb%ddstressEq_ddstress(to_test)
    result2 = cpb%ddstressEq_ddstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
    result1 = cpb%ddstressEq_ddstress(to_test)
    result2 = cpb%ddstressEq_ddstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
    result1 = cpb%ddstressEq_ddstress(to_test)
    result2 = cpb%ddstressEq_ddstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb%ddstressEq_ddstress(to_test)
    result2 = cpb%ddstressEq_ddstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb%ddstressEq_ddstress(to_test)
    result2 = cpb%ddstressEq_ddstress_numeric(to_test)
    passed = result1 .approx. result2
    if (.not. passed) return

end subroutine

subroutine CPB06_test_materials(cpb)
    ! Materials of the tests of the analytical derivatives
    use, intrinsic :: iso_fortran_env
    use muscle_yield_cpb06
    implicit none
    type(CPB06), intent(out) :: cpb(6)

    ! 1: isotropic limit (C = I, k = 0, a = 2), equal to von Mises
    call cpb(1)%init(c11=1D0, c12=0D0, c13=0D0, &
                     c21=0D0, c22=1D0, c23=0D0, &
                     c31=0D0, c32=0D0, c33=1D0, &
                     c44=1D0, c55=1D0, c66=1D0, &
                     k=0D0, a=2D0               &
                     )
    ! 2: Ti-6Al-4V (Tuninetti et al., 2013, Table 2, row Wp = 1.857)
    call cpb(2)%init(c11=1.000D0, c12=-2.373D0, c13=-2.364D0, &
                     c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                     c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                     c44=3.607D0, c55=3.607D0, c66=3.607D0,   &
                     k=-0.136D0, a=2D0                        &
                     )
    ! 3: material 2 with a non-symmetric C1 (c21 and c32 changed) and three different shear
    !    coefficients: the derivatives need C1^T and the C2 entry of each shear
    call cpb(3)%init(c11=1.000D0, c12=-2.373D0, c13=-2.364D0, &
                     c21=-1.500D0, c22=-1.838D0, c23=1.196D0, &
                     c31=-2.364D0, c32=0.800D0, c33=-2.444D0, &
                     c44=3.607D0, c55=2.900D0, c66=4.100D0,   &
                     k=-0.136D0, a=2D0                        &
                     )
    ! 4: material 2 with a = 8: the powers psi**(a-1) and psi**(a-2) are not trivial
    call cpb(4)%init(c11=1.000D0, c12=-2.373D0, c13=-2.364D0, &
                     c21=-2.373D0, c22=-1.838D0, c23=1.196D0, &
                     c31=-2.364D0, c32=1.196D0, c33=-2.444D0, &
                     c44=3.607D0, c55=3.607D0, c66=3.607D0,   &
                     k=-0.136D0, a=8D0                        &
                     )
    ! 5: isotropic with k = 1: psi = |lam| - lam vanishes for the positive principal values
    call cpb(5)%init(c11=1D0, c12=0D0, c13=0D0, &
                     c21=0D0, c22=1D0, c23=0D0, &
                     c31=0D0, c32=0D0, c33=1D0, &
                     c44=1D0, c55=1D0, c66=1D0, &
                     k=1D0, a=2D0               &
                     )
    ! 6: isotropic with k = -1: psi = |lam| + lam vanishes for the negative principal values
    call cpb(6)%init(c11=1D0, c12=0D0, c13=0D0, &
                     c21=0D0, c22=1D0, c23=0D0, &
                     c31=0D0, c32=0D0, c33=1D0, &
                     c44=1D0, c55=1D0, c66=1D0, &
                     k=-1D0, a=2D0              &
                     )
end subroutine CPB06_test_materials

subroutine test_CPB06_gradient_numeric_7(passed)
    ! Gradient against the finite-difference gradient of the base class. Measured differences:
    ! 1.3e-10 or less in the general state and 1.2e-7 in pure shear, where the finite differences
    ! straddle the zero principal value
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0e-6_real64
    type(CPB06) :: cpb(6)
    type(ten_3D2Osym) :: to_test
    type(ten_3D2Osym) :: result1, result2

    call CPB06_test_materials(cpb)

    ! General state: every normal and shear component is non-zero, so the shears and C1^T
    ! act on all the components of the result
    call to_test%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))
    result1 = cpb(2)%dstressEq_dstress(to_test)
    result2 = cpb(2)%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 gradient vs numeric failed: Ti-6Al-4V, general state"
    if (.not. passed) return

    ! Same state, non-symmetric C1
    result1 = cpb(3)%dstressEq_dstress(to_test)
    result2 = cpb(3)%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 gradient vs numeric failed: non-symmetric C1, general state"
    if (.not. passed) return

    ! Same state, a = 8
    result1 = cpb(4)%dstressEq_dstress(to_test)
    result2 = cpb(4)%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 gradient vs numeric failed: a = 8, general state"
    if (.not. passed) return

    ! Same state, k = 1: psi = 0 for the positive principal value
    result1 = cpb(5)%dstressEq_dstress(to_test)
    result2 = cpb(5)%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 gradient vs numeric failed: k = 1, general state"
    if (.not. passed) return

    ! Same state, k = -1: psi = 0 for the negative principal values
    result1 = cpb(6)%dstressEq_dstress(to_test)
    result2 = cpb(6)%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 gradient vs numeric failed: k = -1, general state"
    if (.not. passed) return

    ! Pure shear xy: Sigma has a zero principal value, where h(lam) = (|lam| - k lam)**a
    ! changes branch
    call to_test%init((/0D0, 0D0, 0D0, 200D0, 0D0, 0D0/))
    result1 = cpb(2)%dstressEq_dstress(to_test)
    result2 = cpb(2)%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 gradient vs numeric failed: Ti-6Al-4V, pure shear"
    if (.not. passed) return
end subroutine test_CPB06_gradient_numeric_7

subroutine test_CPB06_gradient_analytical_8(passed)
    ! Gradient against analytical values: von Mises in the isotropic limit, a high-precision
    ! reference and the exact identities df/dsigma : sigma = f and tr(df/dsigma) = 0
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_cpb06
    use muscle_yield_vonmises
    implicit none
    logical, intent(out) :: passed

    real(real64), parameter :: TOL = 1.0e-12_real64
    type(CPB06) :: cpb(6)
    type(VonMises) :: vm
    type(ten_3D2Osym) :: to_test
    type(ten_3D2Osym) :: result1, result2
    real(real64) :: seq

    call CPB06_test_materials(cpb)

    ! Isotropic limit, uniaxial tension along x: Sigma has two equal principal values
    call to_test%init((/300D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = cpb(1)%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 isotropic gradient vs von Mises failed: uniaxial"
    if (.not. passed) return

    ! Isotropic limit, general state (normal and shear components together)
    call to_test%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))
    result1 = cpb(1)%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 isotropic gradient vs von Mises failed: general state"
    if (.not. passed) return

    ! Ti-6Al-4V, general state: reference computed with 60-digit arithmetic (mpmath) from an
    ! independent implementation of Eqs. 8, 9 and 12 of Cazacu et al. (2006)
    result1 = cpb(2)%dstressEq_dstress(to_test)
    call result2%init((/0.7936671101859936737886589D0, -0.5398406780715101404861498D0, &
                        -0.2538264321144835333025091D0, 0.4019790387089366846503347D0, &
                        -0.117991814026224128346693D0, 0.2714960025436936656127475D0/))
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) print*, "CPB06 gradient vs reference failed: Ti-6Al-4V, general state"
    if (.not. passed) return

    ! Non-symmetric C1, general state: homogeneity of degree one and pressure insensitivity
    ! (the gradient is dimensionless, so its trace is compared with TOL directly)
    result1 = cpb(3)%dstressEq_dstress(to_test)
    seq = cpb(3)%stress_eq(to_test)
    passed = abs((result1 .ddot. to_test) - seq) <= TOL*seq &
             .and. abs(result1 .ddot. iden_2O()) <= TOL
    if (.not. passed) print*, "CPB06 gradient identities failed: non-symmetric C1, general state"
    if (.not. passed) return
end subroutine test_CPB06_gradient_analytical_8
