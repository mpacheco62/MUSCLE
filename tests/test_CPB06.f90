program test_CPB06
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

    call test_CPB06_stresseq_derivates_5(passed)
    if (.not. passed) STOP 5

    call test_CPB06_stresseq_derivates2_6(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
    STOP 0
end program test_CPB06

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

