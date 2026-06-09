program test_muscle_yield_vonmises
    implicit none
    
    logical :: passed

    call test_vonMises_stresseq_hydrostatic(passed)
    if (.not. passed) STOP 1

    call test_vonMises_stresseq_simple_tensile(passed)
    if (.not. passed) STOP 2

    call test_vonMises_stresseq_zero_stress(passed)
    if (.not. passed) STOP 3

    call test_vonMises_stresseq_biaxial(passed)
    if (.not. passed) STOP 4

    call test_vonMises_stresseq_shear(passed)
    if (.not. passed) STOP 5

    call test_vonMises_stresseq_derivates(passed)
    if (.not. passed) STOP 6

    call test_vonMises_stresseq_derivates2(passed)
    if (.not. passed) STOP 7

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_vonmises

subroutine test_vonMises_stresseq_hydrostatic(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_vonmises
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: to_test1 

    real(real64) :: result
    real(real64) :: expected_result1 = 0.0D0


    call to_test1%init((/5D0, 5D0, 5D0, 0D0, 0D0, 0D0/))

    result = vm%stress_eq(to_test1)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return
end subroutine

subroutine test_vonMises_stresseq_simple_tensile(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_vonmises
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: to_test2

    real(real64) :: result
    real(real64) :: expected_result2 = 5D0

    call to_test2%init((/5D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    
    result = vm%stress_eq(to_test2)
    passed = (abs(result - expected_result2) < EPS)
    if (.not. passed) return
    
end subroutine

subroutine test_vonMises_stresseq_zero_stress(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_vonmises
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: to_test1 

    real(real64) :: result
    real(real64) :: expected_result1 = 0.0D0


    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))

    result = vm%stress_eq(to_test1)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return
end subroutine

subroutine test_vonMises_stresseq_biaxial(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_vonmises
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: to_test3

    real(real64) :: result
    real(real64) :: expected_result3 = 5D0


    call to_test3%init((/5D0, 5D0, 0D0, 0D0, 0D0, 0D0/))
    result = vm%stress_eq(to_test3)
    passed = (abs(result - expected_result3) < EPS)
    if (.not. passed) return

end subroutine

subroutine test_vonMises_stresseq_shear(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_vonmises
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: to_test4

    real(real64) :: result
    real(real64) :: expected_result4 = 3D0**0.5D0 

    call to_test4%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result = vm%stress_eq(to_test4)
    passed = (abs(result - expected_result4) < EPS)
    if (.not. passed) return
end subroutine

subroutine test_vonMises_stresseq_derivates(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_vonmises
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: to_test

    type(ten_3D2Osym) :: result1, result2
    real(real64), parameter :: TOL = 1D-10

    call to_test%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = vm%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = vm%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    result1 = vm%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result1 = vm%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
    result1 = vm%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
    result1 = vm%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = vm%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) return

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = vm%dstressEq_dstress(to_test)
    result2 = vm%dstressEq_dstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=TOL)
    if (.not. passed) return

end subroutine


subroutine test_vonMises_stresseq_derivates2(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_vonmises
    implicit none
    
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: to_test

    type(ten_3D4O3sym) :: result1, result2

    call to_test%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = vm%ddstressEq_ddstress(to_test)
    result2 = vm%ddstressEq_ddstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=1D-7) 
    if (.not. passed) then
        print*, "1.- Error: ddstressEq_ddstress failed"
        print*, "Numerical:"
        print*, result2
        print*, "Analytical:"
        print*, result1
        print*, "differences:"
        print*, result2-result1
        return
    end if

    call to_test%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = vm%ddstressEq_ddstress(to_test)
    result2 = vm%ddstressEq_ddstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=1D-7) 
    if (.not. passed) then
        print*, "2.- Error: ddstressEq_ddstress failed"
        print*, "Numerical:"
        print*, result2
        print*, "Analytical:"
        print*, result1
        print*, "differences:"
        print*, result2-result1
        return
    end if

    call to_test%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    result1 = vm%ddstressEq_ddstress(to_test)
    result2 = vm%ddstressEq_ddstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=1D-7) 
    if (.not. passed) then
        print*, "3.- Error: ddstressEq_ddstress failed"
        print*, "Numerical:"
        print*, result2
        print*, "Analytical:"
        print*, result1
        print*, "differences:"
        print*, result2-result1
        return
    end if

    call to_test%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result1 = vm%ddstressEq_ddstress(to_test)
    result2 = vm%ddstressEq_ddstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=1D-7) 
    if (.not. passed) then
        print*, "4.- Error: ddstressEq_ddstress failed"
        print*, "Numerical:"
        print*, result2
        print*, "Analytical:"
        print*, result1
        print*, "differences:"
        print*, result2-result1
        return
    end if

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
    result1 = vm%ddstressEq_ddstress(to_test)
    result2 = vm%ddstressEq_ddstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=1D-7) 
    if (.not. passed) then
        print*, "5.- Error: ddstressEq_ddstress failed"
        print*, "Numerical:"
        print*, result2
        print*, "Analytical:"
        print*, result1
        print*, "differences:"
        print*, result2-result1
        return
    end if

    call to_test%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
    result1 = vm%ddstressEq_ddstress(to_test)
    result2 = vm%ddstressEq_ddstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=1D-7) 
    if (.not. passed) then
        print*, "6.- Error: ddstressEq_ddstress failed"
        print*, "Numerical:"
        print*, result2
        print*, "Analytical:"
        print*, result1
        print*, "differences:"
        print*, result2-result1
        return
    end if

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = vm%ddstressEq_ddstress(to_test)
    result2 = vm%ddstressEq_ddstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=1D-7) 
    if (.not. passed) then
        print*, "7.- Error: ddstressEq_ddstress failed"
        print*, "Numerical:"
        print*, result2
        print*, "Analytical:"
        print*, result1
        print*, "differences:"
        print*, result2-result1
        return
    end if

    call to_test%init((/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result1 = vm%ddstressEq_ddstress(to_test)
    result2 = vm%ddstressEq_ddstress_numeric(to_test)
    passed = result1%is_approx(result2, tol=1D-7) 
    if (.not. passed) then
        print*, "8.- Error: ddstressEq_ddstress failed"
        print*, "Numerical:"
        print*, result2
        print*, "Analytical:"
        print*, result1
        print*, "differences:"
        print*, result2-result1
        return
    end if

end subroutine

