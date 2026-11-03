program test_muscle_hard_ludwik
    implicit none
    
    logical :: passed

    call test_ludwik_like_hollomon(passed)
    if (.not. passed) STOP 1

    call test_ludwik(passed)
    if (.not. passed) STOP 2

    call test_dludwik(passed)
    if (.not. passed) STOP 3

    call test_ddludwik(passed)
    if (.not. passed) STOP 4

    print*, "Passed!", passed
    STOP 0
end program test_muscle_hard_ludwik

subroutine test_ludwik_like_hollomon(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_ludwik_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_ludwik_dp) :: ludwik

    real(real64) :: result
    real(real64) :: expected_result1

    ludwik = hard_ludwik_dp(sigma0=0D0, k=1D0, n=1D0)

    expected_result1 = 0D0
    result = hard_stress(ludwik, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = hard_stress(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    ludwik = hard_ludwik_dp(sigma0=0D0, k=1D0, n=2D0)

    expected_result1 = 1D0
    result = hard_stress(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 4D0
    result = hard_stress(ludwik, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 9D0
    result = hard_stress(ludwik, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    ludwik = hard_ludwik_dp(sigma0=0D0, k=2D0, n=1D0)

    expected_result1 = 2D0
    result = hard_stress(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 4D0
    result = hard_stress(ludwik, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 6D0
    result = hard_stress(ludwik, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


!***********************************************************************************    
    ludwik = hard_ludwik_dp(sigma0=0D0, k=2D0, n=2D0)

    expected_result1 = 2D0
    result = hard_stress(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 8D0
    result = hard_stress(ludwik, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 18D0
    result = hard_stress(ludwik, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine


subroutine test_ludwik(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_ludwik_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_ludwik_dp) :: ludwik

    real(real64) :: result
    real(real64) :: expected_result1

    ludwik = hard_ludwik_dp(sigma0=1D0, k=1D0, n=1D0)

    expected_result1 = 1D0
    result = hard_stress(ludwik, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = hard_stress(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    ludwik = hard_ludwik_dp(sigma0=1D0, k=1D0, n=2D0)

    expected_result1 = 1D0
    result = hard_stress(ludwik, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = hard_stress(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 5D0
    result = hard_stress(ludwik, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    ludwik = hard_ludwik_dp(sigma0=1D0, k=2D0, n=1D0)

    expected_result1 = 1D0
    result = hard_stress(ludwik, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 3D0
    result = hard_stress(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 5D0
    result = hard_stress(ludwik, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


!***********************************************************************************    
    ludwik = hard_ludwik_dp(sigma0=1D0, k=2D0, n=2D0)

    expected_result1 = 1D0
    result = hard_stress(ludwik, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 3D0
    result = hard_stress(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 9D0
    result = hard_stress(ludwik, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


end subroutine


subroutine test_dludwik(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_ludwik_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_ludwik_dp) :: ludwik

    real(real64) :: result
    real(real64) :: expected_result1

    ludwik = hard_ludwik_dp(sigma0=1D0, k=1D0, n=1D0)

    expected_result1 = 1D0
    result = hard_dstress_dep(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = hard_dstress_dep(ludwik, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    ludwik = hard_ludwik_dp(sigma0=1D0, k=1D0, n=0.5D0)

    expected_result1 = 0.25D0
    result = hard_dstress_dep(ludwik, 4D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0/6D0
    result = hard_dstress_dep(ludwik, 9D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    ludwik = hard_ludwik_dp(sigma0=1D0, k=2D0, n=0.5D0)

    expected_result1 = 0.5D0
    result = hard_dstress_dep(ludwik, 4D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0/3D0
    result = hard_dstress_dep(ludwik, 9D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine


subroutine test_ddludwik(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_ludwik_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_ludwik_dp) :: ludwik

    real(real64) :: result
    real(real64) :: expected_result1

    ludwik = hard_ludwik_dp(sigma0=1D0, k=1D0, n=1D0)

    expected_result1 = 0D0
    result = hard_ddstress_ddep(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 0D0
    result = hard_ddstress_ddep(ludwik, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    ludwik = hard_ludwik_dp(sigma0=1D0, k=1D0, n=0.5D0)

    expected_result1 = -0.25D0
    result = hard_ddstress_ddep(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = -1D0/32D0
    result = hard_ddstress_ddep(ludwik, 4D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    ludwik = hard_ludwik_dp(sigma0=1D0, k=2D0, n=0.5D0)

    expected_result1 = -0.5D0
    result = hard_ddstress_ddep(ludwik, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = -1D0/16D0
    result = hard_ddstress_ddep(ludwik, 4D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return
end subroutine