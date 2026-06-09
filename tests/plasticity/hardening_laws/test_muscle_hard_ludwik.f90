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
    use muscle_hard_ludwik
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(ludwik_hardening) :: ludwik

    real(real64) :: result
    real(real64) :: expected_result1

    ludwik = ludwik_hardening(sigma0=0D0, k=1D0, n=1D0)

    expected_result1 = 0D0
    result = ludwik%stress(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = ludwik%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    ludwik = ludwik_hardening(sigma0=0D0, k=1D0, n=2D0)

    expected_result1 = 1D0
    result = ludwik%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 4D0
    result = ludwik%stress(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 9D0
    result = ludwik%stress(3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    ludwik = ludwik_hardening(sigma0=0D0, k=2D0, n=1D0)

    expected_result1 = 2D0
    result = ludwik%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 4D0
    result = ludwik%stress(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 6D0
    result = ludwik%stress(3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


!***********************************************************************************    
    ludwik = ludwik_hardening(sigma0=0D0, k=2D0, n=2D0)

    expected_result1 = 2D0
    result = ludwik%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 8D0
    result = ludwik%stress(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 18D0
    result = ludwik%stress(3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine


subroutine test_ludwik(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_ludwik
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(ludwik_hardening) :: ludwik

    real(real64) :: result
    real(real64) :: expected_result1

    ludwik = ludwik_hardening(sigma0=1D0, k=1D0, n=1D0)

    expected_result1 = 1D0
    result = ludwik%stress(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = ludwik%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    ludwik = ludwik_hardening(sigma0=1D0, k=1D0, n=2D0)

    expected_result1 = 1D0
    result = ludwik%stress(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = ludwik%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 5D0
    result = ludwik%stress(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    ludwik = ludwik_hardening(sigma0=1D0, k=2D0, n=1D0)

    expected_result1 = 1D0
    result = ludwik%stress(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 3D0
    result = ludwik%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 5D0
    result = ludwik%stress(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


!***********************************************************************************    
    ludwik = ludwik_hardening(sigma0=1D0, k=2D0, n=2D0)

    expected_result1 = 1D0
    result = ludwik%stress(0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 3D0
    result = ludwik%stress(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 9D0
    result = ludwik%stress(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


!***********************************************************************************    
    ! swift = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    ! result = swift%stress(0.34857844D0)  ! 90MPa
    ! print*, "swift", result
    ! passed = .false.
    ! return

end subroutine


subroutine test_dludwik(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_ludwik
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(ludwik_hardening) :: ludwik

    real(real64) :: result
    real(real64) :: expected_result1

    ludwik = ludwik_hardening(sigma0=1D0, k=1D0, n=1D0)

    expected_result1 = 1D0
    result = ludwik%dstress_dep(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = ludwik%dstress_dep(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    ludwik = Ludwik_hardening(sigma0=1D0, k=1D0, n=0.5D0)

    expected_result1 = 0.25D0
    result = ludwik%dstress_dep(4D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0/6D0
    result = ludwik%dstress_dep(9D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    ludwik = ludwik_hardening(sigma0=1D0, k=2D0, n=0.5D0)

    expected_result1 = 0.5D0
    result = ludwik%dstress_dep(4D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0/3D0
    result = ludwik%dstress_dep(9D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine


subroutine test_ddludwik(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_ludwik
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(ludwik_hardening) :: ludwik

    real(real64) :: result
    real(real64) :: expected_result1

    ludwik = ludwik_hardening(sigma0=1D0, k=1D0, n=1D0)

    expected_result1 = 0D0
    result = ludwik%ddstress_ddep(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 0D0
    result = ludwik%ddstress_ddep(2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    ludwik = ludwik_hardening(sigma0=1D0, k=1D0, n=0.5D0)

    expected_result1 = -0.25D0
    result = ludwik%ddstress_ddep(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = -1D0/32D0
    result = ludwik%ddstress_ddep(4D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    ludwik = ludwik_hardening(sigma0=1D0, k=2D0, n=0.5D0)

    expected_result1 = -0.5D0
    result = ludwik%ddstress_ddep(1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = -1D0/16D0
    result = ludwik%ddstress_ddep(4D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return
end subroutine