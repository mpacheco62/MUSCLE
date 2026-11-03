program test_muscle_hard_swift
    implicit none
    
    logical :: passed

    call test_swift_like_hollomon(passed)
    if (.not. passed) STOP 1

    call test_swift(passed)
    if (.not. passed) STOP 2

    call test_dswift(passed)
    if (.not. passed) STOP 3

    call test_ddswift(passed)
    if (.not. passed) STOP 4

    print*, "Passed!", passed
    STOP 0
end program test_muscle_hard_swift

subroutine test_swift_like_hollomon(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_swift_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_swift_dp) :: swift

    real(real64) :: result
    real(real64) :: expected_result1

    swift = hard_swift_dp(k=1D0, n=1D0, e0=0D0)

    expected_result1 = 0D0
    result = hard_stress(swift, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = hard_stress(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    swift = hard_swift_dp(k=1D0, n=2D0, e0=0D0)

    expected_result1 = 1D0
    result = hard_stress(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 4D0
    result = hard_stress(swift, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 9D0
    result = hard_stress(swift, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    swift = hard_swift_dp(k=2D0, n=1D0, e0=0D0)

    expected_result1 = 2D0
    result = hard_stress(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 4D0
    result = hard_stress(swift, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 6D0
    result = hard_stress(swift, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


!***********************************************************************************    
    swift = hard_swift_dp(k=2D0, n=2D0, e0=0D0)

    expected_result1 = 2D0
    result = hard_stress(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 8D0
    result = hard_stress(swift, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 18D0
    result = hard_stress(swift, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine


subroutine test_swift(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_swift_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_swift_dp) :: swift

    real(real64) :: result
    real(real64) :: expected_result1

    swift = hard_swift_dp(k=1D0, n=1D0, e0=1D0)

    expected_result1 = 1D0
    result = hard_stress(swift, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 2D0
    result = hard_stress(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    swift = hard_swift_dp(k=1D0, n=2D0, e0=1D0)

    expected_result1 = 1D0
    result = hard_stress(swift, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 4D0
    result = hard_stress(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 9D0
    result = hard_stress(swift, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************    
    swift = hard_swift_dp(k=2D0, n=1D0, e0=1D0)

    expected_result1 = 2D0
    result = hard_stress(swift, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 4D0
    result = hard_stress(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 6D0
    result = hard_stress(swift, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


!***********************************************************************************    
    swift = hard_swift_dp(k=2D0, n=2D0, e0=1D0)

    expected_result1 = 2D0
    result = hard_stress(swift, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 8D0
    result = hard_stress(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 18.0D0
    result = hard_stress(swift, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return


!***********************************************************************************    
    ! swift = hard_swift_dp(k=100D0, n=0.1D0, e0=1D-4)
    ! result = hard_stress(swift, 0.34857844D0)  ! 90MPa
    ! print*, "swift", result
    ! passed = .false.
    ! return

end subroutine


subroutine test_dswift(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_swift_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_swift_dp) :: swift

    real(real64) :: result
    real(real64) :: expected_result1

    swift = hard_swift_dp(k=1D0, n=1D0, e0=1D0)

    expected_result1 = 1D0
    result = hard_dstress_dep(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0
    result = hard_dstress_dep(swift, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    swift = hard_swift_dp(k=1D0, n=0.5D0, e0=1D0)

    expected_result1 = 0.25D0
    result = hard_dstress_dep(swift, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0/6D0
    result = hard_dstress_dep(swift, 8D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    swift = hard_swift_dp(k=2D0, n=0.5D0, e0=1D0)

    expected_result1 = 0.5D0
    result = hard_dstress_dep(swift, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 1D0/3D0
    result = hard_dstress_dep(swift, 8D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

end subroutine


subroutine test_ddswift(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_hard_laws, only : hard_stress, hard_dstress_dep, hard_ddstress_ddep, hard_swift_dp
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(hard_swift_dp) :: swift

    real(real64) :: result
    real(real64) :: expected_result1

    swift = hard_swift_dp(k=1D0, n=1D0, e0=1D0)

    expected_result1 = 0D0
    result = hard_ddstress_ddep(swift, 1D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = 0D0
    result = hard_ddstress_ddep(swift, 2D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    swift = hard_swift_dp(k=1D0, n=0.5D0, e0=1D0)

    expected_result1 = -0.25D0
    result = hard_ddstress_ddep(swift, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = -1D0/32D0
    result = hard_ddstress_ddep(swift, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

!***********************************************************************************

    swift = hard_swift_dp(k=2D0, n=0.5D0, e0=1D0)

    expected_result1 = -0.5D0
    result = hard_ddstress_ddep(swift, 0D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return

    expected_result1 = -1D0/16D0
    result = hard_ddstress_ddep(swift, 3D0)
    passed = (abs(result - expected_result1) < EPS)
    if (.not. passed) return
end subroutine