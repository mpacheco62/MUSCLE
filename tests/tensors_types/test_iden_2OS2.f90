module tests_I2OS
    use iso_fortran_env
    use iso_c_binding
    use muscle_tensors
    implicit none
    contains

integer function test_iden_2OS_sum() result(notPassed) bind(C)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    type(iden_2OS) :: to_test1, to_test2
    type(iden_2OS) :: temp
    real(real64), parameter :: eps=1e-15

    notPassed = 0

    call to_test1%init(0D0)
    temp = to_test1 + to_test1
    if (abs(temp%val) > eps) notPassed = 1
    if (notPassed == 1) return
    
    call to_test2%init(1D0)
    temp = (to_test1 + to_test2)
    if (abs(temp%val - 1D0) > eps) notPassed = 1
    if (notPassed == 1) return

    temp = (to_test2 + to_test1)
    if (abs(temp%val - 1D0) > eps) notPassed = 1
    if (notPassed == 1) return

    call to_test1%init(2D0)
    temp = (to_test1 + to_test2) 
    if (abs(temp%val - 3D0) > eps) notPassed = 1
    if (notPassed  == 1) return

    temp = (to_test2 + to_test1) 
    if(abs(temp%val - 3D0) > eps) notPassed = 1
    if (notPassed  == 1) return

end function

! subroutine test_iden_2OS_sub(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     implicit none
!     logical, intent(out) :: passed
!     type(iden_2OS) :: to_test1, to_test2
!     type(iden_2OS) :: temp
!     real(real64), parameter :: eps=1e-15

!     call to_test1%init(0D0)
!     temp = to_test1 - to_test1
!     passed = abs(temp%val) < eps
!     if (.not. passed) return
    
!     call to_test2%init(1D0)
!     temp = to_test2 - to_test1
!     passed = abs(temp%val - 1D0) < eps
!     if (.not. passed) return

!     temp = to_test1 - to_test2
!     passed = abs(temp%val - (-1D0)) < eps
!     if (.not. passed) return

!     call to_test1%init(3D0)
!     temp = to_test1 - to_test2
!     passed = abs(temp%val - 2D0) < eps
!     if (.not. passed) return

!     temp = -to_test1
!     passed = abs(temp%val - (-3D0)) < eps
!     if (.not. passed) return

! end subroutine

! subroutine test_iden_2OS_mul(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     implicit none
!     logical, intent(out) :: passed
!     type(iden_2OS) :: to_test1
!     type(iden_2OS) :: temp
!     real(real64), parameter :: eps=1e-15

!     call to_test1%init(0D0)
!     temp = 10D0*to_test1
!     passed = abs(temp%val) < eps
!     if (.not. passed) return

!     temp = to_test1*10D0
!     passed = abs(temp%val) < eps
!     if (.not. passed) return

!     call to_test1%init(1D0)
!     temp = to_test1*0D0
!     passed = abs(temp%val) < eps
!     if (.not. passed) return

!     temp = 0D0*to_test1
!     passed = abs(temp%val) < eps
!     if (.not. passed) return

!     temp = 2D0*to_test1
!     passed = abs(temp%val-2D0) < eps
!     if (.not. passed) return

!     temp = to_test1*2D0
!     passed = abs(temp%val-2D0) < eps
!     if (.not. passed) return

! end subroutine

! subroutine test_iden_2OS_div(passed)
!     use, intrinsic :: iso_fortran_env
!     use muscle_tensors
!     implicit none
!     logical, intent(out) :: passed
!     type(iden_2OS) :: to_test1
!     type(iden_2OS) :: temp
!     real(real64), parameter :: eps=1e-15

!     call to_test1%init(0D0)
!     temp = to_test1/10D0
!     passed = abs(temp%val) < eps
!     if (.not. passed) return


!     call to_test1%init(2D0)
!     temp = to_test1/1D0
!     passed = abs(temp%val-2D0) < eps
!     if (.not. passed) return

!     temp = to_test1/2D0
!     passed = abs(temp%val-1D0) < eps
!     if (.not. passed) return
! end subroutine

end module
