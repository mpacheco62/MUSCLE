program test_I2OS
    use tensors_types
    implicit none
    logical :: passed

    ! real*8 :: a
    ! call test123(a, 5D0)

    call test_iden_2OS_sum(passed)
    if (.not. passed) STOP 1

    call test_iden_2OS_sub(passed)
    if (.not. passed) STOP 2

    call test_iden_2OS_mul(passed)
    if (.not. passed) STOP 3

    call test_iden_2OS_div(passed)
    if (.not. passed) STOP 4

    ! call test_ten_3D2O_mul(passed)
    ! if (.not. passed) STOP 4

    ! call test_ten_3D2O_div(passed)
    ! if (.not. passed) STOP 5

    ! call test_ten_3D2O_dev(passed)
    ! if (.not. passed) STOP 6

    ! call test_ten_3D2O_ddot(passed)
    ! if (.not. passed) STOP 7

    print*, "Hola!", passed
    STOP 0
end program test_I2OS


subroutine test_iden_2OS_sum(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    type(iden_2OS) :: to_test1, to_test2
    type(iden_2OS) :: temp
    real(real64), parameter :: eps=1e-15

    call to_test1%init(0D0)
    temp = to_test1 + to_test1
    passed = abs(temp%val) < eps
    if (.not. passed) return
    
    call to_test2%init(1D0)
    temp = (to_test1 + to_test2)
    passed = abs(temp%val - 1D0) < eps
    if (.not. passed) return

    temp = (to_test2 + to_test1)
    passed = abs(temp%val - 1D0) < eps
    if (.not. passed) return

    call to_test1%init(2D0)
    temp = (to_test1 + to_test2) 
    passed = abs(temp%val - 3D0) < eps
    if (.not. passed) return

    temp = (to_test2 + to_test1) 
    passed = abs(temp%val - 3D0) < eps
    if (.not. passed) return

end subroutine

subroutine test_iden_2OS_sub(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    type(iden_2OS) :: to_test1, to_test2
    type(iden_2OS) :: temp
    real(real64), parameter :: eps=1e-15

    call to_test1%init(0D0)
    temp = to_test1 - to_test1
    passed = abs(temp%val) < eps
    if (.not. passed) return
    
    call to_test2%init(1D0)
    temp = to_test2 - to_test1
    passed = abs(temp%val - 1D0) < eps
    if (.not. passed) return

    temp = to_test1 - to_test2
    passed = abs(temp%val - (-1D0)) < eps
    if (.not. passed) return

    call to_test1%init(3D0)
    temp = to_test1 - to_test2
    passed = abs(temp%val - 2D0) < eps
    if (.not. passed) return

    temp = -to_test1
    passed = abs(temp%val - (-3D0)) < eps
    if (.not. passed) return

end subroutine

subroutine test_iden_2OS_mul(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    type(iden_2OS) :: to_test1
    type(iden_2OS) :: temp
    real(real64), parameter :: eps=1e-15

    call to_test1%init(0D0)
    temp = 10D0*to_test1
    passed = abs(temp%val) < eps
    if (.not. passed) return

    temp = to_test1*10D0
    passed = abs(temp%val) < eps
    if (.not. passed) return

    call to_test1%init(1D0)
    temp = to_test1*0D0
    passed = abs(temp%val) < eps
    if (.not. passed) return

    temp = 0D0*to_test1
    passed = abs(temp%val) < eps
    if (.not. passed) return

    temp = 2D0*to_test1
    passed = abs(temp%val-2D0) < eps
    if (.not. passed) return

    temp = to_test1*2D0
    passed = abs(temp%val-2D0) < eps
    if (.not. passed) return

end subroutine

subroutine test_iden_2OS_div(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    type(iden_2OS) :: to_test1
    type(iden_2OS) :: temp
    real(real64), parameter :: eps=1e-15

    call to_test1%init(0D0)
    temp = to_test1/10D0
    passed = abs(temp%val) < eps
    if (.not. passed) return


    call to_test1%init(2D0)
    temp = to_test1/1D0
    passed = abs(temp%val-2D0) < eps
    if (.not. passed) return

    temp = to_test1/2D0
    passed = abs(temp%val-1D0) < eps
    if (.not. passed) return
end subroutine
