program test_basic_operations
    use basic_operations
    implicit none
    
    logical :: passed

    ! call test_deviatoric_3x3_1(passed)
    ! if (.not. passed) STOP 1

    ! call test_deviatoric_6_2(passed)
    ! if (.not. passed) STOP 2

    ! call test_dp_sym_3x3_3(passed)
    ! if (.not. passed) STOP 3

    ! call test_dp_sym_6_4(passed)
    ! if (.not. passed) STOP 4

    print*, "Hola!", passed
    STOP 0
end program test_basic_operations

! subroutine test_deviatoric_3x3_1(passed)
!     use, intrinsic :: iso_fortran_env
!     use basic_operations
!     implicit none
    
!     logical, intent(out) :: passed

!     real(real64), parameter :: EPS=1e-10
!     real(real64), dimension(3,3) :: result

!     real(real64), dimension(3,3) :: to_test1 = reshape((/1,0,0, &
!                                                          0,1,0, &
!                                                          0,0,1/), shape(to_test1))

!     real(real64), dimension(3,3) :: expected_result1 = reshape((/0,0,0, &
!                                                                  0,0,0, &
!                                                                  0,0,0/), shape(expected_result1))


!     result = deviatoric(to_test1)
!     passed = compare_mat(result, expected_result1, EPS)
!     if (.not. passed) return
    
!     passed = .not. compare_mat(result, expected_result1+1, EPS)
!     if (.not. passed) return
! end subroutine

! subroutine test_deviatoric_6_2(passed)
!     use, intrinsic :: iso_fortran_env
!     use basic_operations
!     implicit none
    
!     logical, intent(out) :: passed

!     real(real64), parameter :: EPS=1e-10
!     real(real64), dimension(6) :: result

!     real(real64), dimension(6) :: to_test1 = (/1,1,1,0,0,0/)
!     real(real64), dimension(6) :: expected_result1 = (/0,0,0,0,0,0/)
                                                          

!     result = deviatoric(to_test1)
!     passed = compare_vec(result, expected_result1, EPS)
!     if (.not. passed) return
    
!     passed = .not. compare_vec(result, expected_result1+1, EPS)
!     if (.not. passed) return
! end subroutine

! subroutine test_dp_sym_3x3_3(passed)
!     use, intrinsic :: iso_fortran_env
!     use basic_operations
!     implicit none
    
!     logical, intent(out) :: passed

!     real(real64), parameter :: EPS=1e-10

!     real(real64), dimension(3,3) :: to_test1 = reshape((/1,2,3, &
!                                                          2,4,5, &
!                                                          3,5,6/), shape(to_test1))
!     real(real64), dimension(3,3) :: to_test2 = reshape((/11,12,13, &
!                                                          12,14,15, &
!                                                          13,15,16/), shape(to_test2))

!     real(real64) :: result
!     real(real64) :: expected_result = 439


!     result = dc_sym(to_test1, to_test2)
!     passed = (abs(result - expected_result) < EPS)
!     if (.not. passed) return
! end subroutine

! subroutine test_dp_sym_6_4(passed)
!     use, intrinsic :: iso_fortran_env
!     use basic_operations
!     implicit none
    
!     logical, intent(out) :: passed

!     real(real64), parameter :: EPS=1e-10

!     real(real64), dimension(6) :: to_test1 = (/1,4,6,2,3,5/)
!     real(real64), dimension(6) :: to_test2 = (/11,14,16,12,13,15/)

!     real(real64) :: result
!     real(real64) :: expected_result = 439


!     result = dc_sym(to_test1, to_test2)
!     passed = (abs(result - expected_result) < EPS)
!     if (.not. passed) return
! end subroutine