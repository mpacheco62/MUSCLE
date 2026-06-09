program test_muscle_math_operations
    use muscle_tensors
    implicit none
    
    logical :: passed

    call test_ten_3D2Osym_approx(passed)
    if (.not. passed) STOP 1

    call test_ten_3D2Osym_sum(passed)
    if (.not. passed) STOP 2

    call test_ten_3D2Osym_sub(passed)
    if (.not. passed) STOP 3

    call test_ten_3D2Osym_mul(passed)
    if (.not. passed) STOP 4

    call test_ten_3D2Osym_div(passed)
    if (.not. passed) STOP 5

    call test_ten_3D2Osym_dev(passed)
    if (.not. passed) STOP 6

    call test_ten_3D2Osym_ddot(passed)
    if (.not. passed) STOP 7

    call test_ten_3D4O3sym_approx(passed)
    if (.not. passed) STOP 8

    call test_ten_3D4O3sym_sum(passed)
    if (.not. passed) STOP 9
    
    call test_ten_3D4O3sym_sub(passed)
    if (.not. passed) STOP 10

    call test_ten_3D4O3sym_mul(passed)
    if (.not. passed) STOP 11

    call test_ten_3D4O3sym_div(passed)
    if (.not. passed) STOP 12

    call test_ten_3D4O3sym_3D2Osym_ddot(passed)
    if (.not. passed) STOP 13

    call test_I2O_3D2Osym_sum(passed)
    if (.not. passed) STOP 14

    call test_I2O_3D2Osym_sub(passed)
    if (.not. passed) STOP 15

    call test_I2_real64_mul(passed)
    if (.not. passed) STOP 16

    call test_I2_real64_div(passed)
    if (.not. passed) STOP 17

    call test_I2O_3D2Osym_ddot(passed)
    if (.not. passed) STOP 18 

    ! call test_dp_sym_6_4(passed)
    ! if (.not. passed) STOP 5

    STOP 0
end program test_muscle_math_operations

subroutine test_ten_3D2Osym_approx(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2

    call to_test1%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0)
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0)

    passed = to_test1 .approx. to_test2
    if (.not. passed) return

    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6.0001D0)
    passed = .not. (to_test1 .approx. to_test2)
    if (.not. passed) return

    call to_test1%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0)
    call to_test2%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0)
    passed = to_test1 .approx. to_test2
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_sum(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    type(ten_3D2Osym) :: expected

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (to_test1 + to_test1) .approx. to_test1
    if (.not. passed) return
    
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test1 + to_test2) .approx. to_test2
    if (.not. passed) return

    passed = (to_test2 + to_test1) .approx. to_test2
    if (.not. passed) return

    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    call expected%init((/12D0, 14D0, 16D0, 18D0, 20D0, 22D0/))
    passed = (to_test2 + to_test1) .approx. expected
    if (.not. passed) return

    passed = (to_test1 + to_test2) .approx. expected
    if (.not. passed) return

end subroutine

subroutine test_ten_3D2Osym_sub(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    type(ten_3D2Osym) :: expected

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (to_test1 - to_test1) .approx. to_test1
    if (.not. passed) return
    
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test2 - to_test1) .approx. to_test2
    if (.not. passed) return

    passed = (to_test1 - to_test2) .approx. (-to_test2)
    if (.not. passed) return

    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    call expected%init((/10D0, 10D0, 10D0, 10D0, 10D0, 10D0/))
    passed = (to_test1 - to_test2) .approx. expected
    if (.not. passed) return

    passed = (to_test2 - to_test1) .approx. (-expected)
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (10D0*to_test1) .approx. to_test1
    if (.not. passed) return

    passed = (to_test1*10D0) .approx. to_test1
    if (.not. passed) return

    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (0D0*to_test1) .approx. expected
    if (.not. passed) return

    passed = (to_test1*0D0) .approx. expected
    if (.not. passed) return

    call expected%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0/))
    passed = (2D0*to_test1) .approx. expected
    if (.not. passed) return

    passed = (to_test1*2D0) .approx. expected
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (to_test1/2D0) .approx. to_test1
    if (.not. passed) return

    call to_test1%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0/))
    call expected%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    passed = (to_test1/2D0) .approx. expected
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_dev(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (.dev. to_test1) .approx. to_test1
    if (.not. passed) return

    call to_test1%init((/5D0, 5D0, 5D0, 0D0, 0D0, 0D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    passed = (.dev. to_test1) .approx. expected
    if (.not. passed) return

    call to_test1%init((/5D0, 5D0, 5D0, 1D0, 1D0, 1D0/))
    call expected%init((/0D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    passed = (.dev. to_test1) .approx. expected
    if (.not. passed) return

    call to_test1%init((/3D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    call expected%init((/2D0, -1D0, -1D0, 1D0, 1D0, 1D0/))
    passed = (.dev. to_test1) .approx. expected
    if (.not. passed) return
end subroutine

subroutine test_ten_3D2Osym_ddot(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    
    passed = .false.

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    if ((to_test1 .ddot. to_test1) < 1D-10) passed = .true.
    if (.not. passed) return

    passed = .false.
    if ((to_test1 .ddot. to_test2) < 1D-10) passed = .true.
    if (.not. passed) return

    passed = .false.
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    if (abs((to_test1 .ddot. to_test2) - 528D0) < 1D-10) passed = .true.
    if (.not. passed) return

    passed = .false.
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    if (abs((to_test2 .ddot. to_test1) - 528D0) < 1D-10) passed = .true.
    print*, to_test1 .ddot. to_test2
    if (.not. passed) return
end subroutine

subroutine test_ten_3D4O3sym_approx(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D4O3sym) :: to_test1, to_test2

    call to_test1%init((/ 1D0,  2D0,  3D0,  4D0,  5D0,  6D0,  7D0,  8D0,  9D0, 10D0, 11D0,&
                         12D0, 13D0, 14D0, 15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0 /))
    call to_test2%init((/ 1D0,  2D0,  3D0,  4D0,  5D0,  6D0,  7D0,  8D0,  9D0, 10D0, 11D0,&
                         12D0, 13D0, 14D0, 15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0 /))

    passed = to_test1 .approx. to_test2
    if (.not. passed) return

    call to_test2%init((/ 1D0,  2D0,  3D0,  4D0,  5D0,  6D0,  7D0,  8D0,  9D0, 10D0, 11D0,&
                         12D0, 13D0, 14D0, 15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 0D0 /))

    passed = .not. (to_test1 .approx. to_test2)
    if (.not. passed) return
end subroutine

subroutine test_ten_3D4O3sym_sum(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D4O3sym) :: to_test1, to_test2
    type(ten_3D4O3sym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  &
                         /))
    passed = (to_test1 + to_test1) .approx. to_test1
    if (.not. passed) return
    
    call to_test2%init((/ 1D0,  2D0,  3D0,  4D0,  5D0,  6D0,  7D0, &
                          8D0,  9D0, 10D0, 11D0, 12D0, 13D0, 14D0, &
                         15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0  &
                        /))
    passed = (to_test1 + to_test2) .approx. to_test2
    if (.not. passed) return

    passed = (to_test2 + to_test1) .approx. to_test2
    if (.not. passed) return

    call to_test1%init((/101D0, 102D0, 103D0, 104D0, 105D0, 106D0, 107D0, &
                         108D0, 109D0, 110D0, 111D0, 112D0, 113D0, 114D0, &
                         115D0, 116D0, 117D0, 118D0, 119D0, 120D0, 121D0  &
                        /))

    call expected_result%init((/102D0, 104D0, 106D0, 108D0, 110D0, 112D0, 114D0, &
                                116D0, 118D0, 120D0, 122D0, 124D0, 126D0, 128D0, &
                                130D0, 132D0, 134D0, 136D0, 138D0, 140D0, 142D0  &
                                /))
    passed = (to_test2 + to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1 + to_test2) .approx. expected_result
    if (.not. passed) return

end subroutine

subroutine test_ten_3D4O3sym_sub(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D4O3sym) :: to_test1, to_test2
    type(ten_3D4O3sym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  &
                         /))
    passed = (to_test1 - to_test1) .approx. to_test1
    if (.not. passed) return
    
    call to_test2%init((/ 1D0,  2D0,  3D0,  4D0,  5D0,  6D0,  7D0, &
                          8D0,  9D0, 10D0, 11D0, 12D0, 13D0, 14D0, &
                         15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0  &
                        /))
    passed = (to_test2 - to_test1) .approxto_test2
    if (.not. passed) return

    passed = (to_test1 - to_test2) .approx (-to_test2)
    if (.not. passed) return

    call to_test1%init((/101D0, 102D0, 103D0, 104D0, 105D0, 106D0, 107D0, &
                         108D0, 109D0, 110D0, 111D0, 112D0, 113D0, 114D0, &
                         115D0, 116D0, 117D0, 118D0, 119D0, 120D0, 121D0  &
                        /))

    call expected_result%init((/100D0, 100D0, 100D0, 100D0, 100D0, 100D0, 100D0, &
                                100D0, 100D0, 100D0, 100D0, 100D0, 100D0, 100D0, &
                                100D0, 100D0, 100D0, 100D0, 100D0, 100D0, 100D0 &
                                /))
    passed = (to_test1 - to_test2) .approx expected_result
    if (.not. passed) return

    passed = (to_test2 - to_test1) .approx. (-expected_result)
    if (.not. passed) return
end subroutine

subroutine test_ten_3D4O3sym_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D4O3sym) :: to_test1
    type(ten_3D4O3sym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  &
                         /))
    passed = (10D0*to_test1) .approx. to_test1
    if (.not. passed) return

    passed = (to_test1*10D0) .approx. to_test1
    if (.not. passed) return

    call to_test1%init((/ 1D0,  2D0,  3D0,  4D0,  5D0,  6D0,  7D0, &
                          8D0,  9D0, 10D0, 11D0, 12D0, 13D0, 14D0, &
                         15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0  &
                         /))
    call expected_result%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                                0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                                0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  &
                                /))
    passed = (0D0*to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1*0D0) .approx. expected_result
    if (.not. passed) return

    call expected_result%init((/ 2D0,  4D0,  6D0,  8D0, 10D0, 12D0, 14D0, &
                                16D0, 18D0, 20D0, 22D0, 24D0, 26D0, 28D0, &
                                30D0, 32D0, 34D0, 36D0, 38D0, 40D0, 42D0  &
                                /))
    passed = (2D0*to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1*2D0) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_ten_3D4O3sym_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D4O3sym) :: to_test1
    type(ten_3D4O3sym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  &
                         /))
    passed = (to_test1/2D0) .approx. to_test1
    if (.not. passed) return

    call to_test1%init((/ 2D0,  4D0,  6D0,  8D0, 10D0, 12D0, 14D0, &
                         16D0, 18D0, 20D0, 22D0, 24D0, 26D0, 28D0, &
                         30D0, 32D0, 34D0, 36D0, 38D0, 40D0, 42D0  &
                         /))
    call expected_result%init((/ 1D0,  2D0,  3D0,  4D0,  5D0,  6D0,  7D0, &
                                 8D0,  9D0, 10D0, 11D0, 12D0, 13D0, 14D0, &
                                15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0  &
                                /))
    passed = (to_test1/2D0) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_ten_3D4O3sym_3D2Osym_ddot(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D4O3sym) :: to_testO4
    type(ten_3D2Osym) :: to_testO2, expected_result
    
    call to_testO4%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                          0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                          0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0  &
                          /))
    call to_testO2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected_result%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))

    passed = (to_testO4 .ddot. to_testO2) .approx. expected_result
    if (.not. passed) return

    passed = (to_testO2 .ddot. to_testO4) .approx. expected_result
    if (.not. passed) return

    call to_testO4%init((/ 1D0,  2D0,  3D0,  4D0,  5D0,  6D0,  7D0, &
                           8D0,  9D0, 10D0, 11D0, 12D0, 13D0, 14D0, &
                          15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0  &
                          /))
    call expected_result%init((/621D0, 549D0, 465D0, 381D0, 357D0, 417D0/))

    passed = (to_testO2 .ddot. to_testO4) .approx. expected_result
    if (.not. passed) return

    passed = (to_testO4 .ddot. to_testO2) .approx. expected_result
    if (.not. passed) return

end subroutine

subroutine test_I2O_3D2Osym_sum(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected_result%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    passed = (to_test1 + iden_2O()) .approx. expected_result
    if (.not. passed) return
    
    passed = (iden_2O() + to_test1) .approx. expected_result
    if (.not. passed) return

    call to_test1%init((/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected_result%init((/2D0, 2D0, 2D0, 1D0, 1D0, 1D0/))
    passed = (to_test1 + iden_2O()) .approx. expected_result
    if (.not. passed) return
    
    passed = (iden_2O() + to_test1) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_I2O_3D2Osym_sub(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected_result

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected_result%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    passed = (iden_2O() - to_test1) .approx. expected_result
    if (.not. passed) return
    
    passed = (to_test1 - iden_2O()) .approx. (-expected_result)
    if (.not. passed) return

    call to_test1%init((/2D0, 2D0, 2D0, 1D0, 1D0, 1D0/))
    call expected_result%init((/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    passed = (to_test1 - iden_2O()) .approx. expected_result
    if (.not. passed) return
    
    passed = (iden_2O() - to_test1) .approx. (-expected_result)
    if (.not. passed) return
end subroutine

subroutine test_I2_real64_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: expected_result

    call expected_result%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    passed = (1D0*iden_2O()) .approx. expected_result
    if (.not. passed) return
    
    passed = (iden_2O()*1D0) .approx. expected_result
    if (.not. passed) return

    passed = ((-1D0)*iden_2O()) .approx. (-expected_result)
    if (.not. passed) return
    
    passed = (iden_2O()*(-1D0)) .approx. (-expected_result)
    if (.not. passed) return
    
    call expected_result%init((/2D0, 2D0, 2D0, 0D0, 0D0, 0D0/))
    passed = (2D0*iden_2O()) .approx. expected_result
    if (.not. passed) return
    
    passed = (iden_2O()*2D0) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_I2_real64_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: expected_result

    call expected_result%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    passed = (iden_2O()/1D0) .approx. expected_result
    if (.not. passed) return
    
    passed = (iden_2O()/(-1D0)) .approx. (-expected_result)
    if (.not. passed) return
    
    call expected_result%init((/2D0, 2D0, 2D0, 0D0, 0D0, 0D0/))
    passed = (iden_2O()/(0.5D0)) .approx. expected_result
    if (.not. passed) return
    
end subroutine

subroutine test_I2O_3D2Osym_ddot(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed
    type(ten_3D2Osym) :: to_test1
    real(real64), parameter :: EPS=1D-7

    passed = .false.

    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    if (abs(iden_2O() .ddot. to_test1) > EPS ) return
    if (abs(to_test1 .ddot. iden_2O()) > EPS ) return

    call to_test1%init((/0D0, 0D0, 0D0, 10D0, 10D0, 10D0/))
    if (abs(iden_2O() .ddot. to_test1) > EPS ) return
    if (abs(to_test1 .ddot. iden_2O()) > EPS ) return
    
    call to_test1%init((/1D0, 1D0, 1D0, 10D0, 10D0, 10D0/))
    if (abs((iden_2O() .ddot. to_test1) -3.0D0) > EPS ) return
    if (abs((to_test1 .ddot. iden_2O()) -3.0D0) > EPS ) return

    passed = .true.
    
end subroutine