program test_2D4O3sym
    use tensors_types
    implicit none
    
    logical :: passed

    call test_ten_2D4O3sym_approx(passed)
    if (.not. passed) STOP 1

    call test_ten_2D4O3sym_sum(passed)
    if (.not. passed) STOP 2

    call test_ten_2D4O3sym_sub(passed)
    if (.not. passed) STOP 3

    call test_ten_2D4O3sym_mul(passed)
    if (.not. passed) STOP 4

    call test_ten_2D4O3sym_div(passed)
    if (.not. passed) STOP 5

    call test_ten_2D4O3sym_ddot(passed)
    if (.not. passed) STOP 6

    ! TODO tdot!!!!

    STOP 0
end program test_2D4O3sym

subroutine test_ten_2D4O3sym_approx(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1, to_test2

    call to_test1%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    call to_test2%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    passed = to_test1 .approx. to_test2
    if (.not. passed) return

    call to_test2%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10.001D0                           &
                       )
    passed = .not. (to_test1 .approx. to_test2)
    if (.not. passed) return

    call to_test1%init(xxxx=0D0, yyyy=0D0, zzzz=0D0, xyxy=0D0, &
                       xxyy=0D0, yyzz=0D0, zzxy=0D0,           &
                       xxzz=0D0, yyxy=0D0,                     &
                       xxxy=0D0                                &
                       )
    call to_test2%init(xxxx=0D0, yyyy=0D0, zzzz=0D0, xyxy=0D0, &
                       xxyy=0D0, yyzz=0D0, zzxy=0D0,           &
                       xxzz=0D0, yyxy=0D0,                     &
                       xxxy=0D0                                &
                       )
    passed = to_test1 .approx. to_test2
    if (.not. passed) return
end subroutine

subroutine test_ten_2D4O3sym_sum(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1, to_test2
    type(ten_2D4O3sym) :: expected_result

    call to_test1%init(xxxx=0D0, yyyy=0D0, zzzz=0D0, xyxy=0D0, &
                       xxyy=0D0, yyzz=0D0, zzxy=0D0,           &
                       xxzz=0D0, yyxy=0D0,                     &
                       xxxy=0D0                                &
                       )
    passed = (to_test1 + to_test1) .approx. to_test1
    if (.not. passed) return
    
    call to_test2%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    passed = (to_test1 + to_test2) .approx. to_test2
    if (.not. passed) return

    passed = (to_test2 + to_test1) .approx. to_test2
    if (.not. passed) return

    call to_test1%init(xxxx=1D1, yyyy=2D1, zzzz=3D1, xyxy=4D1, &
                       xxyy=5D1, yyzz=6D1, zzxy=7D1,           &
                       xxzz=8D1, yyxy=9D1,                     &
                       xxxy=10D1                               &
                       )
    call expected_result%init(xxxx=11D0, yyyy=22D0, zzzz=33D0, xyxy=44D0, &
                              xxyy=55D0, yyzz=66D0, zzxy=77D0,            &
                              xxzz=88D0, yyxy=99D0,                       &
                              xxxy=110D0                                  &
                              )
    passed = (to_test2 + to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1 + to_test2) .approx. expected_result
    if (.not. passed) return

end subroutine

subroutine test_ten_2D4O3sym_sub(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1, to_test2
    type(ten_2D4O3sym) :: expected_result

    call to_test1%init(xxxx=0D0, yyyy=0D0, zzzz=0D0, xyxy=0D0, &
                       xxyy=0D0, yyzz=0D0, zzxy=0D0,           &
                       xxzz=0D0, yyxy=0D0,                     &
                       xxxy=0D0                                &
                       )
    passed = (to_test1 - to_test1) .approx. to_test1
    if (.not. passed) return
    
    call to_test2%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    passed = (to_test2 - to_test1) .approx. to_test2
    if (.not. passed) return

    passed = (to_test1 - to_test2) .approx. (-to_test2)
    if (.not. passed) return

    call to_test1%init(xxxx=11D0, yyyy=12D0, zzzz=13D0, xyxy=14D0, &
                       xxyy=15D0, yyzz=16D0, zzxy=17D0,            &
                       xxzz=18D0, yyxy=19D0,                       &
                       xxxy=20D0                                   &
                       )
    call expected_result%init(xxxx=10D0, yyyy=10D0, zzzz=10D0, xyxy=10D0, &
                              xxyy=10D0, yyzz=10D0, zzxy=10D0,            &
                              xxzz=10D0, yyxy=10D0,                       &
                              xxxy=10D0                                   &
                              )
    passed = (to_test1 - to_test2) .approx. expected_result
    if (.not. passed) return

    passed = (to_test2 - to_test1) .approx. (-expected_result)
    if (.not. passed) return
end subroutine

subroutine test_ten_2D4O3sym_mul(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1
    type(ten_2D4O3sym) :: expected_result

    call to_test1%init(xxxx=0D0, yyyy=0D0, zzzz=0D0, xyxy=0D0, &
                       xxyy=0D0, yyzz=0D0, zzxy=0D0,           &
                       xxzz=0D0, yyxy=0D0,                     &
                       xxxy=0D0                                &
                       )
    passed = (10D0*to_test1) .approx. to_test1
    if (.not. passed) return

    passed = (to_test1*10D0) .approx. to_test1
    if (.not. passed) return

    call to_test1%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    call expected_result%init(xxxx=0D0, yyyy=0D0, zzzz=0D0, xyxy=0D0, &
                              xxyy=0D0, yyzz=0D0, zzxy=0D0,           &
                              xxzz=0D0, yyxy=0D0,                     &
                              xxxy=0D0                                &
                              )
    passed = (0D0*to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1*0D0) .approx. expected_result
    if (.not. passed) return

    call expected_result%init(xxxx=2D0, yyyy=4D0, zzzz=6D0, xyxy=8D0, &
                              xxyy=10D0, yyzz=12D0, zzxy=14D0,        &
                              xxzz=16D0, yyxy=18D0,                   &
                              xxxy=20D0                               &
                              )
    passed = (2D0*to_test1) .approx. expected_result
    if (.not. passed) return

    passed = (to_test1*2D0) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_ten_2D4O3sym_div(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1
    type(ten_2D4O3sym) :: expected_result

    call to_test1%init(xxxx=0D0, yyyy=0D0, zzzz=0D0, xyxy=0D0, &
                       xxyy=0D0, yyzz=0D0, zzxy=0D0,           &
                       xxzz=0D0, yyxy=0D0,                     &
                       xxxy=0D0                                &
                       )
    passed = (to_test1/2D0) .approx. to_test1
    if (.not. passed) return

    call to_test1%init(xxxx=2D0, yyyy=4D0, zzzz=6D0, xyxy=8D0, &
                       xxyy=10D0, yyzz=12D0, zzxy=14D0,        &
                       xxzz=16D0, yyxy=18D0,                   &
                       xxxy=20D0                               &
                       )
    call expected_result%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                              xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                              xxzz=8D0, yyxy=9D0,                     &
                              xxxy=10D0                               &
                              )
    passed = (to_test1/2D0) .approx. expected_result
    if (.not. passed) return
end subroutine

subroutine test_ten_2D4O3sym_ddot(passed)
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1
    type(ten_2D2Osym) :: to_test2, expected_result
    
    passed = .false.

    call to_test1%init(xxxx=0D0, yyyy=0D0, zzzz=0D0, xyxy=0D0, &
                       xxyy=0D0, yyzz=0D0, zzxy=0D0,           &
                       xxzz=0D0, yyxy=0D0,                     &
                       xxxy=0D0                                &
                       )
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call expected_result%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0)
    if ((to_test1 .ddot. to_test2) .approxexpected_result) passed = .true.
    if (.not. passed) return

    passed = .false.
    if ((to_test2 .ddot. to_test1) .approx. expected_result) passed = .true.
    if (.not. passed) return

    passed = .false.
    call to_test1%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                                &
                       )
    call expected_result%init(xx=115D0, yy=99D0, zz=85D0, xy=81D0)
    if ((to_test1 .ddot. to_test2) .approx. expected_result) passed = .true.
    if (.not. passed) return

    passed = .false.
    if ((to_test2 .ddot. to_test1) .approx. expected_result) passed = .true.
    if (.not. passed) return
end subroutine