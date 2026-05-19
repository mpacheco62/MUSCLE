program test_operator_3d2osym_3d4o3sym
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical :: passed

    print*, "Running tests for operator_3d2osym_3d4o3sym..."

    call test_tdotsym_unary(passed)
    if (.not. passed) stop 1

    call test_tdotsym_binary(passed)
    if (.not. passed) stop 2

    print*, "All operator_3d2osym_3d4o3sym tests passed successfully!"
    stop 0
end program test_operator_3d2osym_3d4o3sym


subroutine test_tdotsym_unary(passed)
    !! Tests the .tdotsym. operator in its "unary-like" form: a .tdotsym. a
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    
    type(ten_3D2Osym) :: a
    type(ten_3D4O3sym) :: result, expected

    ! Initialize tensor 'a'
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, yz=5D0, xz=6D0)

    ! Calculate the expected result of a ⊗ a manually
    expected%vals(1)  = a%vals(1) * a%vals(1) ! 1*1=1
    expected%vals(2)  = a%vals(2) * a%vals(2) ! 2*2=4
    expected%vals(3)  = a%vals(3) * a%vals(3) ! 3*3=9
    expected%vals(4)  = a%vals(4) * a%vals(4) ! 4*4=16
    expected%vals(5)  = a%vals(5) * a%vals(5) ! 5*5=25
    expected%vals(6)  = a%vals(6) * a%vals(6) ! 6*6=36
    expected%vals(7)  = a%vals(1) * a%vals(2) ! 1*2=2
    expected%vals(8)  = a%vals(2) * a%vals(3) ! 2*3=6
    expected%vals(9)  = a%vals(3) * a%vals(4) ! 3*4=12
    expected%vals(10) = a%vals(4) * a%vals(5) ! 4*5=20
    expected%vals(11) = a%vals(5) * a%vals(6) ! 5*6=30
    expected%vals(12) = a%vals(1) * a%vals(3) ! 1*3=3
    expected%vals(13) = a%vals(2) * a%vals(4) ! 2*4=8
    expected%vals(14) = a%vals(3) * a%vals(5) ! 3*5=15
    expected%vals(15) = a%vals(4) * a%vals(6) ! 4*6=24
    expected%vals(16) = a%vals(1) * a%vals(4) ! 1*4=4
    expected%vals(17) = a%vals(2) * a%vals(5) ! 2*5=10
    expected%vals(18) = a%vals(3) * a%vals(6) ! 3*6=18
    expected%vals(19) = a%vals(1) * a%vals(5) ! 1*5=5
    expected%vals(20) = a%vals(2) * a%vals(6) ! 2*6=12
    expected%vals(21) = a%vals(1) * a%vals(6) ! 1*6=6
    
    ! Perform the operation using the overloaded operator
    result = a .tdotsym. a

    passed = result .approx. expected
    if (.not. passed) then
        print*, "Unary-like .tdotsym. (a .tdotsym. a) FAILED", new_line('A'), &
                "Expected:", expected%vals, new_line('A'), &
                "Actual  :", result%vals
        return
    end if

    
    result = .tdotsym. a

    passed = result .approx. expected
    if (.not. passed) then
        print*, "Unary-like .tdotsym. (a .tdotsym. a) FAILED", new_line('A'), &
                "Expected:", expected%vals, new_line('A'), &
                "Actual  :", result%vals
        return
    end if
end subroutine test_tdotsym_unary


subroutine test_tdotsym_binary(passed)
    !! Tests the .tdotsym. operator in its binary form: a .tdotsym. b
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit none
    logical, intent(out) :: passed
    
    type(ten_3D2Osym) :: a, b
    type(ten_3D4O3sym) :: result, expected

    ! Initialize tensors 'a' and 'b'
    call a%init(xx=1D0, yy=2D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call b%init(xx=0D0, yy=0D0, zz=0D0, xy=3D0, yz=0D0, xz=0D0)

    ! Calculate the expected result of 0.5 * (a ⊗ b + b ⊗ a) manually
    ! The only non-zero components will be C_1112, C_2212 and their symmetric counterparts
    expected = 0.0_real64
    
    ! C_1112 = 0.5 * (a_11*b_12 + b_11*a_12) = 0.5 * (1*3 + 0*0) = 1.5
    expected%vals(16) = 1.5D0
    
    ! C_2212 = 0.5 * (a_22*b_12 + b_22*a_12) = 0.5 * (2*3 + 0*0) = 3.0
    expected%vals(13) = 3.0D0
    
    ! Perform the operation using the overloaded operator
    result = a .tdotsym. b

    passed = result .approx. expected
    if (.not. passed) then
        print*, "Binary .tdotsym. (a .tdotsym. b) FAILED", new_line('A'), &
                "Expected:", expected%vals, new_line('A'), &
                "Actual  :", result%vals
        return
    end if
end subroutine test_tdotsym_binary