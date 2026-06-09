program test_3D2O
    use muscle_tensors
    implicit none
    logical :: passed

    ! real*8 :: a
    ! call test123(a, 5D0)
    call test_ten_3D2O_approx(passed)
    if (.not. passed) STOP 1

    call test_ten_3D2O_sum(passed)
    if (.not. passed) STOP 2

    call test_ten_3D2O_sub(passed)
    if (.not. passed) STOP 3

    call test_ten_3D2O_mul(passed)
    if (.not. passed) STOP 4

    call test_ten_3D2O_div(passed)
    if (.not. passed) STOP 5

    call test_ten_3D2O_dev(passed)
    if (.not. passed) STOP 6

    call test_ten_3D2O_ddot(passed)
    if (.not. passed) STOP 7

    call test_ten_3D2O_assign(passed)
    if (.not. passed) STOP 8

    STOP 0
end program test_3D2O


subroutine test_ten_3D2O_approx(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2O) :: to_test1, to_test2

    ! 1.- Test identity for general tensors
    call to_test1%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0, yx=7D0, zx=8D0, zy=9D0)
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0, yx=7D0, zx=8D0, zy=9D0)

    passed = to_test1 .approx. to_test2
    if (.not. passed) then
        print*, "1.- Error: Approx failed for identical general tensors", new_line('A'), &
                "T1:", to_test1, new_line('A'), &
                "T2:", to_test2
        return
    end if

    ! 2.- Test different values (detecting small differences)
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6.0001D0, yx=7D0, zx=8D0, zy=9D0)
    passed = .not. (to_test1 .approx. to_test2)
    if (.not. passed) then
        print*, "2.- Error: Approx failed to detect difference in general tensors", new_line('A'), &
                "T1:", to_test1, new_line('A'), &
                "T2:", to_test2
        return
    end if

    ! 3.- Test zero tensors
    call to_test1%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0, yx=0D0, zx=0D0, zy=0D0)
    call to_test2%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0, yx=0D0, zx=0D0, zy=0D0)
    passed = to_test1 .approx. to_test2
    if (.not. passed) then
        print*, "3.- Error: Approx failed for zero general tensors", new_line('A'), &
                "T1:", to_test1, new_line('A'), &
                "T2:", to_test2
        return
    end if
end subroutine


subroutine test_ten_3D2O_sum(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2O) :: to_test1, to_test2
    type(ten_3D2O) :: expected, result

    ! 1.- Sum of two zero tensors
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = to_test1 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Sum of two zero general tensors failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
    
    ! 2.- Sum with zero tensor (neutral element: 0 + A)
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0/))
    expected = to_test2
    result = to_test1 + to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Sum with zero failed (0 + A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Commutativity with zero (A + 0)
    result = to_test2 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Sum with zero failed (A + 0)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 4.- General sum of two tensors (A + B)
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0, 17D0, 18D0, 19D0/))
    call expected%init((/12D0, 14D0, 16D0, 18D0, 20D0, 22D0, 24D0, 26D0, 28D0/))
    result = to_test2 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: General sum failed (A + B)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 5.- Commutativity of general sum (B + A)
    result = to_test1 + to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: General sum failed (B + A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

end subroutine

subroutine test_ten_3D2O_sub(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2O) :: to_test1, to_test2
    type(ten_3D2O) :: expected, result

    ! 1.- Subtraction of two zero tensors
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = to_test1 - to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Subtraction of two zero general tensors failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
    
    ! 2.- Subtraction with zero tensor (A - 0)
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0/))
    expected = to_test2
    result = to_test2 - to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Subtraction with zero failed (A - 0)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Subtraction from zero (0 - A = -A)
    expected = -to_test2
    result = to_test1 - to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Subtraction from zero failed (0 - A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 4.- General subtraction of two tensors (A - B)
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0, 17D0, 18D0, 19D0/))
    call expected%init((/10D0, 10D0, 10D0, 10D0, 10D0, 10D0, 10D0, 10D0, 10D0/))
    result = to_test1 - to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: General subtraction failed (A - B)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 5.- Opposite subtraction (B - A = -(A - B))
    expected = -expected
    result = to_test2 - to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: General subtraction failed (B - A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
end subroutine

subroutine test_ten_3D2O_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2O) :: to_test1
    type(ten_3D2O) :: expected, result

    ! 1.- Multiplication of zero tensor by scalar (10 * 0)
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = 10D0 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Multiplication of zero general tensor failed (10 * 0)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 2.- Multiplication of zero tensor by scalar (0 * 10)
    result = to_test1 * 10D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Multiplication of zero general tensor failed (0 * 10)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Multiplication by zero scalar (0 * A)
    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result = 0D0 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Multiplication by zero scalar failed (0 * A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 4.- Multiplication by zero scalar (A * 0)
    result = to_test1 * 0D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: Multiplication by zero scalar failed (A * 0)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 5.- General scalar multiplication (2 * A)
    call expected%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0, 14D0, 16D0, 18D0/))
    result = 2D0 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: General scalar multiplication failed (2 * A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 6.- General scalar multiplication (A * 2)
    result = to_test1 * 2D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "6.- Error: General scalar multiplication failed (A * 2)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
end subroutine

subroutine test_ten_3D2O_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2O) :: to_test1
    type(ten_3D2O) :: expected, result

    ! 1.- Division of zero general tensor
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = to_test1 / 2D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Division of zero general tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 2.- General scalar division for general tensor
    call to_test1%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0, 14D0, 16D0, 18D0/))
    call expected%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0/))
    result = to_test1 / 2D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: General scalar division failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
end subroutine

subroutine test_ten_3D2O_dev(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2O) :: to_test1
    type(ten_3D2O) :: expected, result

    ! 1.- Deviator of zero tensor
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = .dev. to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Deviator of zero general tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 2.- Deviator of purely hydrostatic tensor
    call to_test1%init((/5D0, 0D0, 0D0, &
                         0D0, 5D0, 0D0, &
                         0D0, 0D0, 5D0/))
    call expected%init((/0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0, &
                         0D0, 0D0, 0D0/))
    result = .dev. to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Deviator of hydrostatic general tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Deviator with non-symmetric off-diagonal components
    call to_test1%init((/5D0, 1D0, 2D0, &
                         4D0, 5D0, 3D0, &
                         5D0, 6D0, 5D0/))
    call expected%init((/0D0, 1D0, 2D0, &
                         4D0, 0D0, 3D0, &
                         5D0, 6D0, 0D0/))
    result = .dev. to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Deviator with non-symmetric components failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 4.- General deviatoric case (mean stress subtraction)
    call to_test1%init((/3D0, 1D0, 2D0, &
                         4D0, 0D0, 3D0, &
                         5D0, 6D0, 0D0/))
    call expected%init((/2D0,  1D0,  2D0, &
                         4D0, -1D0,  3D0, &
                         5D0,  6D0, -1D0/))
    result = .dev. to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: General deviatoric case failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
end subroutine

subroutine test_ten_3D2O_ddot(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2O) :: to_test1, to_test2
    real(real64) :: result, expected

    ! 1.- Double dot product of zero tensors
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = 0D0
    result = to_test1 .ddot. to_test1
    passed = abs(result - expected) < 1D-10
    if (.not. passed) then
        print*, "1.- Error: Double dot product of zeros failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if
    
    ! 2.- Double dot product with one zero tensor
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0/))
    expected = 0D0
    result = to_test1 .ddot. to_test2
    passed = abs(result - expected) < 1D-10
    if (.not. passed) then
        print*, "2.- Error: Double dot product with zero failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 3.- General double dot product (A : B)
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0, 17D0, 18D0, 19D0/))
    expected = 735D0
    result = to_test1 .ddot. to_test2
    passed = abs(result - expected) < 1D-10
    if (.not. passed) then
        print*, "3.- Error: General double dot product failed (A : B)", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 4.- Commutativity of double dot product (B : A)
    result = to_test2 .ddot. to_test1
    passed = abs(result - expected) < 1D-10
    if (.not. passed) then
        print*, "4.- Error: Double dot product commutativity failed (B : A)", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if
end subroutine

subroutine test_ten_3D2O_assign(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2O) :: to_test1, expected
    real(real64) :: val_to_assign

    ! 1.- Assignment of zero to all 9 components
    val_to_assign = 0D0
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    to_test1 = val_to_assign
    passed = to_test1 .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Assignment of zero failed for ten_3D2O", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", to_test1
        return
    end if

    ! 2.- Assignment of a general value (e.g., 7.7)
    val_to_assign = 7.7D0
    call expected%init((/7.7D0, 7.7D0, 7.7D0, 7.7D0, 7.7D0, 7.7D0, 7.7D0, 7.7D0, 7.7D0/))
    to_test1 = val_to_assign
    passed = to_test1 .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Assignment of scalar value failed for ten_3D2O", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", to_test1
        return
    end if
end subroutine