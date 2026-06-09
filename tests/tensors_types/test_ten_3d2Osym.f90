program test_3D2Osym
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

    call test_ten_3D2Osym_components(passed)
    if (.not. passed) STOP 8

    call test_ten_3D2Osym_dot(passed)
    if (.not. passed) STOP 9

    call test_ten_3D2Osym_square(passed)
    if (.not. passed) STOP 10

    call test_ten_3D2Osym_assign(passed)
    if (.not. passed) STOP 11

    STOP 0
end program test_3D2Osym

subroutine test_ten_3D2Osym_approx(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2

    ! 1.- Test identity
    call to_test1%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0)
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6D0)
    passed = to_test1 .approx. to_test2
    if (.not. passed) then
        print*, "1.- Error: Approx failed for identical tensors", new_line('A'), &
                "T1:", to_test1, new_line('A'), &
                "T2:", to_test2
        return
    end if

    ! 2.- Test different values
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0, xz=5D0, yz=6.0001D0)
    passed = .not. (to_test1 .approx. to_test2)
    if (.not. passed) then
        print*, "2.- Error: Approx failed to detect difference", new_line('A'), &
                "T1:", to_test1, new_line('A'), &
                "T2:", to_test2
        return
    end if
    
    ! 3.- Test zero tensors
    call to_test1%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0)
    call to_test2%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, xz=0D0, yz=0D0)
    passed = to_test1 .approx. to_test2
    if (.not. passed) then
        print*, "3.- Error: Approx failed for zero tensors", new_line('A'), &
                "T1:", to_test1, new_line('A'), &
                "T2:", to_test2
        return
    end if
end subroutine

subroutine test_ten_3D2Osym_sum(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensor_3d2osym
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    type(ten_3D2Osym) :: expected, result

    ! 1.- Sum of two zero tensors
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result =  to_test1 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Sum of two zero tensors failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
    
    ! 2.- Sum with zero tensor (neutral element)
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    expected = to_test2
    result = to_test1 + to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Sum with zero tensor failed (0 + A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Commutativity with zero (A + 0)
    result = to_test2 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Sum with zero tensor failed (A + 0)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 4.- General sum of two tensors
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected%init((/12D0, 14D0, 16D0, 18D0, 20D0, 22D0/))
    result = to_test1 + to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: General sum failed (A + B)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
    
    ! 5.- Commutativity of general sum (B + A)
    result = to_test2 + to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: General sum failed (B + A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

end subroutine

subroutine test_ten_3D2Osym_sub(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    type(ten_3D2Osym) :: expected, result

    ! 1.- Subtraction of two zero tensors
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = to_test1 - to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Subtraction of two zero tensors failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
    
    ! 2.- Subtraction with zero tensor (A - 0)
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
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

    ! 4.- General subtraction of two tensors
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    call expected%init((/10D0, 10D0, 10D0, 10D0, 10D0, 10D0/))
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

subroutine test_ten_3D2Osym_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected, result

    ! 1.- Multiplication of zero tensor by scalar (10 * 0)
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = 10D0 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Multiplication of zero tensor failed (10 * 0)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 2.- Multiplication of zero tensor by scalar (0 * 10)
    result = to_test1 * 10D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Multiplication of zero tensor failed (0 * 10)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Multiplication by zero scalar (0 * A)
    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
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
    call expected%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0/))
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

subroutine test_ten_3D2Osym_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected, result

    ! 1.- Division of zero tensor
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = to_test1 / 2D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Division of zero tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 2.- General scalar division
    call to_test1%init((/2D0, 4D0, 6D0, 8D0, 10D0, 12D0/))
    call expected%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    result = to_test1 / 2D0
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: General scalar division failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
end subroutine

subroutine test_ten_3D2Osym_dev(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    type(ten_3D2Osym) :: expected, result

    ! 1.- Deviator of zero tensor
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = to_test1
    result = .dev. to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Deviator of zero tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 2.- Deviator of purely hydrostatic tensor
    call to_test1%init((/5D0, 5D0, 5D0, 0D0, 0D0, 0D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result = .dev. to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Deviator of hydrostatic tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Deviator with shear components
    call to_test1%init((/5D0, 5D0, 5D0, 1D0, 1D0, 1D0/))
    call expected%init((/0D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    result = .dev. to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Deviator with shear components failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 4.- General deviatoric case
    call to_test1%init((/3D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    call expected%init((/2D0, -1D0, -1D0, 1D0, 1D0, 1D0/))
    result = .dev. to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: General deviatoric case failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if
end subroutine

subroutine test_ten_3D2Osym_ddot(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    real(real64) :: result, expected
    real(real64), parameter :: EPS = 1D-10


    ! 1.- Double dot product of zero tensors
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    expected = 0D0
    result = to_test1 .ddot. to_test1
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "1.- Error: Double dot product of zeros failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 2.- Double dot product with one zero tensor
    call to_test2%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    expected = 0D0
    result = to_test1 .ddot. to_test2
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "2.- Error: Double dot product with zero failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 3.- General double dot product (A .ddot. B)
    call to_test1%init((/11D0, 12D0, 13D0, 14D0, 15D0, 16D0/))
    expected = 528D0
    result = to_test1 .ddot. to_test2
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "3.- Error: General double dot product failed (A : B)", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 4.- Commutativity of double dot product (B .ddot. A)
    result = to_test2 .ddot. to_test1
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "4.- Error: Double dot product commutativity failed (B : A)", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if
end subroutine


subroutine test_ten_3D2Osym_components(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1
    real(real64) :: expected, result
    real(real64), parameter :: EPS = 1D-10
    
    ! Initialize with distinct values for each component
    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))

    ! 1.- Test xx component
    expected = 1D0
    result = to_test1%xx()
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "1.- Error: Component xx() failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 2.- Test yy component
    expected = 2D0
    result = to_test1%yy()
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "2.- Error: Component yy() failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 3.- Test zz component
    expected = 3D0
    result = to_test1%zz()
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "3.- Error: Component zz() failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 4.- Test xy component
    expected = 4D0
    result = to_test1%xy()
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "4.- Error: Component xy() failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 5.- Test yz component
    expected = 5D0
    result = to_test1%yz()
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "5.- Error: Component yz() failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

    ! 6.- Test xz component
    expected = 6D0
    result = to_test1%xz()
    passed = abs(result - expected) < EPS
    if (.not. passed) then
        print*, "6.- Error: Component xz() failed", new_line('A'), &
                "Expected:", expected, " Actual:", result
        return
    end if

end subroutine


subroutine test_ten_3D2Osym_dot(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, to_test2
    type(ten_3D2O) :: expected, result
    
    ! 1.- Product of zero tensors
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Product of zero tensors failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 2.- Product of identity-like tensors (unit diagonal)
    call to_test1%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    call expected%init((/1D0, 0D0, 0D0, 0D0, 1D0, 0D0, 0D0, 0D0, 1D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Product of identity-like tensors failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Product with scaled diagonal
    call to_test1%init((/2D0, 2D0, 2D0, 0D0, 0D0, 0D0/))
    call expected%init((/4D0, 0D0, 0D0, 0D0, 4D0, 0D0, 0D0, 0D0, 4D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Product with scaled diagonal failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 4.- Product of tensors filled with ones
    call to_test1%init((/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected%init((/3D0, 3D0, 3D0, 3D0, 3D0, 3D0, 3D0, 3D0, 3D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: Product of tensors with ones failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 5.- Self-product with general symmetric components
    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected%init((/53D0, 42D0, 44D0, 42D0, 45D0, 49D0, 44D0, 49D0, 70D0/))
    result = to_test1 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "5.- Error: Self-product with general components failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 6.- Product with a zero tensor (neutral element test)
    call to_test2%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result = to_test1 * to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "6.- Error: Product with zero tensor failed (A * 0)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 7.- Commutativity with zero tensor (0 * A)
    result = to_test2 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "7.- Error: Product with zero tensor failed (0 * A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 8.- General dot product (A * B)
    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call to_test2%init((/1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call expected%init((/11D0, 11D0, 14D0, 11D0, 11D0, 14D0, 11D0, 11D0, 14D0/))
    result = to_test1 * to_test2
    passed = result .approx. expected
    if (.not. passed) then
        print*, "8.- Error: General dot product failed (A * B)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 9.- Non-commutativity test (B * A)
    call expected%init((/11D0, 11D0, 11D0, 11D0, 11D0, 11D0, 14D0, 14D0, 14D0/))
    result = to_test2 * to_test1
    passed = result .approx. expected
    if (.not. passed) then
        print*, "9.- Error: General dot product failed (B * A)", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

end subroutine


subroutine test_ten_3D2Osym_square(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, expected, result
    
    ! 1.- Square of a zero tensor
    call to_test1%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result = to_test1%square()
    passed = result .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Square of zero tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 2.- Square of identity-like tensor
    call to_test1%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    call expected%init((/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    result = to_test1%square()
    passed = result .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Square of identity-like tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 3.- Square of scaled diagonal tensor
    call to_test1%init((/2D0, 2D0, 2D0, 0D0, 0D0, 0D0/))
    call expected%init((/4D0, 4D0, 4D0, 0D0, 0D0, 0D0/))
    result = to_test1%square()
    passed = result .approx. expected
    if (.not. passed) then
        print*, "3.- Error: Square of scaled diagonal failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

    ! 4.- Square of a general symmetric tensor
    call to_test1%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    ! Order: xx, yy, zz, xy, xz, yz
    call expected%init((/53D0, 45D0, 70D0, 42D0, 49D0, 44D0/))
    result = to_test1%square()
    passed = result .approx. expected
    if (.not. passed) then
        print*, "4.- Error: Square of general symmetric tensor failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", result
        return
    end if

end subroutine


subroutine test_ten_3D2Osym_assign(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym) :: to_test1, expected
    real(real64) :: val_to_assign

    ! 1.- Assignment of zero
    val_to_assign = 0D0
    call expected%init((/0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    to_test1 = val_to_assign
    passed = to_test1 .approx. expected
    if (.not. passed) then
        print*, "1.- Error: Assignment of zero failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", to_test1
        return
    end if

    ! 2.- Assignment of a general value
    val_to_assign = 5.5D0
    call expected%init((/5.5D0, 5.5D0, 5.5D0, 5.5D0, 5.5D0, 5.5D0/))
    to_test1 = val_to_assign
    passed = to_test1 .approx. expected
    if (.not. passed) then
        print*, "2.- Error: Assignment of scalar value failed", new_line('A'), &
                "Expected:", expected, new_line('A'), &
                "Actual:", to_test1
        return
    end if
end subroutine