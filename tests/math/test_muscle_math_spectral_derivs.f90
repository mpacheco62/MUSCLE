module test_muscle_math_spectral_derivs_mod
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: first_eigval, second_eigval, third_eigval
contains

    pure function first_eigval(x) result(res)
        use, intrinsic :: iso_fortran_env
        use tensors_types, only : ten_3D2Osym
        use muscle_math_operations, only : eigenvals
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res
        real(real64) :: eigenvalues(3)
        eigenvalues = eigenvals(x)
        res = eigenvalues(1)
    end function

    pure function second_eigval(x) result(res)
        use, intrinsic :: iso_fortran_env
        use tensors_types, only : ten_3D2Osym
        use muscle_math_operations, only : eigenvals
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res
        real(real64) :: eigenvalues(3)
        eigenvalues = eigenvals(x)
        res = eigenvalues(2)
    end function


    pure function third_eigval(x) result(res)
        use, intrinsic :: iso_fortran_env
        use tensors_types, only : ten_3D2Osym
        use muscle_math_operations, only : eigenvals
        implicit none
        type(ten_3D2Osym), intent(in) :: x
        real(real64) :: res
        real(real64) :: eigenvalues(3)
        eigenvalues = eigenvals(x)
        res = eigenvalues(3)
    end function
end module test_muscle_math_spectral_derivs_mod


! ********************** PROGRAM TEST ************************************
program test_derivatives
    use muscle_math_derivatives
    implicit none
    
    logical :: passed

    call test_first_derivative(passed)
    if (.not. passed) STOP 1

    call test_second_derivative(passed)
    if (.not. passed) STOP 2

    STOP 0
end program test_derivatives

subroutine test_first_derivative(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_math_derivatives
    use tensors_types
    use test_muscle_math_spectral_derivs_mod
    use muscle_math_operations, only : eigenvals
    use muscle_math_spectral_derivs, only : dEigenvalues_dTensor
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym):: to_test
    type(ten_3D2Osym) :: expected, result, result_tmp(3)
    real(real64) :: eigenvalues(3)

    call to_test%init(vals=(/1D0, 2D0, 3D0, 0D0, 0D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=0.0D0, yy=0.0D0, zz=1.0D0, xy=0.0D0, xz=0.0D0, yz=0.0D0)
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 1 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=0.0D0, yy=1.0D0, zz=0.0D0, xy=0.0D0, xz=0.0D0, yz=0.0D0)
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 2 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=1.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, xz=0.0D0, yz=0.0D0)
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 3 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=1.0D0, yy=0.0D0, zz=0.0D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 4.1 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=0.0D0, yy=0.5D0, zz=0.5D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 5.1 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=0.0D0, yy=0.5D0, zz=0.5D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 6.1 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(xx=0.9330127D0, yy=0.01674682D0, zz=0.05024047D0, &
                      xy=-0.125D0, yz=-0.02900635D0, xz=0.21650635D0)
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=0.9330127D0, yy=0.01674682D0, zz=0.05024047D0, &
                       xy=-0.125D0, yz=-0.02900635D0, xz=0.21650635D0)
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 4.2 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=     0.02406487932339D0, yy=     0.16301789279773D0, zz=     0.81291722787888D0, &
                       xy=    -0.06263390389981D0, yz=     0.36403303903325D0, xz=    -0.13986691885078D0  &
                       ) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 5.2 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=     0.04292241808653D0, yy=     0.82023528279858D0, zz=     0.13684229911488D0, &
                       xy=     0.18763390348657D0, yz=    -0.33502668835975D0, xz=    -0.07663943093820D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 6.2 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(xx=0.9330127018922193D0, yy=0.0167468245269452D0, zz=0.0502404735808355D0, &
                       xy=-0.125D0, yz=-0.0290063509461097D0, xz=0.2165063509461096D0)
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=     0.93301270189222D0, yy=     0.01674682452695D0, zz=     0.05024047358084D0, &
                       xy=    -0.12500000000000D0, yz=    -0.02900635094611D0, xz=     0.21650635094611D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 4.3 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=     0.03349364905389D0, yy=     0.49162658773653D0, zz=     0.47487976320958D0, &
                       xy=     0.06250000000000D0, yz=     0.01450317547305D0, xz=    -0.10825317547305D0  &
                       ) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 5.3 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=     0.03349364905389D0, yy=     0.49162658773653D0, zz=     0.47487976320958D0, &
                       xy=     0.06250000000000D0, yz=     0.01450317547305D0, xz=    -0.10825317547305D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 6.3 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=0.5D0, yy=0.5D0, zz=0.0D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 7 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=0.5D0, yy=0.5D0, zz=0.0D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 8 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=0.0D0, yy=0.0D0, zz=1.0D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 9 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=1D0/3D0, yy=1D0/3D0, zz=1D0/3D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 10 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=1D0/3D0, yy=1D0/3D0, zz=1D0/3D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 11 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=1D0/3D0, yy=1D0/3D0, zz=1D0/3D0, xy=0.0D0, xz=0.0D0, yz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 12 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/0D0, 0D0, 0D0, 2D0, 0D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=0.5D0, yy=0.5D0, zz=0.0D0, &
                       xy=0.5D0, yz=0.0D0, xz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 13 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=0.0D0, yy=0.0D0, zz=1.0D0, &
                       xy=0.0D0, yz=0.0D0, xz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 14 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=0.5D0, yy=0.5D0, zz=0.0D0, &
                       xy=-0.5D0, yz=0.0D0, xz=0.0D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 15 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/0D0, 0D0, 0D0, 2D0, 1D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=0.4D0, yy=0.5D0, zz=0.1D0, &
                       xy=0.4472136D0, yz=0.2236068D0, xz=0.2D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 16 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=0.2D0, yy=0.0D0, zz=0.8D0, &
                       xy=0.0D0, yz=0.0D0, xz=-0.4D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 17 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=0.4D0, yy=0.5D0, zz=0.1D0, &
                       xy=-0.4472136D0, yz=-0.2236068D0, xz=0.2D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 18 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


        !###########################################################################################    
    call to_test%init(vals=(/0D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=1D0/3, yy=1D0/3, zz=1D0/3, &
                       xy=1D0/3, yz=1D0/3, xz=1D0/3) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 19 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=1D0/3, yy=1D0/3, zz=1D0/3, &
                       xy=-1D0/6, yz=-1D0/6, xz=-1D0/6) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 20 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=1D0/3, yy=1D0/3, zz=1D0/3, &
                       xy=-1D0/6, yz=-1D0/6, xz=-1D0/6) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 21 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=0.29360201D0, yy=0.28545584D0, zz=0.42094215D0, &
                       xy=0.28950027D0, yz=0.34664160D0, xz=0.35155293D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 22 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=0.16085805D0, yy=0.71023177D0, zz=0.12891018D0, &
                       xy=-0.33800370D0, yz=-0.30258239D0, xz=0.14400084D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 23 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=0.54553994D0, yy=0.00431239D0, zz=0.45014767D0, &
                       xy=0.04850342D0, yz=-0.04405920D0, xz=-0.49555376D0) 
    passed = result%is_approx(expected, tol=1.0D-7)

    if (.not. passed) print*, "Case 24 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    !###########################################################################################    
    call to_test%init(vals=(/1D0, 0D0, 1D-6, 0D0, 0D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=     1.00000000000000D0, yy=     0.00000000000000D0, zz=     0.00000000000000D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-9)

    if (.not. passed) print*, "Case 25 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=     0.00000000000000D0, yy=     0.00000000000000D0, zz=     1.00000000000000D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-9)

    if (.not. passed) print*, "Case 26 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=     0.00000000000000D0, yy=     1.00000000000000D0, zz=     0.00000000000000D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-9)

    if (.not. passed) print*, "Case 27 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/1D0, 1D0, 1D0 - 1D-5, 0D0, 0D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=     0.50000000000000D0, yy=     0.50000000000000D0, zz=     0.00000000000000D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-5)

    if (.not. passed) print*, "Case 28 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=     0.50000000000000D0, yy=     0.50000000000000D0, zz=     0.00000000000000D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-5)

    if (.not. passed) print*, "Case 29 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=     0.00000000000000D0, yy=     0.00000000000000D0, zz=     1.00000000000000D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-5)

    if (.not. passed) print*, "Case 30 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/1D0, 1D0, 1D0 - 1D-6, 0D0, 0D0, 0D0/))
    eigenvalues = eigenvals(to_test)
    result_tmp = dEigenvalues_dTensor(to_test, eigenvalues)


    result = result_tmp(1) ! Extract the first eigenvalue derivative tensor
    call expected%init(xx=     1D0/3D0, yy=     1D0/3D0, zz=     1D0/3D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-8)

    if (.not. passed) print*, "Case 31 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return



    result = result_tmp(2) ! Extract the second eigenvalue derivative tensor
    call expected%init(xx=     1D0/3D0, yy=     1D0/3D0, zz=     1D0/3D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-8)

    if (.not. passed) print*, "Case 32 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    result = result_tmp(3) ! Extract the third eigenvalue derivative tensor
    call expected%init(xx=     1D0/3D0, yy=     1D0/3D0, zz=     1D0/3D0, &
                       xy=     0.00000000000000D0, yz=     0.00000000000000D0, xz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-8)

    if (.not. passed) print*, "Case 33 Spectral Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return
end subroutine





subroutine test_second_derivative(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_math_derivatives
    use tensors_types
    use test_muscle_math_spectral_derivs_mod
    use muscle_math_operations, only : eigenvals
    use muscle_math_spectral_derivs, only : dEigenvalues_dTensor, d2Eigenvalues_dTensor2
    implicit none
    
    logical, intent(out) :: passed

    type(ten_3D2Osym):: to_test
    type(ten_3D2Osym):: result_first_derivative(3)
    type(ten_3D4O3sym) :: expected, result, result_tmp(3)
    real(real64) :: eigenvalues(3)

    call to_test%init(vals=(/1D0, 2D0, 3D0, 0D0, 0D0, 0D0/))
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.0D0, yzyz= 0.5D0, xzxz= 0.25D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    eigenvalues = eigenvals(to_test)
    result_first_derivative = dEigenvalues_dTensor(to_test, eigenvalues)
    result_tmp = d2Eigenvalues_dTensor2(to_test, eigenvalues, result_first_derivative)
    result = result_tmp(1) ! Extract the first eigenvalue second derivative tensor
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 1 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(2) ! Extract the second eigenvalue second derivative tensor
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.5D0, yzyz= -0.5D0, xzxz= 0.0D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 2 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(3) ! Extract the third eigenvalue second derivative tensor
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= -0.5D0, yzyz= 0.0D0, xzxz= -0.25D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 3 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    

    call to_test%init(vals=(/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.5D0, yzyz= 0.0D0, xzxz= 0.5D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    eigenvalues = eigenvals(to_test)
    result_first_derivative = dEigenvalues_dTensor(to_test, eigenvalues)
    result_tmp = d2Eigenvalues_dTensor2(to_test, eigenvalues, result_first_derivative)
    result = result_tmp(1) ! Extract the first eigenvalue second derivative tensor
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 4.1 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(2) ! Extract the second eigenvalue second derivative tensor
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= -0.25D0, yzyz= 0.0D0, xzxz= -0.25D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 5.1 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(3) ! Extract the third eigenvalue second derivative tensor
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= -0.25D0, yzyz= 0.0D0, xzxz= -0.25D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 6.1 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return


    !###########################################################################################    
    call to_test%init(vals=(/1D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.0D0, yzyz= 0.25D0, xzxz= 0.25D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    eigenvalues = eigenvals(to_test)
    result_first_derivative = dEigenvalues_dTensor(to_test, eigenvalues)
    result_tmp = d2Eigenvalues_dTensor2(to_test, eigenvalues, result_first_derivative)
    result = result_tmp(1) ! Extract the first eigenvalue second derivative tensor
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 7 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(2) ! Extract the second eigenvalue second derivative tensor
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.0D0, yzyz= 0.25D0, xzxz= 0.25D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 8 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(3) ! Extract the third eigenvalue second derivative tensor
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.0D0, yzyz= -0.5D0, xzxz= -0.5D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 9 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/1D0, 1D0, 1D0, 0D0, 0D0, 0D0/))
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.0D0, yzyz= 0.0D0, xzxz= 0.0D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    eigenvalues = eigenvals(to_test)
    result_first_derivative = dEigenvalues_dTensor(to_test, eigenvalues)
    result_tmp = d2Eigenvalues_dTensor2(to_test, eigenvalues, result_first_derivative)
    result = result_tmp(1) ! Extract the first eigenvalue second derivative tensor
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 10 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(2) ! Extract the second eigenvalue second derivative tensor
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.0D0, yzyz= 0.0D0, xzxz= 0.0D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 11 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(3) ! Extract the third eigenvalue second derivative tensor
    call expected%init(xxxx= 0.0D0, yyyy= 0.0D0, zzzz= 0.0D0, &
                       xyxy= 0.0D0, yzyz= 0.0D0, xzxz= 0.0D0, &
                       xxyy= 0.0D0, yyzz= 0.0D0,              &
                       zzxy= 0.0D0, xyyz= 0.0D0, yzxz= 0.0D0, &
                       xxzz= 0.0D0,                           &
                       yyxy= 0.0D0, zzyz= 0.0D0, xyxz= 0.0D0, &
                       xxxy= 0.0D0, yyyz= 0.0D0, zzxz= 0.0D0, &
                       xxyz= 0.0D0, yyxz= 0.0D0, xxxz= 0.0D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 12 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/0D0, 0D0, 0D0, 2D0, 0D0, 0D0/))
    call expected%init(xxxx=     0.12500000000000D0, yyyy=     0.12500000000000D0, zzzz=     0.00000000000000D0, &
                       xyxy=    -0.00000000000000D0, yzyz=     0.12500000000000D0, xzxz=     0.12500000000000D0, &
                       xxyy=    -0.12500000000000D0, yyzz=     0.00000000000000D0,                               &
                       zzxy=     0.00000000000000D0, xyyz=     0.00000000000000D0, yzxz=     0.12500000000000D0, &
                       xxzz=     0.00000000000000D0,                                                             &
                       yyxy=     0.00000000000000D0, zzyz=     0.00000000000000D0, xyxz=     0.00000000000000D0, &
                       xxxy=     0.00000000000000D0, yyyz=     0.00000000000000D0, zzxz=     0.00000000000000D0, &
                       xxyz=     0.00000000000000D0, yyxz=    -0.00000000000000D0, xxxz=    -0.00000000000000D0  &
                       )

    eigenvalues = eigenvals(to_test)
    result_first_derivative = dEigenvalues_dTensor(to_test, eigenvalues)
    result_tmp = d2Eigenvalues_dTensor2(to_test, eigenvalues, result_first_derivative)
    result = result_tmp(1) ! Extract the first eigenvalue second derivative tensor
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 13 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(2) ! Extract the second eigenvalue second derivative tensor
    call expected%init(xxxx=     0.00000000000000D0, yyyy=     0.00000000000000D0, zzzz=     0.00000000000000D0, &
                       xyxy=     0.00000000000000D0, yzyz=     0.00000000000000D0, xzxz=    -0.00000000000000D0, &
                       xxyy=     0.00000000000000D0, yyzz=     0.00000000000000D0,                               &
                       zzxy=     0.00000000000000D0, xyyz=     0.00000000000000D0, yzxz=    -0.25000000000000D0, &
                       xxzz=     0.00000000000000D0,                                                             &
                       yyxy=     0.00000000000000D0, zzyz=     0.00000000000000D0, xyxz=     0.00000000000000D0, &
                       xxxy=     0.00000000000000D0, yyyz=     0.00000000000000D0, zzxz=     0.00000000000000D0, &
                       xxyz=     0.00000000000000D0, yyxz=     0.00000000000000D0, xxxz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 14 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(3) ! Extract the third eigenvalue second derivative tensor
    call expected%init(xxxx=    -0.12500000000000D0, yyyy=    -0.12500000000000D0, zzzz=     0.00000000000000D0, &
                       xyxy=     0.00000000000000D0, yzyz=    -0.12500000000000D0, xzxz=    -0.12500000000000D0, &
                       xxyy=     0.12500000000000D0, yyzz=     0.00000000000000D0,                               &
                       zzxy=     0.00000000000000D0, xyyz=     0.00000000000000D0, yzxz=     0.12500000000000D0, &
                       xxzz=     0.00000000000000D0,                                                             &
                       yyxy=     0.00000000000000D0, zzyz=     0.00000000000000D0, xyxz=     0.00000000000000D0, &
                       xxxy=     0.00000000000000D0, yyyz=     0.00000000000000D0, zzxz=     0.00000000000000D0, &
                       xxyz=     0.00000000000000D0, yyxz=     0.00000000000000D0, xxxz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 15 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/0D0, 0D0, 0D0, 2D0, 1D0, 0D0/))
    call expected%init(xxxx=     0.14310835055999D0, yyyy=     0.11180339887499D0, zzzz=     0.07602631123499D0, &
                       xyxy=     0.02236067977500D0, yzyz=     0.08944271909999D0, xzxz=     0.05813776741499D0, &
                       xxyy=    -0.08944271909999D0, yyzz=    -0.02236067977500D0,                               &
                       zzxy=    -0.04000000000000D0, xyyz=    -0.04472135955000D0, yzxz=     0.06000000000000D0, &
                       xxzz=    -0.05366563145999D0,                                                             &
                       yyxy=    -0.00000000000000D0, zzyz=     0.08000000000000D0, xyxz=    -0.03000000000000D0, &
                       xxxy=     0.04000000000000D0, yyyz=     0.00000000000000D0, zzxz=     0.06260990336999D0, &
                       xxyz=    -0.08000000000000D0, yyxz=    -0.04472135955000D0, xxxz=    -0.01788854382000D0  &
                       )

    eigenvalues = eigenvals(to_test)
    result_first_derivative = dEigenvalues_dTensor(to_test, eigenvalues)
    result_tmp = d2Eigenvalues_dTensor2(to_test, eigenvalues, result_first_derivative)
    result = result_tmp(1) ! Extract the first eigenvalue second derivative tensor
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 16 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(2) ! Extract the second eigenvalue second derivative tensor
    call expected%init(xxxx=     0.00000000000000D0, yyyy=     0.00000000000000D0, zzzz=     0.00000000000000D0, &
                       xyxy=    -0.00000000000000D0, yzyz=     0.00000000000000D0, xzxz=     0.00000000000000D0, &
                       xxyy=     0.00000000000000D0, yyzz=     0.00000000000000D0,                               &
                       zzxy=     0.08000000000000D0, xyyz=    -0.00000000000000D0, yzxz=    -0.12000000000000D0, &
                       xxzz=     0.00000000000000D0,                                                             &
                       yyxy=     0.00000000000000D0, zzyz=    -0.16000000000000D0, xyxz=     0.06000000000000D0, &
                       xxxy=    -0.08000000000000D0, yyyz=    -0.00000000000000D0, zzxz=     0.00000000000000D0, &
                       xxyz=     0.16000000000000D0, yyxz=     0.00000000000000D0, xxxz=     0.00000000000000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 17 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(3) ! Extract the third eigenvalue second derivative tensor
    call expected%init(xxxx=    -0.14310835055999D0, yyyy=    -0.11180339887499D0, zzzz=    -0.07602631123499D0, &
                       xyxy=    -0.02236067977500D0, yzyz=    -0.08944271909999D0, xzxz=    -0.05813776741499D0, &
                       xxyy=     0.08944271909999D0, yyzz=     0.02236067977500D0,                               &
                       zzxy=    -0.04000000000000D0, xyyz=     0.04472135955000D0, yzxz=     0.06000000000000D0, &
                       xxzz=     0.05366563145999D0,                                                             &
                       yyxy=    -0.00000000000000D0, zzyz=     0.08000000000000D0, xyxz=    -0.03000000000000D0, &
                       xxxy=     0.04000000000000D0, yyyz=     0.00000000000000D0, zzxz=    -0.06260990336999D0, &
                       xxyz=    -0.08000000000000D0, yyxz=     0.04472135955000D0, xxxz=     0.01788854382000D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 18 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/0D0, 0D0, 0D0, 1D0, 1D0, 1D0/))
    call expected%init(xxxx=     0.14814814814815D0, yyyy=     0.14814814814815D0, zzzz=     0.14814814814815D0, &
                       xyxy=     0.03703703703704D0, yzyz=     0.03703703703704D0, xzxz=     0.03703703703704D0, &
                       xxyy=    -0.07407407407407D0, yyzz=    -0.07407407407407D0,                               &
                       zzxy=    -0.07407407407407D0, xyyz=    -0.01851851851852D0, yzxz=    -0.01851851851852D0, &
                       xxzz=    -0.07407407407407D0,                                                             &
                       yyxy=     0.03703703703704D0, zzyz=     0.03703703703704D0, xyxz=    -0.01851851851852D0, &
                       xxxy=     0.03703703703704D0, yyyz=     0.03703703703704D0, zzxz=     0.03703703703704D0, &
                       xxyz=    -0.07407407407407D0, yyxz=    -0.07407407407407D0, xxxz=     0.03703703703704D0  &
                       )

    eigenvalues = eigenvals(to_test)
    result_first_derivative = dEigenvalues_dTensor(to_test, eigenvalues)
    result_tmp = d2Eigenvalues_dTensor2(to_test, eigenvalues, result_first_derivative)
    result = result_tmp(1) ! Extract the first eigenvalue second derivative tensor
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 19 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(2) ! Extract the second eigenvalue second derivative tensor
    call expected%init(xxxx=    -0.07407407407407D0, yyyy=    -0.07407407407407D0, zzzz=    -0.07407407407407D0, &
                       xyxy=    -0.01851851851852D0, yzyz=    -0.01851851851852D0, xzxz=    -0.01851851851852D0, &
                       xxyy=     0.03703703703704D0, yyzz=     0.03703703703704D0,                               &
                       zzxy=     0.03703703703704D0, xyyz=     0.00925925925926D0, yzxz=     0.00925925925926D0, &
                       xxzz=     0.03703703703704D0,                                                             &
                       yyxy=    -0.01851851851852D0, zzyz=    -0.01851851851852D0, xyxz=     0.00925925925926D0, &
                       xxxy=    -0.01851851851852D0, yyyz=    -0.01851851851852D0, zzxz=    -0.01851851851852D0, &
                       xxyz=     0.03703703703704D0, yyxz=     0.03703703703704D0, xxxz=    -0.01851851851852D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 20 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(3) ! Extract the third eigenvalue second derivative tensor
    call expected%init(xxxx=    -0.07407407407407D0, yyyy=    -0.07407407407407D0, zzzz=    -0.07407407407407D0, &
                       xyxy=    -0.01851851851852D0, yzyz=    -0.01851851851852D0, xzxz=    -0.01851851851852D0, &
                       xxyy=     0.03703703703704D0, yyzz=     0.03703703703704D0,                               &
                       zzxy=     0.03703703703704D0, xyyz=     0.00925925925926D0, yzxz=     0.00925925925926D0, &
                       xxzz=     0.03703703703704D0,                                                             &
                       yyxy=    -0.01851851851852D0, zzyz=    -0.01851851851852D0, xyxz=     0.00925925925926D0, &
                       xxxy=    -0.01851851851852D0, yyyz=    -0.01851851851852D0, zzxz=    -0.01851851851852D0, &
                       xxyz=     0.03703703703704D0, yyxz=     0.03703703703704D0, xxxz=    -0.01851851851852D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 21 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    !###########################################################################################    
    call to_test%init(vals=(/1D0, 2D0, 3D0, 4D0, 5D0, 6D0/))
    call expected%init(xxxx=     0.02641587237785D0, yyyy=     0.02878290798163D0, zzzz=     0.03102337600254D0, &
                       xyxy=     0.00777796019660D0, yzyz=     0.00752295574214D0, xzxz=     0.00771374176081D0, &
                       xxyy=    -0.01208770217847D0, yyzz=    -0.01669520580316D0,                               &
                       zzxy=    -0.01552987450603D0, xyyz=    -0.00115084570493D0, yzxz=    -0.00682487616557D0, &
                       xxzz=    -0.01432817019938D0,                                                             &
                       yyxy=     0.00863594067902D0, zzyz=     0.00263685626345D0, xyxz=    -0.00235761801587D0, &
                       xxxy=     0.00689393382702D0, yyyz=     0.01060200967800D0, zzxz=     0.00437657393432D0, &
                       xxyz=    -0.01323886594145D0, yyxz=    -0.01420834321562D0, xxxz=     0.00983176928130D0  &
                       )

    eigenvalues = eigenvals(to_test)
    result_first_derivative = dEigenvalues_dTensor(to_test, eigenvalues)
    result_tmp = d2Eigenvalues_dTensor2(to_test, eigenvalues, result_first_derivative)
    result = result_tmp(1) ! Extract the first eigenvalue second derivative tensor
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 22 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(2) ! Extract the second eigenvalue second derivative tensor
    call expected%init(xxxx=     0.07849539359020D0, yyyy=    -0.02565872801556D0, zzzz=     0.04865311490748D0, &
                       xyxy=     0.08414586429932D0, yzyz=     0.07972464915323D0, xzxz=    -0.00729772553812D0, &
                       xxyy=    -0.00209177533359D0, yyzz=     0.02775050334914D0,                               &
                       zzxy=     0.07366840220539D0, xyyz=    -0.08822722981537D0, yzxz=     0.00514141454776D0, &
                       xxzz=    -0.07640361825662D0,                                                             &
                       yyxy=     0.00830324990975D0, zzyz=    -0.06301144889056D0, xyxz=     0.00445544411988D0, &
                       xxxy=    -0.08197165211514D0, yyyz=    -0.02710272980136D0, zzxz=    -0.00702412182514D0, &
                       xxyz=     0.09011417869192D0, yyxz=     0.01456325167107D0, xxxz=    -0.00753912984594D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 23 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

    result = result_tmp(3) ! Extract the third eigenvalue second derivative tensor
    call expected%init(xxxx=    -0.10491126596805D0, yyyy=    -0.00312417996607D0, zzzz=    -0.07967649091001D0, &
                       xyxy=    -0.09192382449592D0, yzyz=    -0.08724760489537D0, xzxz=    -0.00041601622269D0, &
                       xxyy=     0.01417947751205D0, yyzz=    -0.01105529754598D0,                               &
                       zzxy=    -0.05813852769935D0, xyyz=     0.08937807552030D0, yzxz=     0.00168346161781D0, &
                       xxzz=     0.09073178845600D0,                                                             &
                       yyxy=    -0.01693919058877D0, zzyz=     0.06037459262711D0, xyxz=    -0.00209782610401D0, &
                       xxxy=     0.07507771828812D0, yyyz=     0.01650072012336D0, zzxz=     0.00264754789082D0, &
                       xxyz=    -0.07687531275047D0, yyxz=    -0.00035490845545D0, xxxz=    -0.00229263943536D0  &
                       )
    passed = result%is_approx(expected, tol=1.0D-10)

    if (.not. passed) print*, "Case 24 Spectral Second Derivate",  new_line('A'), &
                              "The values obtained is different from the expected one", new_line('A'), &
                              "The values obtained are:", result, new_line('A'), &
                              "The expected are:", expected, new_line('A'), &
                              "The differences:", result - expected, new_line('A')
    if (.not. passed) return

end subroutine