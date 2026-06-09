program test_elasticity_linear
    implicit none
    
    logical :: passed

    call test_tensile_nu0_3D(passed)
    if (.not. passed) STOP 1

    call test_tensile_nu0_2D(passed)
    if (.not. passed) STOP 2

    call test_tensile_3D(passed)
    if (.not. passed) STOP 3

    call test_tensile_2D(passed)
    if (.not. passed) STOP 4

    call test_dsigma_dstrain_3D(passed)
    if (.not. passed) STOP 5

    call test_dsigma_dstrain_2D(passed)
    if (.not. passed) STOP 6

    print*, "Passed!", passed
    STOP 0
end program test_elasticity_linear

subroutine test_tensile_nu0_3D(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use mod_elasticity_linear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Elasticity_linear) :: el
    type(ten_3D2Osym) :: expected_result, result 
    type(ten_3D2Osym) :: strain 

    call el%set_parameters(young=1D0, poisson=0D0)
    
    
    call expected_result%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call strain%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    call strain%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    call strain%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    call strain%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
    call strain%init((/0D0, 0D0, 0D0, 0D0, 1D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
    call strain%init((/0D0, 0D0, 0D0, 0D0, 0D0, 1D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return
end subroutine

subroutine test_tensile_nu0_2D(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use mod_elasticity_linear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Elasticity_linear) :: el
    type(ten_2D2Osym) :: expected_result, result 
    type(ten_2D2Osym) :: strain 

    call el%set_parameters(young=1D0, poisson=0D0)
    
    
    call expected_result%init((/1D0, 0D0, 0D0, 0D0/))
    call strain%init((/1D0, 0D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 1D0, 0D0, 0D0/))
    call strain%init((/0D0, 1D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 0D0, 1D0, 0D0/))
    call strain%init((/0D0, 0D0, 1D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 0D0, 0D0, 1D0/))
    call strain%init((/0D0, 0D0, 0D0, 1D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) print*, "Stress is not equal", new_line('A'),         &
                              "Expected:", expected_result, new_line('A'),  &
                              "Actual Value:", result, new_line('A'),       &
                              "Difference", result - expected_result
    if (.not. passed) return

end subroutine

subroutine test_tensile_3D(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use mod_elasticity_linear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Elasticity_linear) :: el
    type(ten_3D2Osym) :: expected_result, result 
    type(ten_3D2Osym) :: strain 

    call el%set_parameters(young=1D0, poisson=0.3D0)
    
    call expected_result%init((/1D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call strain%init((/1D0, -0.3D0, -0.3D0, 0D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) then
        print*, "1.- Error: Linear tensile result is not expected", new_line('A'), &
                "Expected:", expected_result, new_line('A'), &
                "Actual:", result
        return
    end if
    
    call expected_result%init((/0D0, 1D0, 0D0, 0D0, 0D0, 0D0/))
    call strain%init((/-0.3D0, 1D0, -0.3D0, 0D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) then
        print*, "2.- Error: Linear tensile result is not expected", new_line('A'), &
                "Expected:", expected_result, new_line('A'), &
                "Actual:", result
        return
    end if

    call expected_result%init((/0D0, 0D0, 1D0, 0D0, 0D0, 0D0/))
    call strain%init((/-0.3D0, -0.3D0, 1D0, 0D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) then
        print*, "3.- Error: Linear tensile result is not expected", new_line('A'), &
                "Expected:", expected_result, new_line('A'), &
                "Actual:", result
        return
    end if

    call expected_result%init((/0D0, 0D0, 0D0, 0.7692307692307692D0, 0D0, 0D0/))
    call strain%init((/0D0, 0D0, 0D0, 1D0, 0D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) then
        print*, "4.- Error: Linear tensile result is not expected", new_line('A'), &
                "Expected:", expected_result, new_line('A'), &
                "Actual:", result, new_line('A'), &
                "Diference:", result - expected_result
        return
    end if

end subroutine


subroutine test_tensile_2D(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use mod_elasticity_linear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Elasticity_linear) :: el
    type(ten_2D2Osym) :: expected_result, result 
    type(ten_2D2Osym) :: strain 

    call el%set_parameters(young=1D0, poisson=0.3D0)
    
    
    call expected_result%init((/1D0, 0D0, 0D0, 0D0/))
    call strain%init((/1D0, -0.3D0, -0.3D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) return
    
    call expected_result%init((/0D0, 1D0, 0D0, 0D0/))
    call strain%init((/-0.3D0, 1D0, -0.3D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 0D0, 1D0, 0D0/))
    call strain%init((/-0.3D0, -0.3D0, 1D0, 0D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) return

    call expected_result%init((/0D0, 0D0, 0D0, 0.769230769D0/))
    call strain%init((/0D0, 0D0, 0D0, 1D0/))
    result = el%stress(strain=strain)
    passed = result .approx. expected_result
    if (.not. passed) return

end subroutine


subroutine test_dsigma_dstrain_3D(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use mod_elasticity_linear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Elasticity_linear) :: el
    type(ten_3D4O3sym) :: expected_result, result 
    type(ten_3D2Osym) :: strain

    call el%set_parameters(young=1D0, poisson=0.0D0, calc_tan_3D=.True.)

    call expected_result%init(xxxx=1D0, yyyy=1D0, zzzz=1D0,       & 
                              xyxy=0.5D0, yzyz=0.5D0, xzxz=0.5D0, &
                              xxyy=0D0, yyzz=0D0, zzxy=0D0,       &
                              xyyz=0D0, yzxz=0D0, xxzz=0D0,       &
                              yyxy=0D0, zzyz=0D0, xyxz=0D0,       &
                              xxxy=0D0, yyyz=0D0, zzxz=0D0,       &
                              xxyz=0D0, yyxz=0D0, xxxz=0D0)
    result = el%dstress_dstrain(strain)
    passed = result .approx. expected_result
    if (.not. passed) return


    call el%set_parameters(young=1D0, poisson=0.2D0, calc_tan_3D=.True.)

    call expected_result%init(xxxx=1.111111111111111D0, yyyy=1.111111111111111D0, zzzz=1.111111111111111D0,   & 
                              xyxy=0.4166666666666666D0, yzyz=0.4166666666666666D0, xzxz=0.4166666666666666D0, &
                              xxyy=0.2777777777777777D0, xxzz=0.2777777777777777D0, yyzz=0.2777777777777777D0, &
                              zzxy=0D0, xyyz=0D0, yzxz=0D0,        &
                              yyxy=0D0, zzyz=0D0, xyxz=0D0,       &
                              xxxy=0D0, yyyz=0D0, zzxz=0D0,       &
                              xxyz=0D0, yyxz=0D0, xxxz=0D0)
    result = el%dstress_dstrain(strain)
    passed = result .approx. expected_result
    if (.not. passed) return

end subroutine


subroutine test_dsigma_dstrain_2D(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use mod_elasticity_linear
    implicit none
    
    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(Elasticity_linear) :: el
    type(ten_2D4O3sym) :: expected_result, result 
    type(ten_2D2Osym) :: strain

    call el%set_parameters(young=1D0, poisson=0.0D0, calc_tan_2D=.True.)

    call expected_result%init(xxxx=1D0, yyyy=1D0, zzzz=1D0,       & 
                              xyxy=0.5D0, xxyy=0D0, yyzz=0D0,     &
                              zzxy=0D0, xxzz=0D0, yyxy=0D0,       &
                              xxxy=0D0                            & 
                              )
    result = el%dstress_dstrain(strain)
    passed = result .approx. expected_result
    if (.not. passed) return


    call el%set_parameters(young=1D0, poisson=0.2D0, calc_tan_2D=.True.)

    call expected_result%init(xxxx=1.111111111111111D0, yyyy=1.111111111111111D0, zzzz=1.111111111111111D0,    & 
                              xyxy=0.4166666666666666D0, xxyy=0.2777777777777777D0, xxzz=0.2777777777777777D0, &
                              yyzz=0.2777777777777777D0, zzxy=0D0, yyxy=0D0,                                   &
                              xxxy=0D0                                                                         &
                              )
    result = el%dstress_dstrain(strain)
    passed = result .approx. expected_result
    if (.not. passed) return

end subroutine