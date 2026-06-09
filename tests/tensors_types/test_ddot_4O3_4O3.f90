program test_ddot_4O3_4O3
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_tensor_3d4o3sym
    use muscle_tensor_3d4o2sym
    implicit none
    
    type(ten_3D4O3sym) :: A, B
    type(ten_3D4O2sym) :: res, expected
    logical :: passed
    real(real64), parameter :: tol = 1D-10
    integer :: i, j

    ! 1. Inicializar A y B con valores únicos (1 al 21)
    ! Esto asegura que cada componente aporte un valor distinto al cálculo
    call A%init((/ 1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0, 10D0, 11D0, &
                  12D0, 13D0, 14D0, 15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0 /))
    
    call B%init((/21D0, 20D0, 19D0, 18D0, 17D0, 16D0, 15D0, 14D0, 13D0, &
                  12D0, 11D0, 10D0, 9D0, 8D0, 7D0, 6D0, 5D0, 4D0, 3D0, &
                  2D0, 1D0 /))

    ! 2. Calcular el resultado esperado manualmente (usando la misma fórmula que en la función)
    call expected%init(xxxx= 594D0, xxyy= 885D0, xxzz=1224D0, xxxy=1551D0, xxyz=1626D0, xxxz=1377D0, &
                       yyxx= 555D0, yyyy= 741D0, yyzz=1020D0, yyxy=1320D0, yyyz=1425D0, yyxz=1239D0, &
                       zzxx= 630D0, zzyy= 756D0, zzzz= 891D0, zzxy=1095D0, zzyz=1188D0, zzxz=1050D0, &
                       xyxx= 759D0, xyyy= 858D0, xyzz= 897D0, xyxy= 924D0, xyyz= 951D0, xyxz= 834D0, &
                       yzxx= 966D0, yzyy=1095D0, yzzz=1122D0, yzxy=1083D0, yzyz= 906D0, yzxz= 711D0, &
                       xzxx=1179D0, xzyy=1371D0, xzzz=1446D0, xzxy=1428D0, xzyz=1173D0, xzxz= 777D0)

    ! 3. Operación a probar (tu implementación desenrollada)
    res = A .ddot. B

    ! 4. Verificación exhaustiva
    passed = .true.
    do i = 1, 6
        do j = 1, 6
            if (abs(res%vals(i,j) - expected%vals(i,j)) > tol) then
                print *, "Test FAILED at component (", i, ",", j, ")"
                print *, "Expected:", expected%vals(i,j)
                print *, "Actual  :", res%vals(i,j)
                passed = .false.
            end if
        end do
    end do

    if (passed) then
        print *, "Test PASSED: Double contraction (4th:4th) is mathematically correct."
    else
        stop 1
    end if

end program test_ddot_4O3_4O3