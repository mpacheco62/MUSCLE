! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_tensor_openmp
    use omp_lib
    use muscle_tensors
    implicit none

    logical :: passed

    print*, "========================================================"
    print*, "  RUNNING OPENMP THREAD-SAFETY TEST: Pure Tensor Engine "
    print*, "========================================================"

    call test_openmp_tensor_ops(passed)
    if (.not. passed) stop 1

    print*, "========================================================"
    print*, "  PASSED: Pure Tensor Engine is 100% Thread-Safe!       "
    print*, "========================================================"
    stop 0
end program test_muscle_tensor_openmp


subroutine test_openmp_tensor_ops(passed)
    use omp_lib
    use muscle_tensors
    implicit none

    logical, intent(out) :: passed
    integer, parameter   :: N = 100000
    integer, parameter   :: NUM_THREADS = 4

    type(ten_3D2Osym), allocatable :: a(:)
    type(ten_3D2Osym), allocatable :: seq_res(:)
    type(ten_3D2Osym), allocatable :: par_oper(:), par_fun(:), par_sub(:)
    type(ten_3D2Osym) :: s_loc

    integer :: i
    logical :: err_oper, err_fun, err_sub

    print*, "========================================================="
    print*, "  TEST MWE OPENMP: COMPARACIÓN DE 3 VARIANTES            "
    print*, "========================================================="
    write(*, '(A, I0)') "  Número de elementos : ", N
    write(*, '(A, I0)') "  Hilos OpenMP        : ", NUM_THREADS

    call omp_set_num_threads(NUM_THREADS)

    ! =========================================================================
    ! VERIFICACIÓN EFECTIVA DE OPENMP Y HILOS ACTIVOS
    ! =========================================================================
    print*, ""
    print*, "---------------------------------------------------------"
    print*, "  VERIFICACIÓN DE EJECUCIÓN PARALELA DE HILOS:"
    write(*, '(A, I0)') "   - OMP Max Threads (omp_get_max_threads) : ", omp_get_max_threads()

!$OMP PARALLEL
!$OMP SINGLE
    write(*, '(A, I0)') "   - Hilos REALES activos (omp_get_num_threads): ", omp_get_num_threads()
!$OMP END SINGLE

!$OMP CRITICAL
    write(*, '(A, I0, A)') "     -> [Hilo ", omp_get_thread_num(), "] activo y listo."
!$OMP END CRITICAL
!$OMP END PARALLEL
    print*, "---------------------------------------------------------"
    print*, ""

    allocate(a(N), seq_res(N), par_oper(N), par_fun(N), par_sub(N))

    ! 1. Inicializar datos de entrada con valores únicos
    do i = 1, N
        call a(i)%init(xx=10.0D0 * dble(i)/N, yy=5.0D0 * dble(i)/N, zz=2.0D0 * dble(i)/N, &
                       xy=1.0D0 * dble(i)/N, yz=0.5D0 * dble(i)/N, xz=-0.2D0 * dble(i)/N)
    end do

    ! 2. Referencia Secuencial
    do i = 1, N
        seq_res(i) = .dev. a(i)
    end do


!$OMP PARALLEL DO PRIVATE(i, s_loc) SHARED(a, par_oper, par_fun, par_sub)
    do i = 1, N
        ! Variante 1: Operador Unario .dev.
        par_oper(i) = .dev. a(i)

        ! Variante 2: Método Función
        par_fun(i)  = a(i)%dev()

        ! Variante 3: Subrutina Directa
        call a(i)%calc_dev(par_sub(i))
    end do
!$OMP END PARALLEL DO

    ! 4. Verificación de Resultados punto a punto
    err_oper = .false.
    err_fun  = .false.
    err_sub  = .false.

    passed = .true.

    do i = 1, N
        ! Verificar Variante 1 (Operador .dev.)
        if (any(abs(par_oper(i)%vals - seq_res(i)%vals) > 1.0D-10) .and. .not. err_oper) then
            err_oper = .true.
            passed = .false.
            print*, ""
            print*, "---------------------------------------------------------"
            print*, " [FALLA] VARIANTE 1 (Operador .dev.) falló en índice: ", i
            print*, "   a(i) : ", a(i)%vals
            print*, "   seq  : ", seq_res(i)%vals
            print*, "   par  : ", par_oper(i)%vals
            print*, "   diff : ", par_oper(i)%vals - seq_res(i)%vals
        end if

        ! Verificar Variante 2 (Método Función)
        if (any(abs(par_fun(i)%vals - seq_res(i)%vals) > 1.0D-10) .and. .not. err_fun) then
            err_fun = .true.
            passed = .false.
            print*, ""
            print*, "---------------------------------------------------------"
            print*, " [FALLA] VARIANTE 2 (Método Función) falló en índice: ", i
            print*, "   a(i) : ", a(i)%vals
            print*, "   seq  : ", seq_res(i)%vals
            print*, "   par  : ", par_fun(i)%vals
            print*, "   diff : ", par_fun(i)%vals - seq_res(i)%vals
        end if

        ! Verificar Variante 3 (Subrutina)
        if (any(abs(par_sub(i)%vals - seq_res(i)%vals) > 1.0D-10) .and. .not. err_sub) then
            err_sub = .true.
            passed = .false.
            print*, ""
            print*, "---------------------------------------------------------"
            print*, " [FALLA] VARIANTE 3 (Subrutina) falló en índice: ", i
            print*, "   a(i) : ", a(i)%vals
            print*, "   seq  : ", seq_res(i)%vals
            print*, "   par  : ", par_sub(i)%vals
            print*, "   diff : ", par_sub(i)%vals - seq_res(i)%vals
        end if
    end do

    print*, ""
    print*, "========================================================="
    print*, " RESUMEN DE RESULTADOS:"
    if (.not. err_oper) print*, "   Variante 1 (Operador .dev.) : PASÓ OK "
    if (err_oper)       print*, "   Variante 1 (Operador .dev.) : FALLÓ "

    if (.not. err_fun)  print*, "   Variante 2 (Método Función)  : PASÓ OK "
    if (err_fun)        print*, "   Variante 2 (Método Función)  : FALLÓ "

    if (.not. err_sub)  print*, "   Variante 3 (Subrutina)       : PASÓ OK "
    if (err_sub)        print*, "   Variante 3 (Subrutina)       : FALLÓ "
    print*, "========================================================="

    deallocate(a, seq_res, par_oper, par_fun, par_sub)

end subroutine test_openmp_tensor_ops