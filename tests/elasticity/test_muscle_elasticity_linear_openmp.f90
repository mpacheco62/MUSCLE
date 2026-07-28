! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_elasticity_linear_openmp
    use omp_lib
    use muscle_tensors
    use muscle_elasticity_linear
    implicit none

    logical :: passed

    print*, "========================================================"
    print*, "  RUNNING OPENMP THREAD-SAFETY TEST: Elasticity_linear  "
    print*, "========================================================"

    print*, "1. Testing Elasticity_linear with SHARED object..."
    call test_openmp_elasticity_linear_shared(passed)
    if (.not. passed) stop 1

    print*, "2. Testing Elasticity_linear with THREAD-PRIVATE objects..."
    call test_openmp_elasticity_linear_private(passed)
    if (.not. passed) stop 2

    print*, "========================================================"
    print*, "  PASSED: Elasticity_linear is 100% OpenMP Thread-Safe! "
    print*, "========================================================"
    stop 0
end program test_muscle_elasticity_linear_openmp


! =========================================================================
! 1. TEST CON ELASTICIDAD COMPARTIDA (SHARED)
! =========================================================================
subroutine test_openmp_elasticity_linear_shared(passed)
    use omp_lib
    use muscle_tensors
    use muscle_elasticity_linear
    implicit none

    logical, intent(out) :: passed
    integer, parameter   :: N = 100000
    integer, parameter   :: NUM_THREADS = 4

    type(ten_3D2Osym), allocatable  :: strains(:)
    type(ten_3D2Osym), allocatable  :: stress_seq(:), stress_par(:)
    type(ten_3D4O3sym), allocatable :: tan_seq(:), tan_par(:)
    type(ten_3D2Osym)  :: ten_tmp

    type(Elasticity_linear) :: elas
    integer :: i
    real(8) :: t_start, t_seq, t_par

    passed = .FALSE.

    allocate(strains(N), stress_seq(N), stress_par(N), tan_seq(N), tan_par(N))

    ! Configurar objeto elástico
    call elas%set_parameters(young=210000.0D0, poisson=0.3D0)

    ! Generar 100,000 deformaciones dinámicas
    do i = 1, N
        call strains(i)%init(xx=0.001D0 * sin(dble(i)), &
                             yy=-0.0003D0 * cos(dble(i)), &
                             zz=-0.0003D0 * sin(2.0D0*dble(i)), &
                             xy=0.0005D0 * cos(0.5D0*dble(i)), &
                             yz=0.0001D0 * sin(0.1D0*dble(i)), &
                             xz=-0.0002D0 * cos(0.1D0*dble(i)))
    end do

    ! --- A. EJECUCIÓN SECUENCIAL (1 HILO) ---
    t_start = omp_get_wtime()
    do i = 1, N
        stress_seq(i) = elas%stress(strains(i))
    end do
    t_seq = omp_get_wtime() - t_start

    ! --- B. EJECUCIÓN PARALELA CON OBJETO COMPARTIDO (SHARED) ---
    call omp_set_num_threads(NUM_THREADS)
    t_start = omp_get_wtime()

!$OMP PARALLEL DO PRIVATE(i, ten_tmp) SHARED(elas, strains, stress_par, tan_par)
    do i = 1, N
        stress_par(i) = elas%stress(strains(i))
    end do
!$OMP END PARALLEL DO

    t_par = omp_get_wtime() - t_start

    write(*, '(A, F8.4, A, F8.4, A)') "   Tiempo Secuencial: ", t_seq, "s | Tiempo Paralelo: ", t_par, "s"
    write(*, '(A, F6.2, A)') "   Aceleración (Speedup): ", t_seq / t_par, "x"

    ! --- C. VERIFICACIÓN 1:1 ---
    do i = 1, N
        if (.not. (stress_par(i) .approx. stress_seq(i))) then
            print*, "--------------------------------------------------------"
            print*, "   FAIL: Mismatch in stress_par at index i = ", i
            print*, "   Deformación entrada (strains):", strains(i)%vals
            print*, "   Esfuerzo Secuencial (seq)   :", stress_seq(i)%vals
            print*, "   Esfuerzo Paralelo   (par)   :", stress_par(i)%vals
            print*, "   Diferencia (par - seq)      :", stress_par(i)%vals - stress_seq(i)%vals
            print*, "--------------------------------------------------------"
            return
        end if
        if (.not. (tan_par(i) .approx. tan_seq(i))) then
            print*, "   FAIL: Tangent mismatch in shared elasticity at point ", i
            return
        end if
    end do

    deallocate(strains, stress_seq, stress_par, tan_seq, tan_par)
    passed = .TRUE.
end subroutine test_openmp_elasticity_linear_shared


! =========================================================================
! 2. TEST CON ELASTICIDAD PRIVADA (THREAD-PRIVATE)
! =========================================================================
subroutine test_openmp_elasticity_linear_private(passed)
    use omp_lib
    use muscle_tensors
    use muscle_elasticity_linear
    implicit none

    logical, intent(out) :: passed
    integer, parameter   :: N = 100000
    integer, parameter   :: NUM_THREADS = 4

    type(ten_3D2Osym), allocatable  :: strains(:)
    type(ten_3D2Osym), allocatable  :: stress_seq(:), stress_par(:)
    type(ten_3D4O3sym), allocatable :: tan_seq(:), tan_par(:)

    type(Elasticity_linear) :: elas, elas_loc
    integer :: i
    real(8) :: t_start, t_seq, t_par

    passed = .FALSE.

    allocate(strains(N), stress_seq(N), stress_par(N), tan_seq(N), tan_par(N))

    call elas%set_parameters(young=210000.0D0, poisson=0.3D0)

    do i = 1, N
        call strains(i)%init(xx=0.001D0 * sin(dble(i)), &
                             yy=-0.0003D0 * cos(dble(i)), &
                             zz=-0.0003D0 * sin(2.0D0*dble(i)), &
                             xy=0.0005D0 * cos(0.5D0*dble(i)), &
                             yz=0.0001D0 * sin(0.1D0*dble(i)), &
                             xz=-0.0002D0 * cos(0.1D0*dble(i)))
    end do

    ! --- A. EJECUCIÓN SECUENCIAL ---
    t_start = omp_get_wtime()
    do i = 1, N
        stress_seq(i) = elas%stress(strains(i))
        tan_seq(i)    = elas%dstress_dstrain(strains(i))
    end do
    t_seq = omp_get_wtime() - t_start

    ! --- B. EJECUCIÓN PARALELA CON OBJETOS PRIVADOS POR HILO ---
    call omp_set_num_threads(NUM_THREADS)
    t_start = omp_get_wtime()

!$OMP PARALLEL PRIVATE(i, elas_loc) SHARED(strains, stress_par, tan_par)
    call elas_loc%set_parameters(young=210000.0D0, poisson=0.3D0)

    !$OMP DO
    do i = 1, N
        stress_par(i) = elas_loc%stress(strains(i))
        tan_par(i)    = elas_loc%dstress_dstrain(strains(i))
    end do
    !$OMP END DO
!$OMP END PARALLEL

    t_par = omp_get_wtime() - t_start

    write(*, '(A, F8.4, A, F8.4, A)') "   Tiempo Secuencial: ", t_seq, "s | Tiempo Paralelo: ", t_par, "s"
    write(*, '(A, F6.2, A)') "   Aceleración (Speedup): ", t_seq / t_par, "x"

    ! --- C. VERIFICACIÓN 1:1 ---
    do i = 1, N
        if (.not. (stress_par(i) .approx. stress_seq(i))) then
            print*, "   FAIL: Stress mismatch in private elasticity at point ", i
            print*, "Stress seq:"
            print*, stress_seq(i)
            print*, "Stress par:"
            print*, stress_par(i)
            return
        end if
        if (.not. (tan_par(i) .approx. tan_seq(i))) then
            print*, "   FAIL: Tangent mismatch in private elasticity at point ", i
            return
        end if
    end do

    deallocate(strains, stress_seq, stress_par, tan_seq, tan_par)
    passed = .TRUE.
end subroutine test_openmp_elasticity_linear_private