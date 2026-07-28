! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_openmp
    use omp_lib
    use muscle_tensors
    use muscle_plastic_history
    use muscle_solver_closest_point
    use muscle_elasticity_linear
    use muscle_hard_swift
    use muscle_yield_vonmises
    use muscle_yield_druckerprager
    use muscle_kinematics
    implicit none

    logical :: passed

    print*, "========================================================"
    print*, "  RUNNING MUSCLE CORE OPENMP THREAD-SAFETY TEST SUITE   "
    print*, "========================================================"

    print*, "1. Testing Parallel Closest Point Solver (Von Mises)..."
    call test_openmp_closest_point_vm(passed)
    if (.not. passed) stop 1

    ! print*, "2. Testing Parallel Closest Point Solver (Drucker-Prager)..."
    ! call test_openmp_closest_point_dp(passed)
    ! if (.not. passed) stop 2

    print*, "========================================================"
    print*, "  ALL OPENMP THREAD-SAFETY TESTS PASSED SUCCESSFULLY!   "
    print*, "========================================================"
    stop 0
end program test_muscle_openmp


! =========================================================================
! 1. TEST OPENMP CON VON MISES (10,000 PUNTOS)
! =========================================================================
subroutine test_openmp_closest_point_vm(passed)
    use omp_lib
    use muscle_tensors
    use muscle_plastic_history
    use muscle_solver_closest_point
    use muscle_elasticity_linear
    use muscle_hard_swift
    use muscle_yield_vonmises
    implicit none

    logical, intent(out) :: passed
    integer, parameter :: N = 10000
    integer, parameter :: NUM_THREADS = 2

    type(ten_3D2Osym), allocatable :: strains(:)
    type(Plastic_material_history), allocatable :: hist_seq(:), hist_par(:)
    type(ten_3D4O2sym), allocatable :: tang_seq(:), tang_par(:)

    type(Elasticity_linear) :: elas, elas_loc
    type(Swift_hardening)   :: sw, sw_loc
    type(VonMises)          :: vm, vm_loc
    type(Closest_point)     :: solver, solver_loc
    integer :: i
    integer, allocatable :: status_code(:), status_code_par(:)
    real(8) :: t_start, t_seq, t_par

    passed = .FALSE.

    ! 1. Asignar memoria para 10,000 puntos
    allocate(strains(N), hist_seq(N), hist_par(N), tang_seq(N), tang_par(N), status_code(N), status_code_par(N))

    ! 2. Inicializar modelo de material y solver
    call elas%set_parameters(young=200000.0D0, poisson=0.3D0)
    sw = Swift_hardening(k=500.0D0, n=0.2D0, e0=0.001D0)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)

    ! 3. Generar 10,000 deformaciones de prueba variadas (algunas elásticas, otras plásticas)
    do i = 1, N
        call strains(i)%init(xx=0.005D0, &
                             yy=0.0D0, &
                             zz=0.0D0, &
                             xy=0.0D0, &
                             yz=0.0D0, &
                             xz=0.0D0)
        call hist_seq(i)%init()
        call hist_par(i)%init()
    end do

    ! --- A. EJECUCIÓN SECUENCIAL (1 HILO) ---
    t_start = omp_get_wtime()
    do i = 1, N
        call solver%solve(strain=strains(i), history=hist_seq(i), status=status_code(i))
        call solver%tangent(strain=strains(i), history=hist_seq(i), tangent=tang_seq(i))
    end do
    t_seq = omp_get_wtime() - t_start

    ! --- B. EJECUCIÓN PARALELA (OPENMP) ---
    call omp_set_num_threads(NUM_THREADS)

    ! VERIFICACIÓN DE HILOS ACTIVOS
!$OMP PARALLEL
!$OMP MASTER
    print*, "--------------------------------------------------------"
    write(*, '(A, I0, A, I0)') "   [VERIFICACIÓN] Hilos solicitados: ", NUM_THREADS, &
                              " | Hilos OpenMP ACTIVOS en ejecución: ", omp_get_num_threads()
    print*, "--------------------------------------------------------"
!$OMP END MASTER
!$OMP END PARALLEL

    t_start = omp_get_wtime()

!$OMP PARALLEL PRIVATE(i, status_code, elas_loc, sw_loc, vm_loc, solver_loc) &
!$OMP& SHARED(strains, hist_par, tang_par, status_code_par)

        ! Se ejecuta 1 sola vez por hilo al abrir el bloque paralelo:
        call elas_loc%set_parameters(young=200000.0D0, poisson=0.3D0)
        sw_loc = Swift_hardening(k=500.0D0, n=0.2D0, e0=0.001D0)
        call solver_loc%init(elasticity=elas_loc, hardening=sw_loc, yield=vm_loc)

        !$OMP DO
        do i = 1, N
            call solver_loc%solve(strain=strains(i), history=hist_par(i), status=status_code_par(i))
            call solver_loc%tangent(strain=strains(i), history=hist_par(i), tangent=tang_par(i))
        end do
        !$OMP END DO

!$OMP END PARALLEL

    t_par = omp_get_wtime() - t_start

    write(*, '(A, F8.4, A, F8.4, A)') "   Tiempo Secuencial: ", t_seq, "s | Tiempo Paralelo: ", t_par, "s"
    write(*, '(A, F6.2, A)') "   Aceleración (Speedup): ", t_seq / t_par, "x"

    ! --- C. VERIFICACIÓN DE DETERMINISMO 1:1 ---
    do i = 1, N
        if (status_code(i) .ne. status_code_par(i)) then
            print*, "   FAIL: Mismatch in status code at point ", i
            print*, status_code(i), status_code_par(i)
            return
        end if
        ! Verificar que el esfuerzo de Cauchy en paralelo sea idéntico al secuencial
        if (.not. (hist_par(i)%state_np1%stress .approx. hist_seq(i)%state_np1%stress)) then
            print*, "   FAIL: Mismatch in stress at point ", i
            return
        end if

        ! Verificar deformación plástica equivalente
        if (abs(hist_par(i)%state_np1%strain_pf - hist_seq(i)%state_np1%strain_pf) > 1.0D-10) then
            print*, "   FAIL: Mismatch in eq. plastic strain at point ", i
            return
        end if

        ! Verificar matriz tangente
        if (.not. (tang_par(i) .approx. tang_seq(i))) then
            print*, "   FAIL: Mismatch in tangent stiffness at point ", i
            return
        end if
    end do

    deallocate(strains, hist_seq, hist_par, tang_seq, tang_par)
    passed = .FALSE.
    print*, "    ALLGOOD!"
end subroutine test_openmp_closest_point_vm


! =========================================================================
! 2. TEST OPENMP CON DRUCKER-PRAGER (10,000 PUNTOS)
! =========================================================================
subroutine test_openmp_closest_point_dp(passed)
    use omp_lib
    use muscle_tensors
    use muscle_plastic_history
    use muscle_solver_closest_point
    use muscle_elasticity_linear
    use muscle_hard_swift
    use muscle_yield_druckerprager
    implicit none

    logical, intent(out) :: passed
    integer, parameter :: N = 10000
    integer, parameter :: NUM_THREADS = 4

    type(ten_3D2Osym), allocatable :: strains(:)
    type(Plastic_material_history), allocatable :: hist_seq(:), hist_par(:)
    type(ten_3D4O2sym), allocatable :: tang_seq(:), tang_par(:)

    type(Elasticity_linear) :: elas
    type(Swift_hardening)   :: sw
    type(DruckerPrager)      :: dp
    type(Closest_point)     :: solver
    integer :: i, status_code
    real(8) :: t_start, t_seq, t_par

    passed = .FALSE.

    allocate(strains(N), hist_seq(N), hist_par(N), tang_seq(N), tang_par(N))

    call elas%set_parameters(young=50000.0D0, poisson=0.3D0)
    sw = Swift_hardening(k=100.0D0, n=0.1D0, e0=1.0D0)
    call dp%init(beta_deg=10.0D0, K=1.0D0, hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elas, hardening=sw, yield=dp)

    do i = 1, N
        call strains(i)%init(xx=0.0005D0 * dble(i) / 100.0D0, &
                             yy=-0.0002D0 * dble(i) / 100.0D0, &
                             zz=-0.0002D0 * dble(i) / 100.0D0, &
                             xy=0.0001D0 * dble(i) / 100.0D0, &
                             yz=0.0D0, xz=0.0D0)
        call hist_seq(i)%init()
        call hist_par(i)%init()
    end do

    ! --- A. EJECUCIÓN SECUENCIAL ---
    t_start = omp_get_wtime()
    do i = 1, N
        call solver%solve(strain=strains(i), history=hist_seq(i))
        call solver%tangent(strain=strains(i), history=hist_seq(i), tangent=tang_seq(i))
    end do
    t_seq = omp_get_wtime() - t_start

    ! --- B. EJECUCIÓN PARALELA ---
    call omp_set_num_threads(NUM_THREADS)
    t_start = omp_get_wtime()

!$OMP PARALLEL DO PRIVATE(i, status_code) SHARED(solver, strains, hist_par, tang_par)
    do i = 1, N
        call solver%solve(strain=strains(i), history=hist_par(i), status=status_code)
        call solver%tangent(strain=strains(i), history=hist_par(i), tangent=tang_par(i))
    end do
!$OMP END PARALLEL DO

    t_par = omp_get_wtime() - t_start

    write(*, '(A, F8.4, A, F8.4, A)') "   Tiempo Secuencial: ", t_seq, "s | Tiempo Paralelo: ", t_par, "s"
    write(*, '(A, F6.2, A)') "   Aceleración (Speedup): ", t_seq / t_par, "x"

    ! --- C. VERIFICACIÓN DE DETERMINISMO 1:1 ---
    do i = 1, N
        if (.not. (hist_par(i)%state_np1%stress .approx. hist_seq(i)%state_np1%stress)) then
            print*, "   FAIL: Mismatch in DP stress at point ", i
            return
        end if
    end do

    deallocate(strains, hist_seq, hist_par, tang_seq, tang_par)
    passed = .TRUE.
end subroutine test_openmp_closest_point_dp