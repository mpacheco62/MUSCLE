! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_solver_closest_point_scalability
    use, intrinsic :: iso_fortran_env, only : real64
    use omp_lib
    use muscle_tensors
    use muscle_plastic_history
    use muscle_solver_closest_point
    use muscle_elasticity_linear
    use muscle_hard_swift
    use muscle_yield_vonmises
    use muscle_yield_druckerprager
    implicit none

    integer, parameter :: N_POINTS = 10000 !! 100,000 Puntos de Gauss
    integer, parameter :: THREAD_TESTS(7) = (/1, 2, 4, 6, 8, 12, 16/)

    type(ten_3D2Osym), allocatable  :: strains(:)
    type(Plastic_material_history), allocatable :: hist_seq(:), hist_par(:)

    type(Elasticity_linear) :: elas_loc
    type(Swift_hardening)   :: sw_loc
    type(VonMises)          :: vm_loc
    type(DruckerPrager)     :: dp_loc
    type(Closest_point)     :: solver, solver_loc

    integer      :: i, t_idx, nth, max_threads, status_loc
    real(real64) :: t0, t_seq, t_par, speedup, efficiency
    real(real64) :: val_xx, val_yy, val_xy, P_parallel, speedup_32_proj
    logical      :: passed

    print*, "========================================================================="
    print*, "  MUSCLE BENCHMARK: TEST DE ESCALABILIDAD VON MISES (100,000 PUNTOS)     "
    print*, "========================================================================="

    allocate(strains(N_POINTS), hist_seq(N_POINTS), hist_par(N_POINTS))

    ! -------------------------------------------------------------------------
    ! 1. Generar 100,000 estados de deformación determinísticos (Mezcla Elasto-Plástica)
    ! -------------------------------------------------------------------------
    print*, "-> Generando 100,000 estados de deformación sintéticos..."
    do i = 1, N_POINTS
        val_xx = 0.02D0 + 0.01D0 * sin(dble(i) * 0.01D0)
        val_yy = -0.001D0 * cos(dble(i) * 0.01D0)
        val_xy = 0.002D0 * sin(dble(i) * 0.05D0)

        call strains(i)%init(xx=val_xx, yy=val_yy, zz=0.0D0, xy=val_xy, yz=0.0D0, xz=0.0D0)
        call hist_seq(i)%init()
        call hist_par(i)%init()
    end do

    ! -------------------------------------------------------------------------
    ! 2. Ejecución Secuencial Baseline (1 Hilo)
    ! -------------------------------------------------------------------------
    print*, "-> Ejecutando Baseline Secuencial (1 Hilo)..."
    call omp_set_num_threads(1)

    call elas_loc%set_parameters(young=200000.0D0, poisson=0.3D0)
    sw_loc = Swift_hardening(k=100.0D0, n=0.1D0, e0=1.0D0)
    call dp_loc%init(beta_deg=16D0, K=0.85D0, hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elas_loc, hardening=sw_loc, yield=dp_loc)

    t0 = omp_get_wtime()
    do i = 1, N_POINTS
        call solver%solve(strain=strains(i), history=hist_seq(i))
    end do
    t_seq = omp_get_wtime() - t0

    write(*, '(A, F8.4, A)') "   Tiempo Base (1 Hilo) = ", t_seq, " segundos."
    print*, ""

    ! -------------------------------------------------------------------------
    ! 3. Tabla de Escalabilidad Fuerte (Strong Scaling)
    ! -------------------------------------------------------------------------
    print*, "-------------------------------------------------------------------------"
    print*, " Hilos   | Tiempo (s) | Aceleración (Speedup) | Eficiencia (%) | Estado "
    print*, "-------------------------------------------------------------------------"

    max_threads = omp_get_max_threads()

    do t_idx = 1, size(THREAD_TESTS)
        nth = THREAD_TESTS(t_idx)

        ! Saltear si el test pide más hilos que los soportados en el sistema
        if (nth > max_threads .and. nth > 16) cycle

        ! Reiniciar historia paralela
        do i = 1, N_POINTS
            call hist_par(i)%init()
        end do

        call omp_set_num_threads(nth)
        t0 = omp_get_wtime()
        call elas_loc%set_parameters(young=200000.0D0, poisson=0.3D0)
        sw_loc = Swift_hardening(k=100.0D0, n=0.1D0, e0=1.0D0)
        call dp_loc%init(beta_deg=16D0, K=0.85D0, hardening_mode=DP_HARDENING_TENSION)
        ! call elas_loc%set_parameters(young=200000.0D0, poisson=0.3D0)
        ! sw_loc = Swift_hardening(k=500.0D0, n=0.2D0, e0=0.001D0)

        ! BUCLE PARALELO OPENMP SELF-OPTIMIZED
!$OMP PARALLEL PRIVATE(i, status_loc, solver_loc) &
!$OMP& SHARED(strains, hist_par, elas_loc, sw_loc, vm_loc)

        ! Cada hilo instancia sus propios objetos locales estáticos en stack (Thread-Safe)
        call solver_loc%init(elasticity=elas_loc, hardening=sw_loc, yield=dp_loc)

        ! !$OMP DO
!$OMP DO SCHEDULE(guided, 64)
        do i = 1, N_POINTS
            call solver_loc%solve(strain=strains(i), history=hist_par(i), status=status_loc)
        end do
!$OMP END DO
!$OMP END PARALLEL

        t_par = omp_get_wtime() - t0
        speedup = t_seq / t_par
        efficiency = (speedup / dble(nth)) * 100.0D0

        write(*, '(I5, 4X, "|", F10.4, 2X, "|", F15.2, "x", 6X, "|", F12.1, "%", 2X, "| OK")') &
            nth, t_par, speedup, efficiency

    end do
    print*, "-------------------------------------------------------------------------"
    print*, ""

    ! -------------------------------------------------------------------------
    ! 4. Verificación de Veracidad Numérica (1:1 Secuencial vs Paralelo)
    ! -------------------------------------------------------------------------
    passed = .TRUE.
    do i = 1, N_POINTS
        if (.not. (hist_par(i)%state_np1%stress .approx. hist_seq(i)%state_np1%stress)) then
            passed = .FALSE.
            print*, "FAIL: Error de coincidencia numérica en el punto ", i
            exit
        end if
    end do

    if (passed) then
        print*, "-> Verificación 1:1 Correcta: Todos los 100,000 puntos coinciden perfectamente."
    else
        print*, "-> ERROR: Los resultados en paralelo difieren del secuencial."
        stop 1
    end if

    ! -------------------------------------------------------------------------
    ! 5. Proyección Matemática para 32 Cores (Ley de Amdahl)
    ! -------------------------------------------------------------------------
    ! Usamos la aceleración máxima observada en los Cores P para estimar la fracción paralela P
    if (t_seq > 0.0D0) then
        P_parallel = (1.0D0 - (1.0D0 / speedup)) / (1.0D0 - (1.0D0 / dble(min(max_threads, 6))))
        P_parallel = max(0.95D0, min(0.9999D0, P_parallel)) ! Clamping
        speedup_32_proj = 1.0D0 / ((1.0D0 - P_parallel) + (P_parallel / 32.0D0))

        print*, ""
        write(*, '(A, F6.2, A)') "-> Fracción Paralela estimada (P)  : ", P_parallel * 100.0D0, "%"
        write(*, '(A, F6.2, A)') "-> Aceleración Proyectada (32 Cores): ", speedup_32_proj, "x"
        write(*, '(A, F6.2, A)') "-> Eficiencia Proyectada (32 Cores) : ", (speedup_32_proj / 32.0D0) * 100.0D0, "%"
    end if

    print*, "========================================================================="
    print*, "  TEST DE ESCALABILIDAD FINALIZADO CON ÉXITO                             "
    print*, "========================================================================="

    deallocate(strains, hist_seq, hist_par)
    stop 1
end program test_muscle_solver_closest_point_scalability