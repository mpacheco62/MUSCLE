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

    type(ten_3D2Osym), allocatable :: a(:), b(:)
    type(ten_3D2Osym), allocatable :: dev_seq(:), dev_par(:)
    real(8), allocatable           :: ddot_seq(:), ddot_par(:)

    type(ten_3D2Osym):: tmp_dev

    integer :: i
    real(8) :: t_start, t_seq, t_par, hydro

    passed = .FALSE.

    allocate(a(N), b(N), dev_seq(N), dev_par(N), ddot_seq(N), ddot_par(N))

    do i = 1, N
        call a(i)%init(xx=10.0D0 * dble(i)/N, yy=5.0D0 * dble(i)/N, zz=2.0D0 * dble(i)/N, &
                       xy=1.0D0 * dble(i)/N, yz=0.5D0 * dble(i)/N, xz=-0.2D0 * dble(i)/N)
        call b(i)%init(xx=1.0D0, yy=-0.5D0, zz=-0.5D0, xy=0.2D0, yz=0.1D0, xz=0.0D0)
    end do

    ! --- A. SECUENCIAL ---
    t_start = omp_get_wtime()
    do i = 1, N
        dev_seq(i)  = .dev. a(i)
        ddot_seq(i) = a(i) .ddot. b(i)
    end do
    t_seq = omp_get_wtime() - t_start

    ! --- B. PARALELO ---
    call omp_set_num_threads(NUM_THREADS)
    t_start = omp_get_wtime()

!$OMP PARALLEL DO PRIVATE(i, tmp_dev) SHARED(a, b, dev_par, ddot_par)
    do i = 1, N
        dev_par(i)  = .dev. a(i)
        ddot_par(i) = a(i) .ddot. b(i)
    end do
!$OMP END PARALLEL DO

    t_par = omp_get_wtime() - t_start

    write(*, '(A, F8.4, A, F8.4, A)') "   Tiempo Secuencial: ", t_seq, "s | Tiempo Paralelo: ", t_par, "s"

    ! --- C. VERIFICACIÓN ---
    do i = 1, N
        if (.not. (dev_par(i) .approx. dev_seq(i))) then
            print*, "   FAIL: Tensor .dev. mismatch at index ", i
            print*, "   a(i)   :", a(i)%vals
            print*, "   seq    :", dev_seq(i)%vals
            print*, "   par    :", dev_par(i)%vals
            print*, "   diff    :", dev_seq(i)%vals - dev_par(i)%vals
            return
        end if

        if (abs(ddot_par(i) - ddot_seq(i)) > 1.0D-10) then
            print*, "   FAIL: Tensor .ddot. mismatch at index ", i
            print*, "   seq    :", ddot_seq(i)
            print*, "   par    :", ddot_par(i)
            return
        end if
    end do

    deallocate(a, b, dev_seq, dev_par, ddot_seq, ddot_par)
    passed = .TRUE.
end subroutine test_openmp_tensor_ops