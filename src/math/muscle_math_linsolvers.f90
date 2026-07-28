! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_math_linsolvers
  implicit none
  private
  public :: gauss_solve, gauss_seidel
  public :: mklsolve

contains

  !=========================================================
  !  gauss_solve: Eliminación de Gauss con pivoteo parcial
  !
  !  Resuelve A x = b (A densa) de forma directa.
  !  Modifica A y b internamente, así que normalmente
  !  conviene pasar copias si quieres conservarlos.
  !=========================================================
  subroutine gauss_solve(A, b, x, info)
    implicit none
    real(8), intent(inout) :: A(:,:)   ! Matriz del sistema (se sobreescribe)
    real(8), intent(inout) :: b(:)     ! RHS (se sobreescribe)
    real(8), intent(out)   :: x(:)     ! Solución
    integer, intent(out)   :: info     ! =0 OK, ≠0 hay problema

    integer :: n, k, i, j, pivot_row
    real(8) :: maxval_col, factor
    real(8), allocatable :: Awork(:,:), bwork(:)
    real(8) :: tmp_r

    n = size(A,1)
    if (size(A,2) /= n .or. size(b,1) /= n .or. size(x,1) /= n) then
       info = -1
       return
    end if

    ! Trabajamos con copias para no destruir los originales
    allocate(Awork(n,n), bwork(n))
    Awork = A
    bwork = b

    info = 0

    ! Fase de eliminación
    do k = 1, n-1
       ! Pivoteo parcial: buscar fila con máximo valor absoluto
       maxval_col = 0.0d0
       pivot_row  = k
       do i = k, n
          if (abs(Awork(i,k)) > maxval_col) then
             maxval_col = abs(Awork(i,k))
             pivot_row  = i
          end if
       end do

       ! Si el pivote es ~0, matriz casi singular
       if (maxval_col == 0.0d0) then
          info = 1
          deallocate(Awork, bwork)
          return
       end if

       ! Intercambiar filas si hace falta
       if (pivot_row /= k) then
          do j = 1, n
             tmp_r = Awork(k,j)
             Awork(k,j) = Awork(pivot_row,j)
             Awork(pivot_row,j) = tmp_r
          end do

          tmp_r = bwork(k)
          bwork(k) = bwork(pivot_row)
          bwork(pivot_row) = tmp_r
       end if

       ! Eliminación hacia adelante
       do i = k+1, n
          factor = Awork(i,k) / Awork(k,k)
          Awork(i,k) = 0.0d0
          do j = k+1, n
             Awork(i,j) = Awork(i,j) - factor * Awork(k,j)
          end do
          bwork(i) = bwork(i) - factor * bwork(k)
       end do
    end do

    ! Comprobar último pivote
    if (Awork(n,n) == 0.0d0) then
       info = 1
       deallocate(Awork, bwork)
       return
    end if

    ! Sustitución hacia atrás
    x(n) = bwork(n) / Awork(n,n)
    do i = n-1, 1, -1
       tmp_r = bwork(i)
       do j = i+1, n
          tmp_r = tmp_r - Awork(i,j) * x(j)
       end do
       x(i) = tmp_r / Awork(i,i)
    end do

    deallocate(Awork, bwork)
  end subroutine gauss_solve

  !=========================================================
  !  gauss_seidel: método iterativo de Gauss–Seidel
  !
  !  Útil para jugar con convergencia, no es la opción
  !  "eficiente" para un 7x7 denso, pero sirve didácticamente.
  !=========================================================
  subroutine gauss_seidel(A, b, x, tol, maxit, its, info)
    implicit none
    real(8), intent(in)    :: A(:,:), b(:)
    real(8), intent(inout) :: x(:)        ! x inicial y luego solución
    real(8), intent(in)    :: tol         ! tolerancia (p.ej. 1.0d-10)
    integer, intent(in)    :: maxit       ! máx. iteraciones
    integer, intent(out)   :: its         ! iteraciones realizadas
    integer, intent(out)   :: info        ! 0=OK, 1=no converge

    integer :: n, i, j, k
    real(8) :: sigma, res_norm

    n = size(A,1)
    if (size(A,2) /= n .or. size(b,1) /= n .or. size(x,1) /= n) then
       info = -1
       its  = 0
       return
    end if

    info = 1
    its  = 0

    do k = 1, maxit
       ! Un barrido de Gauss-Seidel
       do i = 1, n
          sigma = 0.0d0
          do j = 1, n
             if (j /= i) sigma = sigma + A(i,j)*x(j)
          end do
		  if (A(i,i) == 0.0d0) then
			info = 2   ! código especial: diagonal nula
			return
		  end if

          x(i) = (b(i) - sigma) / A(i,i)
       end do

       ! Norma del residuo (simple) para criterio de parada
       res_norm = 0.0d0
       do i = 1, n
          sigma = 0.0d0
          do j = 1, n
             sigma = sigma + A(i,j)*x(j)
          end do
          res_norm = res_norm + (sigma - b(i))**2
       end do
       res_norm = sqrt(res_norm)

       its = k
       if (res_norm < tol) then
          info = 0
          exit
       end if
    end do

  end subroutine gauss_seidel


!===========================================================
  ! test_mkl:
  !   Misma interfaz que gauss_seidel:
  !     A(:,:)  -> matriz del sistema (no se modifica)
  !     b(:)    -> RHS original (no se modifica)
  !     x(:)    -> solución (salida)
  !     tol     -> no usado (solo para compatibilidad)
  !     maxit   -> no usado (solo para compatibilidad)
  !     its     -> número "de iteraciones" (lo ponemos a 1)
  !     info    -> código de salida de dgesv
  !
  !   Internamente:
  !     - Hace copias locales de A y b (A_loc, b_loc)
  !     - Llama a dgesv (MKL/LAPACK)
  !     - Devuelve x = solución
  !===========================================================
  subroutine mklsolve(A, b, x, tol, maxit, its, info)
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    real(8), intent(in)    :: A(:,:), b(:)
    real(8), intent(out)   :: x(:)
    real(8), intent(in)    :: tol
    integer, intent(in)    :: maxit
    integer, intent(out)   :: its, info

    integer :: n, lda, ldb, nrhs
    real(8), allocatable :: A_loc(:,:), b_loc(:,:)
    integer, allocatable :: ipiv(:)
    integer :: i

    ! Interfaz explícita de dgesv
    interface
       subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         integer :: n, nrhs, lda, ldb, info
         integer :: ipiv(*)
         real(8) :: a(lda,*), b(ldb,*)
       end subroutine dgesv

#ifdef USE_MKL
       subroutine mkl_set_num_threads_local(nth) bind(C, name="mkl_set_num_threads_local")
         import :: c_int
         integer(c_int), value :: nth
       end subroutine mkl_set_num_threads_local
#endif
    end interface

#ifdef USE_MKL
    ! Solo se ejecuta si estamos compilando/enlazando con Intel MKL
    call mkl_set_num_threads_local(1_c_int)
#endif

    n = size(A,1)
    if (size(A,2) /= n .or. size(b) /= n .or. size(x) /= n) then
       info = -1
       its  = 0
       return
    end if

    lda  = n
    ldb  = n
    nrhs = 1

    allocate(A_loc(n,n), b_loc(n, nrhs), ipiv(n))

    ! Copiamos A y b a los arreglos de trabajo (dgesv los sobreescribe)
    A_loc = A
    do i = 1, n
       b_loc(i,1) = b(i)
    end do

    call dgesv(n, nrhs, A_loc, lda, ipiv, b_loc, ldb, info)

    if (info == 0) then
       do i = 1, n
          x(i) = b_loc(i,1)
       end do
       its = 1   ! "una iteración" ficticia (llamada directa)
    else
       its = 0
    end if

    deallocate(A_loc, b_loc, ipiv)
  end subroutine mklsolve

end module muscle_math_linsolvers
