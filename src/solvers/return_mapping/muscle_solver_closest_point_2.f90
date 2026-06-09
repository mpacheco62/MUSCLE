! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco@usach.cl>

module muscle_solver_closest_point_2
	use muscle_math_linsolvers ! con mkl
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_base, only : Base_hardening_laws

    use muscle_yield_base, only : Base_yield_critera
    use muscle_elasticity_base, only : Base_elasticity

    private

    integer, parameter, public :: STATUS_UNDEFINED=-1
    integer, parameter, public :: STATUS_ELASTIC_CASE=0
    integer, parameter, public :: STATUS_CONVERGED=1
    integer, parameter, public :: STATUS_NONCONVERGED=2
    integer, parameter, public :: STATUS_ITER_CONVERGED=10
    integer, parameter, public :: STATUS_ITER_NONCONVERGED=11

    type, public :: Closest_point_2_data
        real(real64), private :: strain_pf_init=0D0, strain_pf_iter=0D0
        type(ten_3D2Osym), private:: strain_p_init = ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))
        type(ten_3D2Osym), private :: strain_p_iter=ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))
        real(real64), private :: dgamma=0D0, ddgamma=0D0
        type(ten_3D2Osym), private :: stress = ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))
        type(integer), private :: status=STATUS_UNDEFINED, iters=0
        contains
        procedure, public :: init => closest_point_data_init
        procedure, public :: set => closest_point_data_set
        procedure, public :: get => closest_point_data_get
    end type Closest_point_2_data

    type, public :: Closest_point_2
        class(Base_elasticity), allocatable :: elasticity
        class(Base_hardening_laws), allocatable :: hardening
        class(Base_yield_critera), allocatable :: yield

        integer :: iter_nw=200
        
        contains
        procedure, public :: init => closest_point_init
        procedure, private :: iter => closest_point_iter
        procedure, public :: solve => closest_point_solve
    end type Closest_point_2



    public :: closest_point2

    contains
        subroutine closest_point_data_init(self, strain_pf, strain_p)
            implicit none
            class(Closest_point_2_data), intent(inout) :: self
            type(ten_3D2Osym), intent(in) :: strain_p
            real(real64), intent(in) :: strain_pf
            self%strain_pf_init = strain_pf
            self%strain_p_init = strain_p 

            self%strain_pf_iter = strain_pf
            self%strain_p_iter = strain_p

            self%stress%vals = 0D0
            self%dgamma = 0D0
            self%ddgamma = 0D0
            self%status = STATUS_UNDEFINED
            self%iters = 0
        end subroutine closest_point_data_init

        subroutine closest_point_data_set(self,                            &
                                          strain_pf_init, strain_pf_iter,  &
                                          strain_p_init, strain_p_iter,    &
                                          dgamma, ddgamma,                 &
                                          stress                           &
                                          )
            implicit none
            class(Closest_point_2_data), intent(inout) :: self
            real(real64), intent(in), optional :: strain_pf_init, strain_pf_iter
            type(ten_3D2Osym), intent(in), optional :: strain_p_init, strain_p_iter
            real(real64), intent(in), optional :: dgamma, ddgamma
            type(ten_3D2Osym), intent(in), optional :: stress


            if(present(strain_pf_init)) self%strain_pf_init=strain_pf_init
            if(present(strain_pf_iter)) self%strain_pf_iter=strain_pf_iter
            if(present(strain_p_init)) self%strain_p_init=strain_p_init
            if(present(strain_p_iter)) self%strain_p_iter=strain_p_iter
            if(present(dgamma)) self%dgamma=dgamma
            if(present(ddgamma)) self%ddgamma=ddgamma
            if(present(stress)) self%stress=stress
        end subroutine closest_point_data_set

        subroutine closest_point_data_get(self,                            &
                                          strain_pf_init, strain_pf,  &
                                          strain_p_init, strain_p,    &
                                          dgamma, ddgamma,                 &
                                          stress,                          &
                                          status, iters                     &
                                          )
            implicit none
            class(Closest_point_2_data), intent(in) :: self
            real(real64), intent(out), optional :: strain_pf_init, strain_pf
            type(ten_3D2Osym), intent(out), optional :: strain_p_init, strain_p
            real(real64), intent(out), optional :: dgamma, ddgamma
            type(ten_3D2Osym), intent(out), optional :: stress
            integer, intent(out), optional :: status, iters

            if(present(strain_pf_init)) strain_pf_init=self%strain_pf_init
            if(present(strain_pf)) strain_pf=self%strain_pf_iter
            if(present(strain_p_init)) strain_p_init=self%strain_p_init
            if(present(strain_p)) strain_p=self%strain_p_iter
            if(present(dgamma)) dgamma=self%dgamma
            if(present(ddgamma)) ddgamma=self%ddgamma
            if(present(stress)) stress=self%stress
            if(present(status)) status=self%status
            if(present(iters)) iters=self%iters
        end subroutine closest_point_data_get

        subroutine closest_point_init(self, elasticity, hardening, yield, iter_nw)
            implicit none
            class(Closest_point_2), intent(inout) :: self
            class(Base_elasticity), intent(in) :: elasticity
            class(Base_hardening_laws), intent(in) :: hardening
            class(Base_yield_critera), intent(in) :: yield
            integer, optional, intent(in) :: iter_nw
            self%elasticity = elasticity
            self%hardening = hardening
            self%yield = yield

            self%iter_nw = 200
            if (present(iter_nw)) self%iter_nw = iter_nw
        end subroutine closest_point_init
		
		
        subroutine closest_point_iter(self, strain, data)
            implicit none
            class(Closest_point_2), intent(in) :: self
            type(ten_3D2Osym), intent(in) :: strain
            type(Closest_point_2_data), intent(inout) :: data

            real(real64), parameter :: TOL=1D-5

            ! closest_point_data
            real(real64) :: dgamma, ddgamma, strain_pf, strain_pf_init
            type(ten_3D2Osym) :: strain_p, strain_p_init
			
			!---
			type(ten_3D2Osym) :: e_basis, tmpC, tmpXi   ! <--- NUEVOS !A11 [TODO] A:B para 4to orden.
			type(ten_3D2Osym) :: tmpA21 !A21 [TODO] sobrecargar operadores %vals da problemas
			!--- 
			
            type(ten_3D2Osym) :: df, residual1, dstrain_p, stress
            type(ten_3D4O3sym) :: elas_tan, EEinv, ddf
            real(real64) :: hard, dhard, f, norm_res
            real(real64) :: residual(7)
            ! type(ten_3D4O2sym) :: test
			
			!! MKL - PARAMETERS
			integer, parameter :: n = 7
			real(real64) :: AA(6,6)
			real(real64) :: A(n,n), b(n), x(n)
			real(real64) :: tolmkl
			integer :: maxit, its, info
			
            dgamma = data%dgamma
            ddgamma = data%ddgamma
            strain_pf = data%strain_pf_iter
            strain_pf_init = data%strain_pf_init
            strain_p = data%strain_p_iter
            strain_p_init = data%strain_p_init
            data%iters = data%iters + 1

		    !C
            elas_tan = self%elasticity%dstress_dstrain(strain-strain_p)  ! constant
			
			! print*, "!*** Start algorithm ***"			
            !*** Start algorithm ***
            stress = self%elasticity%stress(strain-strain_p)
            hard = self%hardening%stress(strain_pf)
			!! (Sigma_eq - Y) debe ser menor o igual 0
			!! Rf_n+1		
            f = self%yield%stress_eq(stress) - hard
			
			!! DY/Desp_p
			! A22
            dhard = self%hardening%dstress_dep(strain_pf)

            !! DSigma_eq/DSigma
			df = self%yield%dstressEq_dstress(stress)				
			!! D2sigma_eq/D2Sigma
			ddf = self%yield%ddstressEq_ddstress(stress)
			
			! Rep_n+1
            residual1 = strain_p_init - strain_p + (dgamma*df)
            
			!! CONVERGENCY CONDITION
            residual(1:6) = residual1%vals**2
            norm_res = (residual(1)+residual(2)+residual(3)+residual(4)+residual(5)+residual(6))**0.5					
            if ((abs(f) .lt. tol) .and. (norm_res .lt. tol)) then  ! Converged
                data%status = STATUS_ITER_CONVERGED
                data%stress = stress
                return
            end if
				
			EEinv = ((.inv. elas_tan) + dgamma*ddf)
			
			! =========================
			! A11: bloque 6x6
			! A11 * Δε^p = - EEinv : (elas_tan : Δε^p)
			! =========================
			!print*,"A11"
			!print*, - ( EEinv .ddot. elas_tan)
			!A(1:6,1:6) = - tmpA11%vals	 ! ESTO NO FUNCIONA
			
			do its = 1, 6
				! 1) Vector base en Voigt: [0,0,0,0,0,0]^T salvo 1 en la componente k
				e_basis = ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))
				e_basis%vals(its) = 1.0_real64

				! 2) tmpC = C : e_basis  (4º con 2º -> 2º)
				tmpC  = elas_tan .ddot. e_basis
	
				! 3) tmpXi = EEinv : tmpC  (4º con 2º -> 2º)
				tmpXi = EEinv .ddot. tmpC

				! 4) Columna k de A11 = - tmpXi
				A(1:6, its) = - tmpXi%vals
			end do
			!!!! END A11
			
			!A12
			A(1:6,7) = df%vals

			!A21
			tmpA21= ( -elas_tan .ddot. df)
			A(7,1:6) = tmpA21%vals

			!A22
			A(7,7) = -dhard

			!b	
			b(1:6) = - residual1%vals
			b(7)   = - f
			
			!!!! SOLVE LINEAL SYSTEM W/ MKL
			tolmkl   = 1.0d-12
			maxit = 1000

			call mklsolve(A, b, x, tolmkl, maxit, its, info)
			
			if (info == 0) then
				! print *, "Solucion x:"
				! print '(3(f12.6,1x))', x
				! print *, "its =", its
				
				! 1) Pasar x al formato tensorial
				dstrain_p = ten_3D2Osym(x(1:6))   ! constructor de tu tipo
				ddgamma   = x(7)

				! 2) Actualizar variables internas de esta iteración
				dgamma   = dgamma + ddgamma              ! Δγ_{k+1} = Δγ_k + Δ²γ
				strain_p = strain_p + dstrain_p          ! ε^p_{k+1} = ε^p_k + Δε^p
				strain_pf = strain_pf_init + dgamma      ! o strain_pf = strain_pf + ddgamma

				! 3) Guardar en la estructura data
				data%dgamma         = dgamma
				data%ddgamma        = ddgamma
				data%strain_pf_iter = strain_pf
				data%strain_p_iter  = strain_p
				
			else
				print *, "Error en test_mkl, info =", info
			end if
			!!! END - PRUEBAS MKL			
			
            data%status = STATUS_ITER_NONCONVERGED  !  not converged
        end subroutine closest_point_iter

        subroutine closest_point_solve(self, strain, data)
            implicit none
            class(Closest_point_2), intent(in) :: self
            type(ten_3D2Osym), intent(in) :: strain
            type(Closest_point_2_data), intent(inout) :: data

            real(real64), parameter :: TOL2=1D-16
            type(ten_3D2Osym) :: stress, strain_p
            real(real64) :: strain_pf, hard, f
            real(real64) :: stress_eq

            integer :: i
			
			!!! SEPARAR EL ESTADO ANTERIOR 
            strain_p = data%strain_p_init
            strain_pf = data%strain_pf_init

            !*** Check Elastic Case ***
            stress = self%elasticity%stress(strain-strain_p)
            hard = self%hardening%stress(strain_pf)
            stress_eq = self%yield%stress_eq(stress)
            f = stress_eq - hard
            
			!! CASO ELASTICO
			if (stress_eq/hard - 1D0 .le. -TOL2) then  ! Elastic Case			
				! print*,"Elastic Case"
                data%status = STATUS_ELASTIC_CASE				
                data%stress = stress
				return
            end if
			
			!!! ACA COMIENZA LA ITERACION DE NEWTON RAPHSON
            do i=1,self%iter_nw 
				! print*,i
                call self%iter(strain, data)
                if (data%status .eq. STATUS_ITER_CONVERGED) then  ! converged
                    data%status = STATUS_CONVERGED
                    return
                end if
            end do
            data%status = STATUS_NONCONVERGED  ! not converged
        end subroutine closest_point_solve

        !subroutine closest_point2(strain, elasticity, hardening, yield, strain_pf, strain_p, status, stress)
		subroutine closest_point2(strain, elasticity, hardening, yield, stress, strain_pf, strain_p, status)
            implicit none
            class(ten_3D2Osym), intent(in) :: strain
            class(Base_elasticity), intent(in) :: elasticity
            class(Base_hardening_laws), intent(in) :: hardening
            class(Base_yield_critera), intent(in) :: yield
            real(real64), intent(inout) :: strain_pf
            type(ten_3D2Osym), intent(inout) :: strain_p
            type(ten_3D2Osym), intent(out) :: stress
            logical, intent(out) :: status

            type(Closest_point_2_data) :: data
            type(Closest_point_2) :: solver
						
            call solver%init(elasticity=elasticity, hardening=hardening, yield=yield)
            call data%init(strain_pf=strain_pf, strain_p=strain_p)
            call solver%solve(strain, data)
            call data%get(stress = stress, strain_pf=strain_pf, strain_p=strain_p)
            status = .false.

        end subroutine closest_point2


end module