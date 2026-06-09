module mod_closest_point_time
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_viscoplastic_laws, only : Base_viscoplastic_law
    use mod_yield_criteria, only : Base_yield_critera
    use mod_base_elasticity, only : Base_elasticity

    private

    integer, parameter, public :: STATUS_UNDEFINED=-1
    integer, parameter, public :: STATUS_ELASTIC_CASE=0
    integer, parameter, public :: STATUS_CONVERGED=1
    integer, parameter, public :: STATUS_NONCONVERGED=2
    integer, parameter, public :: STATUS_ITER_CONVERGED=10
    integer, parameter, public :: STATUS_ITER_NONCONVERGED=11

    type, public :: Closest_point_data
        real(real64), private :: strain_pf_init=0D0, strain_pf_iter=0D0
        type(ten_3D2Osym), private:: strain_p_init = ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))
        type(ten_3D2Osym), private :: strain_p_iter=ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))
        real(real64), private :: dlambda=0D0, ddlambda=0D0
        type(ten_3D2Osym), private :: stress = ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))
        type(integer), private :: status=STATUS_UNDEFINED, iters=0
        real(real64), private :: dt=0D0
        real(real64), private :: omega=1D0
        contains
        procedure, public :: init => closest_point_data_init
        procedure, public :: set => closest_point_data_set
        procedure, public :: get => closest_point_data_get
    end type Closest_point_data

    type, public :: Closest_point
        class(Base_elasticity), allocatable :: elasticity
        class(Base_viscoplastic_law), allocatable :: hardening
        class(Base_yield_critera), allocatable :: yield
        integer, private :: iter_nw = 200
        
        contains
        procedure, public :: init => closest_point_init
        procedure, private :: iter => closest_point_iter
        procedure, public :: solve => closest_point_solve
        procedure, public :: tangent => closest_point_tangent
        procedure, public :: tangent_numerical => closest_point_tangent_numerical
    end type Closest_point

    public :: closest_point2

    contains
        pure subroutine closest_point_data_init(self, strain_pf, strain_p, dt)
            implicit none
            class(Closest_point_data), intent(inout) :: self
            type(ten_3D2Osym), intent(in) :: strain_p
            real(real64), intent(in) :: strain_pf
            real(real64), intent(in) :: dt
            self%strain_pf_init = strain_pf
            self%strain_p_init = strain_p 
            self%dt = dt

            self%strain_pf_iter = strain_pf
            self%strain_p_iter = strain_p

            self%stress%vals = 0D0
            self%dlambda = 0D0
            self%ddlambda = 0D0
            self%status = STATUS_UNDEFINED
            self%iters = 0
        end subroutine closest_point_data_init

        pure subroutine closest_point_data_set(self,                            &
                                               strain_pf_init, strain_pf_iter,  &
                                               strain_p_init, strain_p_iter,    &
                                               dlambda, ddlambda,                 &
                                               stress,                          &
                                               dt                               &
                                               )
            implicit none
            class(Closest_point_data), intent(inout) :: self
            real(real64), intent(in), optional :: strain_pf_init, strain_pf_iter
            type(ten_3D2Osym), intent(in), optional :: strain_p_init, strain_p_iter
            real(real64), intent(in), optional :: dlambda, ddlambda
            type(ten_3D2Osym), intent(in), optional :: stress
            real(real64), intent(in), optional :: dt


            if(present(strain_pf_init)) self%strain_pf_init=strain_pf_init
            if(present(strain_pf_iter)) self%strain_pf_iter=strain_pf_iter
            if(present(strain_p_init)) self%strain_p_init=strain_p_init
            if(present(strain_p_iter)) self%strain_p_iter=strain_p_iter
            if(present(dlambda)) self%dlambda=dlambda
            if(present(ddlambda)) self%ddlambda=ddlambda
            if(present(stress)) self%stress=stress
            if(present(dt)) self%dt=dt
        end subroutine closest_point_data_set

        pure subroutine closest_point_data_get(self,                            &
                                               strain_pf_init, strain_pf,       &
                                               strain_p_init, strain_p,         &
                                               dlambda, ddlambda,                 &
                                               stress,                          &
                                               status, iters,                   &
                                               dt                               &
                                               )
            implicit none
            class(Closest_point_data), intent(in) :: self
            real(real64), intent(out), optional :: strain_pf_init, strain_pf
            type(ten_3D2Osym), intent(out), optional :: strain_p_init, strain_p
            real(real64), intent(out), optional :: dlambda, ddlambda
            type(ten_3D2Osym), intent(out), optional :: stress
            integer, intent(out), optional :: status, iters
            real(real64), intent(out), optional :: dt

            if(present(strain_pf_init)) strain_pf_init=self%strain_pf_init
            if(present(strain_pf)) strain_pf=self%strain_pf_iter
            if(present(strain_p_init)) strain_p_init=self%strain_p_init
            if(present(strain_p)) strain_p=self%strain_p_iter
            if(present(dlambda)) dlambda=self%dlambda
            if(present(ddlambda)) ddlambda=self%ddlambda
            if(present(stress)) stress=self%stress
            if(present(status)) status=self%status
            if(present(iters)) iters=self%iters
            if(present(dt)) dt=self%dt
        end subroutine closest_point_data_get

        pure subroutine closest_point_init(self, elasticity, hardening, yield, iter_nw)
            implicit none
            class(Closest_point), intent(inout) :: self
            class(Base_elasticity), intent(in) :: elasticity
            class(Base_viscoplastic_law), intent(in) :: hardening
            class(Base_yield_critera), intent(in) :: yield
            integer, optional, intent(in) :: iter_nw
            self%elasticity = elasticity
            self%hardening = hardening
            self%yield = yield
            self%iter_nw = 200
            if (present(iter_nw)) self%iter_nw = iter_nw
        end subroutine closest_point_init

        pure subroutine closest_point_iter(self, strain, data)
            implicit none
            class(Closest_point), intent(in) :: self
            type(ten_3D2Osym), intent(in) :: strain
            type(Closest_point_data), intent(inout) :: data

            real(real64), parameter :: TOL=1D-5

            ! closest_point_data
            real(real64) :: dlambda, ddlambda, strain_pf, strain_pf_init
            type(ten_3D2Osym) :: strain_p, strain_p_init

            type(ten_3D2Osym) :: df, residual1, dstrain_p, stress
            type(ten_3D4O3sym) :: elas_tan, hess, ddf
            real(real64) :: hard, dhard, f, norm_res, dt
            real(real64) :: residual(7)
            ! type(ten_3D4O2sym) :: test
            real(real64) :: strain_pf_dot

            dlambda = data%dlambda
            ddlambda = data%ddlambda
            strain_pf = data%strain_pf_iter
            strain_pf_init = data%strain_pf_init
            strain_p = data%strain_p_iter
            strain_p_init = data%strain_p_init
            data%iters = data%iters + 1
            dt = data%dt

            strain_pf_dot = (strain_pf - strain_pf_init)/dt

            elas_tan = self%elasticity%dstress_dstrain(strain-strain_p)  ! constant

            !*** Start algorithm ***
            stress = self%elasticity%stress(strain-strain_p)
            hard = self%hardening%flow_stress(strain_pf, strain_pf_dot)
            f = self%yield%stress_eq(stress) - hard

            dhard = self%hardening%dstress_dep(strain_pf, strain_pf_dot, dt)
            df = self%yield%dstressEq_dstress(stress)

            residual1 = strain_p_init - strain_p + (dlambda*df)
            ! residual2 = strain_pf_init - strain_pf + dlambda

            residual(1:6) = residual1%vals**2
            norm_res = (residual(1)+residual(2)+residual(3)+residual(4)+residual(5)+residual(6))**0.5

            if ((abs(f) .lt. tol) .and. (norm_res .lt. tol)) then  ! Converged
                data%status = STATUS_ITER_CONVERGED
                data%stress = stress
                return
            end if

            ddf = self%yield%ddstressEq_ddstress(stress)
            hess = .inv. ((.inv. elas_tan) + dlambda*ddf)

            ! Creo que residual2 siempre es cero en estos casos
            !ddlambda = (f-(df .ddot. hess .ddot. residual1) - dhard*residual2)/((df .ddot. hess .ddot. df) + dhard)
            ddlambda = (f-(df .ddot. hess .ddot. residual1))/((df .ddot. hess .ddot. df) + dhard)
            dlambda = dlambda + ddlambda*data%omega
            if (dlambda < 0D0) then 
                dlambda = 0D0  ! Ensure non-negative plastic multiplier
                data%omega = data%omega * 0.75D0  ! Reduce relaxation factor if we hit the boundary
            end if
            if (data%iters .eq. 10) then
                data%omega = data%omega * 0.75D0  ! Reduce relaxation factor if we have many iterations
            else if (data%iters .eq. 20) then
                data%omega = data%omega * 0.75D0  ! Further reduce relaxation factor if we still have many iterations
            else if (data%iters .eq. 30) then
                data%omega = data%omega * 0.75D0  ! Further reduce relaxation factor if we still have many iterations
            else if (data%iters .eq. 40) then
                data%omega = data%omega * 0.75D0  ! Further reduce relaxation factor if we still have many iterations
            end if
            strain_pf = strain_pf_init + dlambda
                ! strain_p = strain_p_init + dlambda*df
            dstrain_p = (((.inv. elas_tan) .ddot. hess) .ddot. (residual1 + ddlambda*df))
            strain_p = strain_p + dstrain_p
            data%dlambda = dlambda
            data%ddlambda = ddlambda
            data%strain_pf_iter = strain_pf
            data%strain_p_iter = strain_p
            data%status = STATUS_ITER_NONCONVERGED  !  not converged
        end subroutine closest_point_iter

        pure subroutine closest_point_solve(self, strain, data)
            implicit none
            class(Closest_point), intent(in) :: self
            type(ten_3D2Osym), intent(in) :: strain
            type(Closest_point_data), intent(inout) :: data

            real(real64), parameter :: TOL2=1D-8
            type(ten_3D2Osym) :: stress, strain_p
            real(real64) :: strain_pf, hard, f, dt
            real(real64) :: strain_pf_dot
            real(real64) :: stress_eq

            integer :: i


            strain_p = data%strain_p_init
            strain_pf = data%strain_pf_init
            dt = data%dt

            strain_pf_dot = (strain_pf - data%strain_pf_init)/dt

            !*** Check Elastic Case ***
            stress = self%elasticity%stress(strain-strain_p)
            hard = self%hardening%flow_stress(strain_pf, strain_pf_dot)
            stress_eq = self%yield%stress_eq(stress)
            f = stress_eq - hard
            
            if (stress_eq/hard - 1D0 .le. TOL2) then  ! Elastic Case
                data%status = STATUS_ELASTIC_CASE
                data%stress = stress
                return
            end if

            do i=1,self%iter_nw
                call self%iter(strain, data)
                if (data%status .eq. STATUS_ITER_CONVERGED) then  ! converged
                    data%status = STATUS_CONVERGED
                    return
                end if
            end do
            data%status = STATUS_NONCONVERGED  ! not converged
        end subroutine closest_point_solve



        pure subroutine closest_point_tangent(self, strain, data, tangent)
            implicit none
            class(Closest_point), intent(in) :: self
            type(ten_3D2Osym), intent(in) :: strain
            type(Closest_point_data), intent(inout) :: data
            type(ten_3D4O2sym), intent(out) :: tangent

            ! closest_point_data
            real(real64) :: dlambda, strain_pf
            type(ten_3D2Osym) :: strain_p

            type(ten_3D2Osym) :: df, stress
            type(ten_3D4O3sym) :: elas_tan, hess, ddf
            real(real64) :: hard, dhard, f
            real(real64) :: strain_pf_dot, dt
            ! type(ten_3D4O2sym) :: test

            if (data%status .eq. STATUS_ELASTIC_CASE) then  ! Elastic Case
                tangent = self%elasticity%dstress_dstrain(strain-strain_p)
                return
            end if


            dlambda = data%dlambda
            strain_pf = data%strain_pf_iter
            strain_p = data%strain_p_iter
            dt = data%dt

            strain_pf_dot = (strain_pf - data%strain_pf_init)/dt
            elas_tan = self%elasticity%dstress_dstrain(strain-strain_p)  ! constant

            !*** Start algorithm ***
            stress = self%elasticity%stress(strain-strain_p)
            hard = self%hardening%flow_stress(strain_pf, strain_pf_dot)
            f = self%yield%stress_eq(stress) - hard

            dhard = self%hardening%dstress_dep(strain_pf, strain_pf_dot, dt)
            df = self%yield%dstressEq_dstress(stress)

            ddf = self%yield%ddstressEq_ddstress(stress)
            hess = .inv. ((.inv. elas_tan) + dlambda*ddf)

            tangent = hess - (.tdotsym. (df .ddot. hess))/((df .ddot. hess .ddot. df) + dhard)
        end subroutine closest_point_tangent

        subroutine closest_point_tangent_numerical(self, strain, data, tangent)
            use muscle_math_derivatives
            implicit none
            class(Closest_point), intent(in) :: self
            type(ten_3D2Osym), intent(in) :: strain
            type(Closest_point_data), intent(in) :: data
            type(ten_3D4O2sym), intent(out) :: tangent


            tangent = derivative(wrapper, strain)
            
            contains
            pure function wrapper(x)
                implicit none
                type(ten_3D2Osym), intent(in) :: x
                type(ten_3D2Osym) :: wrapper
                type(Closest_point_data) :: data_local
                
                data_local = data
                data_local%status = STATUS_UNDEFINED

                call self%solve(x, data_local)
                wrapper = data_local%stress
            end function wrapper

        end subroutine closest_point_tangent_numerical   


end module