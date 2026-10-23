! Numerical CPU reference data for future CPU/GPU comparisons.
! Uses the public MUSCLE API; does not modify or replace library algorithms.
program cpu_reference
    use iso_fortran_env, only: real64
    use ieee_arithmetic, only: ieee_is_finite
    use muscle_tensors
    use muscle_elasticity_linear, only: Elasticity_linear
    use muscle_hard_swift, only: Swift_hardening
    use muscle_yield_vonmises, only: VonMises
    use muscle_plastic_history, only: Plastic_material_history
    use muscle_solver_closest_point, only: Closest_point, STATUS_CONVERGED, STATUS_ELASTIC_CASE
    implicit none
    integer, parameter :: npoints=64, nsteps=6
    real(real64), parameter :: load_factors(nsteps)=[0.0d0,0.25d0,0.75d0,0.5d0,1.0d0,0.0d0]
    type(Elasticity_linear) :: elasticity
    type(Swift_hardening) :: hardening
    type(VonMises) :: yield
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history(npoints), before
    type(ten_3D2Osym) :: strain, flow, gradient
    type(ten_3D4O2sym) :: tangent
    real(real64) :: values(70), base(6), sign_load, dgamma, yield_residual, flow_residual, scale
    integer :: unit, ios, point, step, load, status, iterations, n_elastic, n_plastic
    character(len=1024) :: destination
    character(len=11) :: load_name

    call get_command_argument(1, destination)
    if (len_trim(destination)==0) error stop 'Expected output CSV path.'
    open(newunit=unit, file=trim(destination), status='new', action='write', iostat=ios)
    if (ios/=0) error stop 'Cannot create reference CSV (will not overwrite).'
    call elasticity%set_parameters(young=200000d0, poisson=0.3d0)
    hardening=Swift_hardening(k=100d0, n=0.1d0, e0=1d0)
    call solver%init(elasticity=elasticity, hardening=hardening, yield=yield)
    write(unit,'(a)',advance='no') 'criterion,load,point,step,status,iterations'
    call columns('strain',6)
    call columns('state_n_stress',6)
    call columns('state_n_strain_p',6)
    call columns('state_n_strain_pf',1)
    call columns('state_np1_stress',6)
    call columns('state_np1_strain_p',6)
    call columns('state_np1_strain_pf',1)
    call columns('tangent_column_major',36)
    write(unit,'(a)') ',yield_residual,flow_residual'
    n_elastic=0
    n_plastic=0
    do load=1,2
        history=Plastic_material_history()
        sign_load=merge(1d0,-1d0,load==1)
        load_name=merge('tension    ','compression',load==1)
        do step=1,nsteps
            do point=1,npoints
                scale=1d0+0.1d0*real(point,real64)/real(npoints,real64)
                base=0d0
                select case(mod(point,4))
                case(0)
                    ! Zero path, including a zero-stress elastic tangent.
                case(1)
                    base=[1d-5,-3d-6,-3d-6,0d0,0d0,0d0]
                case(2)
                    base=[0.005d0,-0.0015d0,-0.0015d0,0.0002d0,0d0,0d0]
                case(3)
                    base=[0d0,0d0,0d0,0.002d0,0.0003d0,-0.0001d0]
                end select
                call strain%init(sign_load*load_factors(step)*scale*base)
                before=history(point)
                call solver%solve(strain,history(point),status,iterations)
                if (status/=STATUS_CONVERGED .and. status/=STATUS_ELASTIC_CASE) &
                    error stop 'Reference case did not converge.'
                if (any(history(point)%state_n%stress%vals/=before%state_n%stress%vals) .or. &
                    any(history(point)%state_n%strain_p%vals/=before%state_n%strain_p%vals) .or. &
                    history(point)%state_n%strain_pf/=before%state_n%strain_pf) &
                    error stop 'solve changed the committed history.'
                call solver%tangent(strain,history(point),tangent)
                dgamma=history(point)%state_np1%strain_pf-history(point)%state_n%strain_pf
                gradient%vals=0d0
                ! The von Mises gradient is undefined at zero stress; elastic
                ! steps have dgamma=0 and do not require that derivative.
                if (status==STATUS_CONVERGED) &
                    gradient=yield%dstressEq_dstress(history(point)%state_np1%stress)
                flow=history(point)%state_n%strain_p-history(point)%state_np1%strain_p+dgamma*gradient
                yield_residual=yield%stress_eq(history(point)%state_np1%stress) &
                               -hardening%stress(history(point)%state_np1%strain_pf)
                flow_residual=sqrt(sum(flow%vals**2))
                values=[strain%vals,history(point)%state_n%stress%vals,history(point)%state_n%strain_p%vals, &
                        history(point)%state_n%strain_pf,history(point)%state_np1%stress%vals, &
                        history(point)%state_np1%strain_p%vals,history(point)%state_np1%strain_pf, &
                        reshape(tangent%vals,[36]),yield_residual,flow_residual]
                if (.not.all(ieee_is_finite(values))) error stop 'Non-finite reference value.'
                write(unit,'(a,",",a,4(",",i0),70(",",es25.16e3))') &
                    'vm',trim(load_name),point,step,status,iterations,values
                if (status==STATUS_ELASTIC_CASE) n_elastic=n_elastic+1
                if (status==STATUS_CONVERGED) n_plastic=n_plastic+1
                call history(point)%commit()
            end do
        end do
    end do
    close(unit)
    print '(a,i0,a,i0,a,i0)', 'Reference rows: ',npoints*nsteps*2, &
        '; elastic: ',n_elastic,'; converged plastic: ',n_plastic
contains
    subroutine columns(prefix,n)
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: n
        integer :: i
        do i=1,n
            write(unit,'(",",a,"_",i0)',advance='no') prefix,i
        end do
    end subroutine columns
end program cpu_reference
