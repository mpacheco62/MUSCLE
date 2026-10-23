! Minimal reproduction of the nonconvergent upstream benchmark, first point.
program reproduce_dp_tension
    use iso_fortran_env, only: real64
    use muscle_tensors
    use muscle_elasticity_linear, only: Elasticity_linear
    use muscle_hard_swift, only: Swift_hardening
    use muscle_yield_druckerprager, only: DruckerPrager, DP_HARDENING_TENSION
    use muscle_plastic_history, only: Plastic_material_history
    use muscle_solver_closest_point, only: Closest_point, STATUS_NONCONVERGED
    implicit none
    type(Elasticity_linear) :: elasticity
    type(Swift_hardening) :: hardening
    type(DruckerPrager) :: yield
    type(Closest_point) :: solver
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain
    real(real64) :: residual
    integer :: status, iterations
    call elasticity%set_parameters(young=200000d0,poisson=0.3d0)
    hardening=Swift_hardening(k=100d0,n=0.1d0,e0=1d0)
    call yield%init(beta_deg=16d0,K=0.85d0,hardening_mode=DP_HARDENING_TENSION)
    call solver%init(elasticity=elasticity,hardening=hardening,yield=yield)
    call strain%init(xx=0.02d0+0.01d0*sin(0.01d0),yy=-0.001d0*cos(0.01d0), &
                     zz=0d0,xy=0.002d0*sin(0.05d0),yz=0d0,xz=0d0)
    call history%init()
    call solver%solve(strain,history,status,iterations)
    residual=yield%stress_eq(history%state_np1%stress)-hardening%stress(history%state_np1%strain_pf)
    print '(a,i0)', 'status=',status
    print '(a,i0)', 'reported_iterations=',iterations
    print '(a,6es25.16e3)', 'strain=',strain%vals
    print '(a,6es25.16e3)', 'stress=',history%state_np1%stress%vals
    print '(a,es25.16e3)', 'yield_residual=',residual
    print '(a,es25.16e3)', 'equivalent_plastic_strain=',history%state_np1%strain_pf
    if (status/=STATUS_NONCONVERGED) error stop 'Observed behavior differs from recorded known issue.'
    print '(a)', 'KNOWN_NONCONVERGENCE_REPRODUCED; this is not a physically validated result.'
end program reproduce_dp_tension
