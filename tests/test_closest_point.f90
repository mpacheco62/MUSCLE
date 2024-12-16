program test_closest_point
    implicit none
    logical :: passed

    call test_closest_point_vonmises(passed)
    if (.not. passed) STOP 1


    print*, "Passed!", passed
    STOP 0
end program test_closest_point

subroutine test_closest_point_vonmises(passed)
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening, only : Swift_hardening
    use mod_vonMises, only : VonMises
    use mod_elasticity_linear, only : Elasticity_linear
    use mod_closest_point, only : closest_point
    implicit none

    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed

    type(VonMises) :: vm
    type(ten_3D2Osym) :: strain, strain_p, stress
    type(Swift_hardening) :: sw
    type(Elasticity_linear) :: elas
    real(real64) :: strain_pf
    logical :: error

    passed = .False.
    strain_pf = 0D0
    error = .False.
    

    call strain%init(xx=0.15D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call closest_point(strain=strain, elasticity=elas, hardening=sw, yield=vm,    &
                       stress=stress, strain_pf=strain_pf, strain_p=strain_p, error=error)

    write(*,*) "stress"
    write(*,*) stress
    write(*,*) 

    write(*,*) "strain_pf:", strain_pf
    write(*,*) 

    write(*,*) "strain_p"
    write(*,*) strain_p
    write(*,*) 

    write(*,*) "Error:", error
    ! result = vm%stress_eq(to_test1)
    ! passed = (abs(result - expected_result1) < EPS)
    ! if (.not. passed) return
end subroutine
