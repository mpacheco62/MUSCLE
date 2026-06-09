program test_closest_point_2
    implicit none
    logical :: passed

    ! =============================
    ! Swift hardening tests (1–3)
    ! =============================
    call test_closest_point_vonmises_uniaxial_tensile_sw(passed)
    if (.not. passed) STOP 1
    print*, "END - test_closest_point_vonmises_uniaxial_tensile_sw - test 1"

    call test_closest_point_vonmises_zero_strain_sw(passed)
    if (.not. passed) STOP 2
    print*, "END - test_closest_point_vonmises_zero_strain_sw - test 2"

    call test_closest_point_vonmises_elastic_strain_sw(passed)
    if (.not. passed) STOP 3
    print*, "END - test_closest_point_vonmises_elastic_strain_sw - test 3"

    ! =============================
    ! Johnson–Cook (viscoplastic) tests (4–6)
    ! =============================
    call test_closest_point_vonmises_uniaxial_tensile_jc(passed)
    if (.not. passed) STOP 4
    print*, "END - test_closest_point_vonmises_uniaxial_tensile_jc - test 4"

    call test_closest_point_vonmises_zero_strain_jc(passed)
    if (.not. passed) STOP 5
    print*, "END - test_closest_point_vonmises_zero_strain_jc - test 5"

    ! call test_closest_point_vonmises_elastic_strain_jc(passed)
    ! if (.not. passed) STOP 6
    ! print*, "END - test_closest_point_vonmises_elastic_strain_jc - test 6"

    print*, "All tests passed!", passed

end program test_closest_point_2

!=====================================================================
! 1) Swift: Tracción uniaxial con plastificación
!=====================================================================
subroutine test_closest_point_vonmises_uniaxial_tensile_sw(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening,      only : Swift_hardening
    use mod_vonMises,             only : VonMises
    use muscle_elasticity_linear,    only : Elasticity_linear
    use mod_closest_point_2
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed
    type(Closest_point_2_data) :: data
    type(Closest_point_2)      :: solver
    type(VonMises)             :: vm
    type(ten_3D2Osym)          :: strain, strain_p, stress
    type(Swift_hardening)      :: sw
    type(Elasticity_linear)    :: elas
    real(real64)               :: strain_pf
    logical                    :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64)      :: expected_strain_effective
    integer           :: status, iters

    passed    = .False.
    strain_pf = 0.0d0
    error     = .False.

    print*, "test_closest_point_vonmises_uniaxial_tensile_sw - test 1"

    call expected_stress%init(xx=90D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    expected_strain_effective = 0.3485784401D0
    call expected_strain_plastic%init(xx=expected_strain_effective,        &
                                      yy=-0.5D0*expected_strain_effective, &
                                      zz=-0.5D0*expected_strain_effective, &
                                      xy=0D0, yz=0D0, xz=0D0)

    call strain%init(xx=0.43857844D0, yy=-0.20128922D0, zz=-0.20128922D0, &
                     xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call data%init(strain_pf=strain_pf, strain_p=strain_p)
    call solver%init(elasticity=elas, hardening=sw, yield=vm)
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress,        &
                  strain_pf=strain_pf,  &
                  strain_p=strain_p,    &
                  status=status,        &
                  iters=iters)

    passed = stress%is_approx(expected_stress, tol=1D-5)
    if (.not. passed) then
        print*, "Stress is not equal", new_line('A'),          &
                "Expected:", expected_stress, new_line('A'),   &
                "Actual  :", stress, new_line('A'),            &
                "Difference:", stress - expected_stress
        return
    end if

    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) then
        print*, "Effective plastic strain is not equal", new_line('A'), &
                "Expected:", expected_strain_effective, new_line('A'),  &
                "Actual  :", strain_pf, new_line('A'),                  &
                "Difference:", strain_pf - expected_strain_effective
        return
    end if

    passed = strain_p%is_approx(expected_strain_plastic, tol=1D-5)
    if (.not. passed) then
        print*, "Plastic strain is not equal", new_line('A'),  &
                "Expected:", expected_strain_plastic, new_line('A'), &
                "Actual  :", strain_p, new_line('A'),           &
                "Difference:", strain_p - expected_strain_plastic
        return
    end if

    passed = iters .le. 6
    if (.not. passed) then
        print*, "Iterations are greater than expected", new_line('A'), &
                "Iters:", iters, " Maximum:", 6, new_line('A')
        return
    end if

end subroutine test_closest_point_vonmises_uniaxial_tensile_sw

!=====================================================================
! 2) Swift: Deformación nula → respuesta elástica trivial
!=====================================================================
subroutine test_closest_point_vonmises_zero_strain_sw(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening,      only : Swift_hardening
    use mod_vonMises,             only : VonMises
    use muscle_elasticity_linear,    only : Elasticity_linear
    use mod_closest_point_2,      only : closest_point2
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(VonMises)          :: vm
    type(ten_3D2Osym)       :: strain, strain_p, stress
    type(Swift_hardening)   :: sw
    type(Elasticity_linear) :: elas
    real(real64)            :: strain_pf
    logical                 :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64)      :: expected_strain_effective

    passed    = .False.
    strain_pf = 0D0
    error     = .False.

    print*, "test_closest_point_vonmises_zero_strain_sw - test 2"

    expected_strain_effective = 0.0D0
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)

    call strain%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call closest_point2(strain, elas, sw, vm, stress, strain_pf, strain_p, error)

    passed = stress .approx. expected_stress
    if (.not. passed) then
        print*, "Stress is not equal", new_line('A'),          &
                "Expected:", expected_stress, new_line('A'),   &
                "Actual  :", stress, new_line('A'),            &
                "Difference:", stress - expected_stress
        return
    end if

    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) then
        print*, "Effective plastic strain is not equal", new_line('A'), &
                "Expected:", expected_strain_effective, new_line('A'),  &
                "Actual  :", strain_pf, new_line('A'),                  &
                "Difference:", strain_pf - expected_strain_effective
        return
    end if

    passed = strain_p .approx. expected_strain_plastic
    if (.not. passed) then
        print*, "Plastic strain is not equal", new_line('A'),  &
                "Expected:", expected_strain_plastic, new_line('A'), &
                "Actual  :", strain_p, new_line('A'),           &
                "Difference:", strain_p - expected_strain_plastic
        return
    end if

end subroutine test_closest_point_vonmises_zero_strain_sw

!=====================================================================
! 3) Swift: Deformación pequeña → respuesta puramente elástica
!=====================================================================
subroutine test_closest_point_vonmises_elastic_strain_sw(passed)
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening,      only : Swift_hardening
    use mod_vonMises,             only : VonMises
    use muscle_elasticity_linear,    only : Elasticity_linear
    use mod_closest_point_2,      only : closest_point2
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(VonMises)          :: vm
    type(ten_3D2Osym)       :: strain, strain_p, stress
    type(Swift_hardening)   :: sw
    type(Elasticity_linear) :: elas
    real(real64)            :: strain_pf
    logical                 :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64)      :: expected_strain_effective

    passed    = .False.
    strain_pf = 0D0
    error     = .False.

    print*, "test_closest_point_vonmises_elastic_strain_sw - test 3"

    expected_strain_effective = 0.0D0
    call expected_stress%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)

    call strain%init(xx=0.001D0, yy=-0.0003D0, zz=-0.0003D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)

    call closest_point2(strain, elas, sw, vm, stress, strain_pf, strain_p, error)

    passed = stress .approx. expected_stress
    if (.not. passed) then
        print*, "Stress is not equal", new_line('A'),          &
                "Expected:", expected_stress, new_line('A'),   &
                "Actual  :", stress, new_line('A'),            &
                "Difference:", stress - expected_stress
        return
    end if

    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) then
        print*, "Effective plastic strain is not equal", new_line('A'), &
                "Expected:", expected_strain_effective, new_line('A'),  &
                "Actual  :", strain_pf, new_line('A'),                  &
                "Difference:", strain_pf - expected_strain_effective
        return
    end if

    passed = strain_p .approx. expected_strain_plastic
    if (.not. passed) then
        print*, "Plastic strain is not equal", new_line('A'),  &
                "Expected:", expected_strain_plastic, new_line('A'), &
                "Actual  :", strain_p, new_line('A'),           &
                "Difference:", strain_p - expected_strain_plastic
        return
    end if

end subroutine test_closest_point_vonmises_elastic_strain_sw

!=====================================================================
! 4) Johnson–Cook (viscoplastic): tracción uniaxial (chequeo de magnitud)
!=====================================================================
subroutine test_closest_point_vonmises_uniaxial_tensile_jc(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening,   only : Swift_hardening
    use mod_JC_viscoplastic,   only : JC_viscoplastic
    implicit none

    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed
    type(Swift_hardening), target :: sw
    type(JC_viscoplastic) :: jc
    real(real64) :: ep, epd
    real(real64) :: expected, actual

    passed = .False.

    print*, "test_closest_point_vonmises_uniaxial_tensile_jc - test 4"

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    jc = JC_viscoplastic(C=0.05D0, epdmax=0.1D0)
    jc%hard_law = sw  ! Error in gfortran-12 if not assigned here

    ep  = 0.2D0
    epd = 1.0D0

    expected = sw%stress(ep) * (1.0D0 + jc%C * log(epd / jc%epdmax))
    actual   = jc%flow_stress(ep, epd)

    passed = abs(actual - expected) < EPS
    if (.not. passed) then
        print*, "JC flow stress mismatch", new_line('A'), &
                "Expected:", expected, new_line('A'),     &
                "Actual  :", actual, new_line('A'),       &
                "Difference:", actual - expected
        return
    end if

end subroutine test_closest_point_vonmises_uniaxial_tensile_jc

!=====================================================================
! 5) Johnson–Cook (viscoplastic): deformación nula con tasa de referencia
!=====================================================================
subroutine test_closest_point_vonmises_zero_strain_jc(passed)
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening,   only : Swift_hardening
    use mod_JC_viscoplastic,   only : JC_viscoplastic
    implicit none

    real(real64), parameter :: EPS=1e-10
    logical, intent(out) :: passed
    type(Swift_hardening), target :: sw
    type(JC_viscoplastic) :: jc
    real(real64) :: ep, epd
    real(real64) :: expected, actual

    passed = .False.

    print*, "test_closest_point_vonmises_zero_strain_jc - test 5"

    sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
    jc = JC_viscoplastic(C=0.05, epdmax=0.1)
    jc%hard_law = sw  ! Error in gfortran-12 if not assigned here

    ep  = 0.0D0
    epd = jc%epdmax

    expected = sw%stress(ep)
    actual   = jc%flow_stress(ep, epd)

    passed = abs(actual - expected) < EPS
    if (.not. passed) then
        print*, "JC flow stress mismatch at zero strain", new_line('A'), &
                "Expected:", expected, new_line('A'),                   &
                "Actual  :", actual, new_line('A'),                     &
                "Difference:", actual - expected
        return
    end if

end subroutine test_closest_point_vonmises_zero_strain_jc

! !=====================================================================
! ! 6) Johnson–Cook (viscoplastic): tasa nula → respuesta cero ! NO DEBERÍA SER CERO la tasa
! !=====================================================================
! subroutine test_closest_point_vonmises_elastic_strain_jc(passed)
!     use, intrinsic :: iso_fortran_env, only : real64
!     use mod_swift_hardening,   only : Swift_hardening
!     use mod_JC_viscoplastic,   only : JC_viscoplastic
!     implicit none

!     real(real64), parameter :: EPS=1e-12
!     logical, intent(out) :: passed
!     type(Swift_hardening), target :: sw
!     type(JC_viscoplastic) :: jc
!     real(real64) :: ep, epd
!     real(real64) :: actual

!     passed = .False.

!     print*, "test_closest_point_vonmises_elastic_strain_jc - test 6"

!     sw = Swift_hardening(k=100D0, n=0.1D0, e0=1D-4)
!     jc%hard_law => sw
!     jc%C      = 0.05D0
!     jc%epdmax = 0.1D0

!     ep  = 0.001D0
!     epd = 0.0D0

!     actual = jc%flow_stress(ep, epd)
!     passed = abs(actual) < EPS
!     if (.not. passed) then
!         print*, "JC flow stress not zero for null rate", new_line('A'), &
!                 "Actual  :", actual
!         return
!     end if

! end subroutine test_closest_point_vonmises_elastic_strain_jc
