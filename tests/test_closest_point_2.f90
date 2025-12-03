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
    ! Johnson–Cook tests (4–6)
    ! =============================
    call test_closest_point_vonmises_uniaxial_tensile_jc(passed)
    if (.not. passed) STOP 4
    print*, "END - test_closest_point_vonmises_uniaxial_tensile_jc - test 4"

    call test_closest_point_vonmises_zero_strain_jc(passed)
    if (.not. passed) STOP 5
    print*, "END - test_closest_point_vonmises_zero_strain_jc - test 5"

    call test_closest_point_vonmises_elastic_strain_jc(passed)
    if (.not. passed) STOP 6
    print*, "END - test_closest_point_vonmises_elastic_strain_jc - test 6"

    print*, "All tests passed!", passed

end program test_closest_point_2

!=====================================================================
! 1) Swift: Tracción uniaxial con plastificación
!=====================================================================
subroutine test_closest_point_vonmises_uniaxial_tensile_sw(passed)
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening,      only : Swift_hardening
    use mod_vonMises,             only : VonMises
    use mod_elasticity_linear,    only : Elasticity_linear
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

    passed = stress .isequal. expected_stress
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

    passed = strain_p .isequal. expected_strain_plastic
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
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening,      only : Swift_hardening
    use mod_vonMises,             only : VonMises
    use mod_elasticity_linear,    only : Elasticity_linear
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

    passed = stress .isequal. expected_stress
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

    passed = strain_p .isequal. expected_strain_plastic
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
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_swift_hardening,      only : Swift_hardening
    use mod_vonMises,             only : VonMises
    use mod_elasticity_linear,    only : Elasticity_linear
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

    passed = stress .isequal. expected_stress
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

    passed = strain_p .isequal. expected_strain_plastic
    if (.not. passed) then
        print*, "Plastic strain is not equal", new_line('A'),  &
                "Expected:", expected_strain_plastic, new_line('A'), &
                "Actual  :", strain_p, new_line('A'),           &
                "Difference:", strain_p - expected_strain_plastic
        return
    end if

end subroutine test_closest_point_vonmises_elastic_strain_sw

!=====================================================================
! 4) Johnson–Cook: Tracción uniaxial con plastificación (check cualitativo)
!=====================================================================
subroutine test_closest_point_vonmises_uniaxial_tensile_jc(passed)
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_JohnsonCook_full_visco_hardening, only : JohnsonCook_full_visco_hardening
    use mod_vonMises,                         only : VonMises
    use mod_elasticity_linear,                only : Elasticity_linear
    use mod_closest_point_2
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(Closest_point_2_data)           :: data
    type(Closest_point_2)                :: solver
    type(VonMises)                       :: vm
    type(ten_3D2Osym)                    :: strain, strain_p, stress, stress_elastic
    type(JohnsonCook_full_visco_hardening) :: jc
    type(Elasticity_linear)              :: elas
    real(real64)                         :: strain_pf
    integer                              :: status, iters
    real(real64)                         :: seq_trial, seq_return

    ! Parámetros JC elegidos para que haya plastificación (yield bajo)
    real(real64), parameter :: A_mat      = 10.0d0
    real(real64), parameter :: B_mat      = 100.0d0
    real(real64), parameter :: n_mat      = 0.1d0
    real(real64), parameter :: C_mat      = 0.0d0
    real(real64), parameter :: epdot0_mat = 1.0d0
    real(real64), parameter :: m_mat      = 0.0d0
    real(real64), parameter :: Troom_mat  = 293.15d0
    real(real64), parameter :: Tmelt_mat  = 1500.0d0

    passed    = .False.
    strain_pf = 0.0d0

    print*, "test_closest_point_vonmises_uniaxial_tensile_jc - test 4"

    call strain%init(xx=0.43857844D0, yy=-0.20128922D0, zz=-0.20128922D0, &
                     xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)

    ! Construir JC
    jc = JohnsonCook_full_visco_hardening(A=A_mat, B=B_mat, n=n_mat, C=C_mat, &
                                          epdot0=epdot0_mat, m=m_mat,         &
                                          Troom=Troom_mat, Tmelt=Tmelt_mat)
    jc%epdot_current = epdot0_mat
    jc%T_current     = Troom_mat

    ! Inicializar solver
    call data%init(strain_pf=strain_pf, strain_p=strain_p)
    call solver%init(elasticity=elas, hardening=jc, yield=vm)

    ! Tensión elástica trial
    stress_elastic = elas%stress(strain-strain_p)
    seq_trial      = vm%stress_eq(stress_elastic)

    ! Ejecutar return mapping
    call solver%solve(strain=strain, data=data)
    call data%get(stress=stress, strain_pf=strain_pf, strain_p=strain_p, &
                  status=status, iters=iters)

    seq_return = vm%stress_eq(stress)

    ! Checks cualitativos:
    ! 1) Se desarrolló deformación plástica
    if (strain_pf <= 0.0d0) then
        print*, "JC uniaxial tensile: no plastic strain developed."
        return
    end if

    ! 2) Tensión equivalente retornada <= tensión equivalente trial (proyección a la superficie)
    if (seq_return > seq_trial + 1d-8) then
        print*, "JC uniaxial tensile: seq_return > seq_trial, algo raro."
        print*, "seq_trial =", seq_trial, " seq_return =", seq_return
        return
    end if

    ! 3) Iteraciones razonables
    if (iters > 20) then
        print*, "JC uniaxial tensile: demasiadas iteraciones:", iters
        return
    end if

    passed = .True.

end subroutine test_closest_point_vonmises_uniaxial_tensile_jc

!=====================================================================
! 5) Johnson–Cook: Deformación nula → respuesta trivial
!=====================================================================
subroutine test_closest_point_vonmises_zero_strain_jc(passed)
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_JohnsonCook_full_visco_hardening, only : JohnsonCook_full_visco_hardening
    use mod_vonMises,                         only : VonMises
    use mod_elasticity_linear,                only : Elasticity_linear
    use mod_closest_point_2,                  only : closest_point2
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(VonMises)                       :: vm
    type(ten_3D2Osym)                    :: strain, strain_p, stress
    type(JohnsonCook_full_visco_hardening) :: jc
    type(Elasticity_linear)              :: elas
    real(real64)                         :: strain_pf
    logical                              :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64)      :: expected_strain_effective

    ! Parámetros JC con yield alto
    real(real64), parameter :: A_mat      = 1000.0d0
    real(real64), parameter :: B_mat      = 0.0d0
    real(real64), parameter :: n_mat      = 1.0d0
    real(real64), parameter :: C_mat      = 0.0d0
    real(real64), parameter :: epdot0_mat = 1.0d0
    real(real64), parameter :: m_mat      = 1.0d0
    real(real64), parameter :: Troom_mat  = 293.15d0
    real(real64), parameter :: Tmelt_mat  = 1500.0d0

    passed    = .False.
    strain_pf = 0D0
    error     = .False.

    print*, "test_closest_point_vonmises_zero_strain_jc - test 5"

    expected_strain_effective = 0.0D0
    call expected_stress%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)

    call strain%init(xx=0.0D0, yy=0.0D0, zz=0.0D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)

    jc = JohnsonCook_full_visco_hardening(A=A_mat, B=B_mat, n=n_mat, C=C_mat, &
                                          epdot0=epdot0_mat, m=m_mat,         &
                                          Troom=Troom_mat, Tmelt=Tmelt_mat)
    jc%epdot_current = epdot0_mat
    jc%T_current     = Troom_mat

    call closest_point2(strain, elas, jc, vm, stress, strain_pf, strain_p, error)

    passed = stress .isequal. expected_stress
    if (.not. passed) then
        print*, "Stress is not equal (JC zero strain)", new_line('A'),    &
                "Expected:", expected_stress, new_line('A'),              &
                "Actual  :", stress, new_line('A'),                       &
                "Difference:", stress - expected_stress
        return
    end if

    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) then
        print*, "Effective plastic strain is not zero (JC zero strain)", &
                new_line('A'), "Expected:", expected_strain_effective,   &
                new_line('A'), "Actual  :", strain_pf
        return
    end if

    passed = strain_p .isequal. expected_strain_plastic
    if (.not. passed) then
        print*, "Plastic strain is not zero (JC zero strain)", new_line('A'), &
                "Expected:", expected_strain_plastic, new_line('A'),          &
                "Actual  :", strain_p
        return
    end if

end subroutine test_closest_point_vonmises_zero_strain_jc

!=====================================================================
! 6) Johnson–Cook: Deformación pequeña → respuesta puramente elástica
!=====================================================================
subroutine test_closest_point_vonmises_elastic_strain_jc(passed)
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_JohnsonCook_full_visco_hardening, only : JohnsonCook_full_visco_hardening
    use mod_vonMises,                         only : VonMises
    use mod_elasticity_linear,                only : Elasticity_linear
    use mod_closest_point_2,                  only : closest_point2
    implicit none

    real(real64), parameter :: EPS=1e-8
    logical, intent(out) :: passed

    type(VonMises)                       :: vm
    type(ten_3D2Osym)                    :: strain, strain_p, stress
    type(JohnsonCook_full_visco_hardening) :: jc
    type(Elasticity_linear)              :: elas
    real(real64)                         :: strain_pf
    logical                              :: error

    type(ten_3D2Osym) :: expected_stress, expected_strain_plastic
    real(real64)      :: expected_strain_effective

    ! Parámetros JC con yield alto (elástico)
    real(real64), parameter :: A_mat      = 1000.0d0
    real(real64), parameter :: B_mat      = 0.0d0
    real(real64), parameter :: n_mat      = 1.0d0
    real(real64), parameter :: C_mat      = 0.0d0
    real(real64), parameter :: epdot0_mat = 1.0d0
    real(real64), parameter :: m_mat      = 1.0d0
    real(real64), parameter :: Troom_mat  = 293.15d0
    real(real64), parameter :: Tmelt_mat  = 1500.0d0

    passed    = .False.
    strain_pf = 0D0
    error     = .False.

    print*, "test_closest_point_vonmises_elastic_strain_jc - test 6"

    expected_strain_effective = 0.0D0
    call expected_stress%init(xx=1D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call expected_strain_plastic%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)

    call strain%init(xx=0.001D0, yy=-0.0003D0, zz=-0.0003D0, xy=0D0, yz=0D0, xz=0D0)
    call strain_p%init(xx=0D0, yy=0D0, zz=0D0, xy=0D0, yz=0D0, xz=0D0)
    call elas%set_parameters(young=1000D0, poisson=0.3D0)

    jc = JohnsonCook_full_visco_hardening(A=A_mat, B=B_mat, n=n_mat, C=C_mat, &
                                          epdot0=epdot0_mat, m=m_mat,         &
                                          Troom=Troom_mat, Tmelt=Tmelt_mat)
    jc%epdot_current = epdot0_mat
    jc%T_current     = Troom_mat

    call closest_point2(strain, elas, jc, vm, stress, strain_pf, strain_p, error)

    passed = stress .isequal. expected_stress
    if (.not. passed) then
        print*, "Stress is not equal (JC elastic strain)", new_line('A'), &
                "Expected:", expected_stress, new_line('A'),              &
                "Actual  :", stress, new_line('A'),                       &
                "Difference:", stress - expected_stress
        return
    end if

    passed = abs(strain_pf - expected_strain_effective) < EPS
    if (.not. passed) then
        print*, "Effective plastic strain is not zero (JC elastic strain)", &
                new_line('A'), "Expected:", expected_strain_effective,      &
                new_line('A'), "Actual  :", strain_pf
        return
    end if

    passed = strain_p .isequal. expected_strain_plastic
    if (.not. passed) then
        print*, "Plastic strain is not zero (JC elastic strain)", new_line('A'), &
                "Expected:", expected_strain_plastic, new_line('A'),             &
                "Actual  :", strain_p
        return
    end if

end subroutine test_closest_point_vonmises_elastic_strain_jc
