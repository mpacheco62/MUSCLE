! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module test_perzyna_joint_tangent_helpers
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_vonmises, only : VonMises
    implicit none
    private
    public :: Joint_vonmises

    type, extends(VonMises) :: Joint_vonmises
    contains
        procedure :: stress_eq => shifted_stress
        procedure :: stressEq_and_derivatives => joint_vonmises_evaluation
        procedure, nopass :: joint_evaluation => joint_true
    end type Joint_vonmises

contains

    pure function shifted_stress(self, stress) result(res)
        ! Deliberately different separate route, used to detect an incorrect dispatch.
        implicit none
        class(Joint_vonmises), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        real(real64) :: res
        res = self%VonMises%stress_eq(stress) + 20D0
    end function shifted_stress

    pure subroutine joint_vonmises_evaluation(self, stress, f, df, ddf)
        ! Unshifted parent evaluation used by the joint tangent test below.
        implicit none
        class(Joint_vonmises), intent(in) :: self
        type(ten_3D2Osym), intent(in) :: stress
        real(real64), intent(out) :: f
        type(ten_3D2Osym), intent(out) :: df
        type(ten_3D4O3sym), intent(out) :: ddf
        f = self%VonMises%stress_eq(stress)
        df = self%VonMises%dstressEq_dstress(stress)
        ddf = self%VonMises%ddstressEq_ddstress(stress)
    end subroutine joint_vonmises_evaluation

    pure function joint_true() result(res)
        ! Request the joint binding in the tangent dispatch test.
        implicit none
        logical :: res
        res = .true.
    end function joint_true
end module test_perzyna_joint_tangent_helpers

program test_muscle_solver_closest_point_perzyna_joint_tangent
    implicit none
    logical :: passed

    call test_perzyna_tangent_joint_route_mixed(passed)
    if (.not. passed) STOP 1
    call test_perzyna_numerical_tangent_elastic_mixed(passed)
    if (.not. passed) STOP 2

    print*, "Passed!", passed
    STOP 0
end program test_muscle_solver_closest_point_perzyna_joint_tangent

subroutine test_perzyna_tangent_joint_route_mixed(passed)
    ! The separate route differs on purpose; all normal and shear components are coupled.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_vonmises, only : VonMises
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point, STATUS_CONVERGED
    use muscle_test_utils, only : check_approx, check_true
    use test_perzyna_joint_tangent_helpers, only : Joint_vonmises
    implicit none
    logical, intent(out) :: passed
    ! Same parent criterion operations. Measured: 0 relative difference.
    real(real64), parameter :: TOL_ROUTE = 0D0
    ! Central differences of the converged stress. Measured: 2.79e-9 relative.
    real(real64), parameter :: TOL_TANGENT = 1D-7
    type(Closest_point) :: solver, reference
    type(Perzyna_overstress) :: vp
    type(Elasticity_linear) :: elas
    type(VonMises) :: vm
    type(Joint_vonmises) :: joint
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain
    type(ten_3D4O2sym) :: tangent, expected, numerical
    integer :: status

    passed = .true.
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    vp%hard_law = Bilinear_hardening(y0=50D0, k=100D0)
    vp%eta = 0.02D0
    vp%n = 2D0
    call strain%init(xx=0.09D0, yy=-0.025D0, zz=-0.015D0, xy=0.025D0, yz=-0.015D0, xz=0.01D0)
    call reference%init(elas, vp, vm, tol=1D-12)
    call solver%init(elas, vp, joint, tol=1D-12)
    call history%init()
    call reference%solve(strain, history, 0.1D0, status)
    call check_true(passed, "Reference plastic step converged", status == STATUS_CONVERGED)
    call check_true(passed, "Reference plastic increment", history%state_np1%strain_pf > 0D0)
    call reference%tangent(strain, history, 0.1D0, expected)
    call solver%tangent(strain, history, 0.1D0, tangent)
    call reference%tangent_numerical(strain, history, 0.1D0, numerical)
    call check_approx(passed, "Joint tangent uses the unshifted parent", tangent, expected, tol=TOL_ROUTE)
    call check_approx(passed, "Mixed plastic tangent", tangent, numerical, tol=TOL_TANGENT)
end subroutine test_perzyna_tangent_joint_route_mixed

subroutine test_perzyna_numerical_tangent_elastic_mixed(passed)
    ! Validate the numerical reference against Hooke's known tangent with all six components.
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_yield_vonmises, only : VonMises
    use muscle_hard_bilinear, only : Bilinear_hardening
    use muscle_vp_perzyna, only : Perzyna_overstress
    use muscle_elasticity_linear, only : Elasticity_linear
    use muscle_plastic_history, only : Plastic_material_history
    use muscle_solver_closest_point_perzyna, only : Closest_point
    use muscle_test_utils, only : check_approx
    implicit none
    logical, intent(out) :: passed
    ! Central differences of linear elasticity. Measured: 3.45e-13 relative.
    real(real64), parameter :: TOL = 1D-8
    type(Closest_point) :: solver
    type(Perzyna_overstress) :: vp
    type(Elasticity_linear) :: elas
    type(VonMises) :: vm
    type(Plastic_material_history) :: history
    type(ten_3D2Osym) :: strain
    type(ten_3D4O2sym) :: numerical, expected

    passed = .true.
    call elas%set_parameters(young=1000D0, poisson=0.3D0)
    vp%hard_law = Bilinear_hardening(y0=50D0, k=100D0)
    call solver%init(elas, vp, vm)
    call strain%init(xx=9D-4, yy=-2.5D-4, zz=-1.5D-4, xy=2.5D-4, yz=-1.5D-4, xz=1D-4)
    call history%init()
    expected = elas%dstress_dstrain(strain)
    call solver%tangent_numerical(strain, history, 0.1D0, numerical)
    call check_approx(passed, "Numerical reference agrees with Hooke", numerical, expected, tol=TOL)
end subroutine test_perzyna_numerical_tangent_elastic_mixed
