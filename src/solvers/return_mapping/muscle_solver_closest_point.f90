! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_solver_closest_point
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_base, only : Base_hardening_laws
    use muscle_yield_base, only : Base_yield_critera
    use muscle_elasticity_base, only : Base_elasticity
    use muscle_plastic_history, only : Plastic_material_history, Plastic_state

    private

    integer, parameter, public :: STATUS_UNDEFINED=-1
    integer, parameter, public :: STATUS_ELASTIC_CASE=0
    integer, parameter, public :: STATUS_CONVERGED=1
    integer, parameter, public :: STATUS_NONCONVERGED=2
    integer, parameter, public :: STATUS_ITER_CONVERGED    = 10
    integer, parameter, public :: STATUS_ITER_NONCONVERGED = 11


    type, public :: Closest_point
        class(Base_elasticity), allocatable     :: elasticity
        class(Base_hardening_laws), allocatable :: hardening
        class(Base_yield_critera), allocatable  :: yield
        integer, private                        :: iter_nw = 200
        
        contains
        procedure, public :: init => closest_point_init
        procedure, public :: solve => closest_point_solve
        procedure, public :: iter => closest_point_iter
        procedure, public :: tangent => closest_point_tangent
        procedure, public :: tangent_numerical => closest_point_tangent_numerical
    end type Closest_point

    public :: closest_point2

    contains

        pure subroutine closest_point_init(self, elasticity, hardening, yield, iter_nw)
            implicit none
            class(Closest_point), intent(inout) :: self
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

        include 'closest_point_methods.inc'

end module