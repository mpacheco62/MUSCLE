! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program druckerprager_private_params
    !! Must NOT compile (ctest `plasticity-yield_criteria-druckerprager_private_params` expects the
    !! compiler error): the derived parameters f_factor and I1_factor of the Drucker-Prager kernel type
    !! are private (D-017.10b), so the type can only be built with `yield_druckerprager_init`.
    use muscle_kinds, only : dp
    use muscle_yield_laws, only : yield_druckerprager_dp
    implicit none
    type(yield_druckerprager_dp) :: c

    c%f_factor = 1.0_dp
end program druckerprager_private_params
