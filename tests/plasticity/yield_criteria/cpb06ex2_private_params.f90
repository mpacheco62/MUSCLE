! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program cpb06ex2_private_params
    !! Must NOT compile (ctest `plasticity-yield_criteria-cpb06ex2_private_params` expects the compiler
    !! error): the derived parameters L, Lp and B of the CPB06ex2 kernel type are private (D-017.10b),
    !! so the type can only be built with `yield_cpb06ex2_init`.
    use muscle_kinds, only : dp
    use muscle_yield_laws, only : yield_cpb06ex2_dp
    implicit none
    type(yield_cpb06ex2_dp) :: c

    c%B = 1.0_dp
end program cpb06ex2_private_params
