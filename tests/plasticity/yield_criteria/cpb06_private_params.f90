! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program cpb06_private_params
    !! Must NOT compile (ctest `plasticity-yield_criteria-cpb06_private_params` expects the compiler
    !! error): the derived parameters B and L of the CPB06 kernel type are private (D-017.10b), so the
    !! type can only be built with `yield_cpb06_init`, not with its structure constructor.
    use muscle_kinds, only : dp
    use muscle_yield_laws, only : yield_cpb06_dp
    implicit none
    type(yield_cpb06_dp) :: c

    c = yield_cpb06_dp(C1=0.0_dp, C2=0.0_dp, k=0.0_dp, a=2.0_dp, B=1.0_dp, L=0.0_dp)
end program cpb06_private_params
