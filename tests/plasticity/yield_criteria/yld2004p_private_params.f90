! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program yld2004p_private_params
    !! Must NOT compile (ctest `plasticity-yield_criteria-yld2004p_private_params` expects the compiler
    !! error): the derived parameters Lp and Lpp of the kernel type of Yld2004 in polynomial form (D-055) are
    !! private (D-017.10b), so the type can only be built with `yield_yld2004p_init`.
    use muscle_kinds, only : dp
    use muscle_yield_laws, only : yield_yld2004p_dp
    implicit none
    type(yield_yld2004p_dp) :: c

    c%Lp = 0.0_dp
end program yld2004p_private_params
