! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program linear_private_params
    !! Must NOT compile (ctest `elasticity-linear_private_params` expects the compiler error): the
    !! Lame parameters and cached tangents of the linear elastic kernel type are private (D-017.10b),
    !! so the type can only be built with `elas_linear_init`.
    use muscle_kinds, only : dp
    use muscle_elas_laws, only : elas_linear_dp
    implicit none
    type(elas_linear_dp) :: p

    p%mu = 1.0_dp
end program linear_private_params
