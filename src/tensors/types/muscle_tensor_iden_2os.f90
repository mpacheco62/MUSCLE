! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_iden_2os
    !! Compatibility module: historical name `iden_2OS` (= `iden_2OS_dp`) and the operations of
    !! [[muscle_tensor_iden]], where the scaled identities are defined for every kind (D-003, D-014).
    !!
    !! Kept while the operator modules and user code refer to it; new code should
    !! `use muscle_tensors`.
    use muscle_tensor_iden, only: iden_2OS => iden_2OS_dp
    use muscle_tensor_iden, only: operator(+), operator(-), operator(*), operator(/), &
                                  operator(.approx.)
#ifdef ENABLE_UDTIO
    use muscle_tensor_io, only: write(formatted)
#endif
    implicit none
    public
end module muscle_tensor_iden_2os
