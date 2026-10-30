! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_iden_4o3ts
    !! Compatibility module: historical name `iden_4O3TS` (= `iden_4O3TS_dp`) and the operations of
    !! [[muscle_tensor_iden]], where the scaled identities are defined for every kind (D-003, D-014).
    !!
    !! Kept while the operator modules and user code refer to it; new code should
    !! `use muscle_tensors`.
    use muscle_tensor_iden, only: iden_4O3TS => iden_4O3TS_dp
    use muscle_tensor_iden, only: operator(*), operator(/)
#ifdef ENABLE_UDTIO
    use muscle_tensor_io, only: write(formatted)
#endif
    implicit none
    public
end module muscle_tensor_iden_4o3ts
