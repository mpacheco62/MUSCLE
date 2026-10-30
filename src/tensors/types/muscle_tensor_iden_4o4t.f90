! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_iden_4o4t
    !! Compatibility module: the identity `iden_4O4T` (one type for every kind, D-014) and the
    !! operations of [[muscle_tensor_iden]], where the identities are defined.
    !!
    !! Kept while the operator modules and user code refer to it; new code should
    !! `use muscle_tensors`.
    use muscle_tensor_iden, only: iden_4O4T
    use muscle_tensor_iden, only: operator(*), operator(/)
#ifdef ENABLE_UDTIO
    use muscle_tensor_io, only: write(formatted)
#endif
    implicit none
    public
end module muscle_tensor_iden_4o4t
