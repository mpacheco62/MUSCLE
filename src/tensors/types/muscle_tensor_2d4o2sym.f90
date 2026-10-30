! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_2d4o2sym
    !! Compatibility module: historical name `ten_2D4O2sym` (= `ten_2D4O2sym_dp`) and the
    !! operations of [[muscle_tensor_4o2sym]], where the types are defined for every kind (D-003).
    !!
    !! Kept while the operator modules and user code refer to the unsuffixed name; new code should
    !! `use muscle_tensors`.
    use muscle_tensor_4o2sym, only: ten_2D4O2sym => ten_2D4O2sym_dp
    use muscle_tensor_4o2sym, only: operator(.approx.), operator(+), operator(-), operator(*), &
                                    operator(/), operator(.inv.), assignment(=), inverse
#ifdef ENABLE_UDTIO
    use muscle_tensor_io, only: write(formatted)
#endif
    implicit none
    public
end module muscle_tensor_2d4o2sym
