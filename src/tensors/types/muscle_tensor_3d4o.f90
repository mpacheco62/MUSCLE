! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_3d4o
    !! Compatibility module: historical name `ten_3D4O` (= `ten_3D4O_dp`) and the
    !! operations of [[muscle_tensor_4o]], where the types are defined for every kind (D-003).
    !!
    !! Kept while the operator modules and user code refer to the unsuffixed name; new code should
    !! `use muscle_tensors`.
    use muscle_tensor_4o, only: ten_3D4O => ten_3D4O_dp
    use muscle_tensor_4o, only: operator(.approx.), operator(+), operator(-), operator(*), &
                                operator(/), operator(.inv.), assignment(=), inverse
#ifdef ENABLE_UDTIO
    use muscle_tensor_io, only: write(formatted)
#endif
    implicit none
    public
end module muscle_tensor_3d4o
