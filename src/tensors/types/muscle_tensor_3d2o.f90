! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_3d2o
    !! Compatibility module: historical name `ten_3D2O` (= `ten_3D2O_dp`) and the
    !! operations of [[muscle_tensor_2o]], where the types are defined for every kind (D-003).
    !!
    !! Kept while the operator modules and user code refer to the unsuffixed name; new code should
    !! `use muscle_tensors`.
    use muscle_tensor_2o, only: ten_3D2O => ten_3D2O_dp
    use muscle_tensor_2o, only: operator(.approx.), operator(+), operator(-), operator(*), &
                                operator(/), operator(.dev.), operator(.ddot.), operator(.inv.), &
                                assignment(=), inverse, transpose
#ifdef ENABLE_UDTIO
    use muscle_tensor_io, only: write(formatted)
#endif
    implicit none
    public
end module muscle_tensor_3d2o
