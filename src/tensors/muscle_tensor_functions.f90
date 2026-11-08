! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_functions
    !! Free (non type-bound) tensor functions of every family, in one `use`.
    !!
    !! These generics are the kernel-friendly API: `type` arguments, `pure`, callable from GPU
    !! regions. They are not exported by [[muscle_tensors]] so that they do not clash with user
    !! variables named `norm`, `det`, ...; import them explicitly, renaming if needed:
    !!
    !! ```fortran
    !! use muscle_tensors
    !! use muscle_tensor_functions, only: norm, inverse
    !! use muscle_tensor_functions, only: tdet => det
    !! ```
    !!
    !! - `norm(a)`: every family.
    !! - `is_approx(a, b [, tol])`, `call inverse(a, a_inv, ok [, cond])`: every family.
    !! - `det(a)`: `2Osym`, `2O`. `square(a)`: `2Osym`. `transpose(a)`: `2O` (extends the intrinsic).
    !! - `convert_3sym(a)`: `4O2sym` → `4O3sym`. `convert_2sym(a)`: `4O` → `4O2sym`.
    use muscle_tensor_2osym, only: norm, det, square, is_approx, inverse
    use muscle_tensor_2o, only: norm, det, transpose, is_approx, inverse
    use muscle_tensor_4o3sym, only: norm, is_approx, inverse
    use muscle_tensor_4o2sym, only: norm, convert_3sym, is_approx, inverse
    use muscle_tensor_4o, only: norm, convert_2sym, is_approx, inverse
    implicit none
    private

    public :: norm, det, square, transpose, is_approx, inverse, convert_3sym, convert_2sym
end module muscle_tensor_functions
