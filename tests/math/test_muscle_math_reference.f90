! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_math_reference
    !! Minimal reference test of the math modules (plan R): the three results whose mutations no unit
    !! test can detect, compared within the algebra tolerance (D-021) with
    !! `tests/reference/math/math_reference_dp.txt` (see its header for why each one is kept):
    !! `eigenvals` with the middle eigenvalue near the mean of the others, and the eigenvectors that
    !! `spectral_decomposition` selects in two exact ties. Everything else of these modules is checked
    !! against independent truths (math_spectral_cases, math_spectral_kinds, math_lu_solve, ...).
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym
    use muscle_math_eigen, only : eigenvals, spectral_decomposition
    use muscle_test_reference, only : reference_init, reference_check, reference_finish
    implicit none

    type(ten_3D2Osym) :: t

    call reference_init('math/math_reference_dp.txt')

    t%vals = [5.0_dp, 3.0_dp + 1.0e-11_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call reference_check('input eigenvals middle near mean | case 1', t%vals)
    call reference_check('eigenvals middle near mean | case 1', eigenvals(t))

    t%vals = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    call reference_check('input middle = mean | case 1', t%vals)
    call check_vectors('spectral_decomposition V middle = mean | case 1', t)

    t%vals = [2.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    call reference_check('input double exact | case 1', t%vals)
    call check_vectors('spectral_decomposition V double exact | case 1', t)

    call reference_finish()

contains

    subroutine check_vectors(key, t)
        !! Checks the eigenvectors of `spectral_decomposition`.
        character(len=*), intent(in) :: key
            !! Record key.
        type(ten_3D2Osym), intent(in) :: t
            !! Tensor.
        real(dp) :: lam(3), v(3,3)
        call spectral_decomposition(t, lam, v)
        call reference_check(key, reshape(v, [9]))
    end subroutine check_vectors

end program test_math_reference
