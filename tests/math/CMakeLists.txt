add_muscle_fypp_test(muscle_math_derivatives test_muscle_math_derivatives.fypp)
add_muscle_fypp_test(muscle_math_inverses test_muscle_math_inverses.fypp)
# Small fixed-size LU solvers for kernels (muscle_math_lu_solve), sp and dp
add_muscle_fypp_test(math_lu_solve test_muscle_math_lu_solve.fypp)
# Finite-difference macros of tests/common/muscle_fd.fypp (a dependency of every test, tests/CMakeLists.txt)
# against analytical derivatives
add_muscle_fypp_test(math_fd_macros test_muscle_fd_macros.fypp)

# Minimal reference test (plan R): the results whose mutations no unit test detects (eigenvals with the
# middle eigenvalue near the mean, the eigenvectors chosen at exact ties), within the algebra tolerance
# (D-021) against tests/reference/math/math_reference_dp.txt (regenerate.sh there)
add_muscle_test(math_reference test_muscle_math_reference.f90)
muscle_reference_test(math-math_reference)

# Spectral routines (eigenvals, spectral_decomposition, eigenvalue derivatives, spectral
# gradient/Hessian) in sp and dp at stress scales 1e-3..2.1e5 with distinct, nearly repeated and
# repeated eigenvalues, against independent dp references (D-016.6)
add_muscle_fypp_test(math_spectral_kinds test_muscle_math_spectral_derivs_kinds.fypp)
# eigenvals with exactly repeated (largest or smallest), triple and nearly repeated roots, exact
# references
add_muscle_fypp_test(math_eigenvals_repeated test_muscle_math_spectral_derivs_eigenvals.fypp)
# Special cases of the spectral routines with closed-form references: eigenvals with the middle root
# at the mean (exact), the branches of spectral_decomposition (exactly repeated pairs, ties, diagonal
# tensors), spectral_hessian of tr(S^4)/4 + det(S) with distinct, double and triple roots, and
# spectral_hessian against its documented formula for arbitrary phi_a, phi_ab
add_muscle_fypp_test(math_spectral_cases test_muscle_math_spectral_derivs_cases.fypp)
# Abort on a division by zero, overflow or invalid operation in the library: gfortran enables the
# traps from the options of the main program
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_math_spectral_kinds PRIVATE -ffpe-trap=zero,overflow,invalid)
    target_compile_options(test_math_eigenvals_repeated PRIVATE -ffpe-trap=zero,overflow,invalid)
    target_compile_options(test_math_spectral_cases PRIVATE -ffpe-trap=zero,overflow,invalid)
endif()
# Math kernel routines (eigenvals, spectral_decomposition, dEigenvalues_dTensor, FD_GRAD, LU6SOLVE)
# called from a pure point procedure inside an OpenMP target region, sp and dp, bit for bit
# against the CPU call
add_muscle_fypp_test(math_kernels_target test_muscle_math_target.fypp)
muscle_target_test_options(test_math_kernels_target)
