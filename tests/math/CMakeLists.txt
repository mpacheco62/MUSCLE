add_muscle_test(muscle_math_operations test_muscle_math_operations.f90)
add_muscle_fypp_test(muscle_math_derivatives test_muscle_math_derivatives.fypp)
add_muscle_test(muscle_math_spectral_derivs test_muscle_math_spectral_derivs.f90)
add_muscle_fypp_test(muscle_math_inverses test_muscle_math_inverses.fypp)
# Small fixed-size LU solvers for kernels (muscle_math_lu_solve), sp and dp
add_muscle_fypp_test(math_lu_solve test_math_lu_solve.fypp)
# Finite-difference macros of include/muscle_fd.fypp (in MUSCLE_FYPP_DEPENDS) against analytical
# derivatives
add_muscle_fypp_test(math_fd_macros test_math_fd_macros.fypp)

# Golden (characterization) tests of the public routines of the math modules, double precision,
# bit for bit against references captured from the original code (tests/reference/math/,
# regenerate.sh there). The finite differences (derivatives, operations) are the macros of
# include/muscle_fd.fypp, which call the test functions directly: the references were captured
# through a procedure argument, so the test functions must not be inlined into the macro stencils
# (with inlining, -O3 vectorization changes the last bits, measured 4e-11 relative, D-016.9).
# -fno-inline is accepted by gfortran and ifx (IntelLLVM; ifx also takes -inline-level=0).
foreach(_module IN ITEMS derivatives operations spectral_derivs linsolvers)
    if (_module MATCHES "^(derivatives|operations)$")
        add_muscle_fypp_test(math_golden_${_module} test_math_golden_${_module}.fypp)
        if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Intel")
            target_compile_options(test_math_golden_${_module} PRIVATE -fno-inline)
        endif()
    else()
        add_muscle_test(math_golden_${_module} test_math_golden_${_module}.f90)
    endif()
    set_tests_properties(math-math_golden_${_module} PROPERTIES
                         ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")
endforeach()

# Spectral routines (eigenvals, spectral_decomposition, eigenvalue derivatives, spectral
# gradient/Hessian) in sp and dp at stress scales 1e-3..2.1e5 with distinct, nearly repeated and
# repeated eigenvalues, against independent dp references (D-016.6)
add_muscle_fypp_test(math_spectral_kinds test_math_spectral_kinds.fypp)
# eigenvals with exactly repeated (largest or smallest), triple and nearly repeated roots, exact
# references
add_muscle_fypp_test(math_eigenvals_repeated test_math_eigenvals_repeated.fypp)
# Abort on a division by zero, overflow or invalid operation in the library: gfortran enables the
# traps from the options of the main program
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_math_spectral_kinds PRIVATE -ffpe-trap=zero,overflow,invalid)
    target_compile_options(test_math_eigenvals_repeated PRIVATE -ffpe-trap=zero,overflow,invalid)
endif()
# Math kernel routines (eigenvals, spectral_decomposition, dEigenvalues_dTensor, FD_GRAD, LU6SOLVE)
# called from a pure point procedure inside an OpenMP target region, sp and dp, bit for bit
# against the CPU call
add_muscle_fypp_test(math_kernels_target test_math_kernels_target.fypp)
muscle_target_test_options(test_math_kernels_target)
