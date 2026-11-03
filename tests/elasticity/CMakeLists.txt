# Physical checks of the linear model through the kernel API (Hooke law, 3D and 2D, OpenMP)
add_muscle_test(linear test_muscle_elasticity_linear.f90)
add_muscle_test(linear_openmp test_muscle_elasticity_linear_openmp.f90)

# Golden (characterization) test of the elasticity family (every function of the linear model
# through the kernel API, 3D and 2D), double precision, bit for bit against
# tests/reference/elasticity/ (regenerate.sh there), captured from the original OO code.
add_muscle_test(elas_golden test_elas_golden.f90)
set_tests_properties(elasticity-elas_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the elastic family (D-017.12) in every kind and dimension: the linear model through
# the generics of muscle_elasticity_laws against Hooke law (dp reference), 2D vs 3D embedding of
# every function, the generics against the dp specific functions (bit for bit), and pure @:DEVICE()
# point procedures
add_muscle_fypp_test(elas_kinds test_elas_kinds.fypp)

# Compile-fail check (D-017.10b): every component of the kernel type elas_linear_<k> is derived
# (private), so it can only be built with elas_linear_init. Passes when the compiler rejects
# linear_private_params.f90 with its "PRIVATE component" error; any other error (a missing module, a
# fatal error) fails it. Only with gfortran (the message of other compilers is not known here; the
# file name itself contains "private", so a loose pattern would always pass).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_test(NAME elasticity-linear_private_params
             COMMAND ${CMAKE_Fortran_COMPILER} -fsyntax-only
                     -I$<TARGET_PROPERTY:muscle,Fortran_MODULE_DIRECTORY>
                     ${CMAKE_CURRENT_SOURCE_DIR}/linear_private_params.f90)
    set_tests_properties(elasticity-linear_private_params PROPERTIES
                         PASS_REGULAR_EXPRESSION "PRIVATE component"
                         FAIL_REGULAR_EXPRESSION "Cannot open module|Fatal Error")
endif()
