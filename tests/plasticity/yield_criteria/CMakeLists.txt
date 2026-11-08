# Numeric derivatives of every criterion (yield_test_numeric: FD_GRAD/FD_HESS of yield_stress_eq,
# the steps of the retired OO base), for the physical checks of the criteria below
muscle_fypp_preprocess(YIELD_TEST_NUMERIC_SRC yield_test_numeric.fypp)
add_library(muscle_yield_test_numeric STATIC ${YIELD_TEST_NUMERIC_SRC})
target_link_libraries(muscle_yield_test_numeric PUBLIC muscle PRIVATE compile_options)
set_target_properties(muscle_yield_test_numeric PROPERTIES
                      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/tests)
target_include_directories(muscle_yield_test_numeric PUBLIC ${CMAKE_BINARY_DIR}/modules/tests)

# Physical checks of the criteria through the kernel API (closed forms, r-values and shear of the
# anisotropic criteria, normalization, analytical derivatives against the numeric ones)
foreach(_muscle_yield CPB06 HILL48 VONMISES DRUCKERPRAGER YLD2004 CPB06EX2)
    string(TOLOWER ${_muscle_yield} _muscle_yield_file)
    add_muscle_test(${_muscle_yield} test_muscle_yield_${_muscle_yield_file}.f90)
    target_link_libraries(test_${_muscle_yield} PRIVATE muscle_yield_test_numeric)
endforeach()

# Golden (characterization) test of every yield criterion (stress_eq, dstress_eq, ddstress_eq,
# stress_eq_and_derivatives, joint_evaluation through the kernel API, constructors), double precision,
# within the tolerance of its category (D-021) against tests/reference/plasticity/yield_criteria/ (regenerate.sh there), captured from
# the original OO code. Includes the zero tensor and states with repeated principal values.
add_muscle_fypp_test(yield_golden test_yield_golden.fypp)
set_tests_properties(plasticity-yield_criteria-yield_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the yield criteria (D-017) in every kind: free functions through the generics of
# muscle_yield_criteria against closed forms and central differences (dp) and against dp (sp), zero
# derivatives at a null equivalent stress (D-017.9), default joint evaluation, finite-difference
# finite-difference derivatives (FD_GRAD/FD_HESS of the stress), and a pure @:DEVICE() point procedure
add_muscle_fypp_test(yield_kinds test_yield_kinds.fypp)

# Compile-fail checks (D-017.10b): the kernel types with derived parameters (CPB06: B, L;
# Drucker-Prager: f_factor, I1_factor; CPB06ex2: L, Lp, B; Yld2004: Lp, Lpp) have them private, so
# they cannot be set outside the model module (only by yield_<m>_init). Each test passes when the
# compiler rejects <m>_private_params.f90 with its "PRIVATE component" error; any other error (a
# missing module, a fatal error) fails it. Only with gfortran (the message of other compilers is
# not known here; the file name itself contains "private", so a loose pattern would always pass).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    foreach(_muscle_model cpb06 druckerprager cpb06ex2 yld2004)
        add_test(NAME plasticity-yield_criteria-${_muscle_model}_private_params
                 COMMAND ${CMAKE_Fortran_COMPILER} -fsyntax-only
                         -I$<TARGET_PROPERTY:muscle,Fortran_MODULE_DIRECTORY>
                         ${CMAKE_CURRENT_SOURCE_DIR}/${_muscle_model}_private_params.f90)
        set_tests_properties(plasticity-yield_criteria-${_muscle_model}_private_params PROPERTIES
                             PASS_REGULAR_EXPRESSION "PRIVATE component"
                             FAIL_REGULAR_EXPRESSION "Cannot open module|Fatal Error")
    endforeach()
endif()
