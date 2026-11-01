add_muscle_test(CPB06 test_muscle_yield_cpb06.f90)
add_muscle_test(HILL48 test_muscle_yield_hill48.f90)
add_muscle_test(VONMISES test_muscle_yield_vonmises.f90)
add_muscle_test(DRUCKERPRAGER test_muscle_yield_druckerprager.f90)
add_muscle_test(YLD2004 test_muscle_yield_yld2004.f90)
add_muscle_test(CPB06EX2 test_muscle_yield_cpb06ex2.f90)

# Golden (characterization) test of every yield criterion (stress_eq, dstressEq_dstress,
# ddstressEq_ddstress, stressEq_and_derivatives, joint_evaluation through the polymorphic base,
# the finite-difference procedures of the base, constructors), double precision, bit for bit
# against tests/reference/plasticity/yield_criteria/ (regenerate.sh there), captured from the
# original OO code. Includes the zero tensor and states with repeated principal values.
add_muscle_test(yield_golden test_yield_golden.f90)
set_tests_properties(plasticity-yield_criteria-yield_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the yield criteria (D-017) in every kind: free functions through the generics of
# muscle_yield_criteria against closed forms and central differences (dp) and against dp (sp), zero
# derivatives at a null equivalent stress (D-017.9), default joint evaluation, finite-difference
# derivatives of the contract (MODEL_FD in the stress) against the numeric procedures of the OO
# base, historical OO types against the dp functions, and a pure @:DEVICE() point procedure
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
