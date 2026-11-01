add_muscle_test(kinematics test_muscle_kinematics.f90)

# Golden (characterization) test of the kinematics family (every public and type-bound procedure of
# Base_kinematics, Base_F_kinematics and the six concrete kinematics through the polymorphic bases, on
# deterministic deformation gradients and incremental sequences), double precision, bit for bit against
# tests/reference/kinematics/ (regenerate.sh there), captured from the original OO code.
add_muscle_test(kin_golden test_kin_golden.f90)
set_tests_properties(kinematics-kin_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the kinematics family (D-017.14) in every kind: the dp kernel against the historical
# OO types (bit for bit), sp against dp, the status of invalid deformation gradients (D-017.14a),
# objectivity, and pure @:DEVICE() point procedures in an OpenMP target region (host fallback only,
# as math-math_kernels_target: -foffload=disable, -fno-inline for one copy of the point procedures)
add_muscle_fypp_test(kin_kinds test_kin_kinds.fypp)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_kin_kinds PRIVATE -foffload=disable -fno-inline)
    target_link_options(test_kin_kinds PRIVATE -foffload=disable)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(test_kin_kinds PRIVATE -fno-inline)
endif()
