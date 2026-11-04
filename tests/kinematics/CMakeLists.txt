# Physical checks of the six kinematics through the kernel API (strain measures, Cauchy stress,
# push-forward of stress and tangent, conversion roundtrips, inverse quantities of F)
add_muscle_test(kinematics test_muscle_kinematics.f90)

# Golden (characterization) test of the kinematics family (the state and the family and
# muscle_kin_finite functions of the six kinematics through the kernel API, on deterministic
# deformation gradients and incremental sequences), double precision, bit for bit against
# tests/reference/kinematics/ (regenerate.sh there), captured from the original OO code.
add_muscle_fypp_test(kin_golden test_kin_golden.fypp)
set_tests_properties(kinematics-kin_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the kinematics family (D-017.14) in every kind: sp against dp, the status of invalid
# deformation gradients (D-017.14a), objectivity, and pure @:DEVICE() point procedures in an OpenMP
# target region (muscle_target_test_options, tests/CMakeLists.txt)
add_muscle_fypp_test(kin_kinds test_kin_kinds.fypp)
muscle_target_test_options(test_kin_kinds)
