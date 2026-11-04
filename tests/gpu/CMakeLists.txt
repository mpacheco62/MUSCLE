# Target tests of the model families (testing.md, GPU): the generics of every model of a family in an
# OpenMP target region against the CPU, one test per family and kind (add_muscle_fypp_target_tests,
# tests/CMakeLists.txt; the cases are in target_cases.py). The kinematics: kinematics-kin_kinds.
add_muscle_fypp_target_tests(family_target test_family_target.fypp family
                             DEPENDS ${BASE_TESTS_DIR}/materials/mat_golden_common.fypp
                                     ${BASE_TESTS_DIR}/materials/golden_params.py)
