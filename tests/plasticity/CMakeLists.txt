add_subdirectory(hardening_laws)
add_subdirectory(viscoplastic_laws)
add_subdirectory(yield_criteria)

add_muscle_test(plastic_history test_muscle_plastic_history.f90)

# Golden (characterization) test of muscle_plastic_history (defaults, init with every combination of
# optional arguments, commit, rollback, unpack_from_fea/pack_to_fea and multi-step transactional
# sequences), double precision, bit for bit against tests/reference/plasticity/ (regenerate.sh
# there), captured from the original OO code.
add_muscle_test(plastic_history_golden test_plastic_history_golden.f90)
set_tests_properties(plasticity-plastic_history_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel state of plasticity (muscle_plastic_state) in every kind: FE array layout, dp against the
# historical Plastic_material_history (bit for bit), sp against the sp-rounded dp state and a pure
# @:DEVICE() point procedure
add_muscle_fypp_test(plastic_state_kinds test_plastic_state_kinds.fypp)
