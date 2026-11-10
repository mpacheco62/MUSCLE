add_subdirectory(hardening_laws)
add_subdirectory(viscoplastic_laws)
add_subdirectory(yield_criteria)

# Kernel state of plasticity (muscle_plastic_state) in every kind: defaults and a pure @:DEVICE() point
# procedure (the FE history array layout is that of each material, tests/materials)
add_muscle_fypp_test(plastic_state_kinds test_muscle_plastic_state_kinds.fypp)
