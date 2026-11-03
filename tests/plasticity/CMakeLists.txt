add_subdirectory(hardening_laws)
add_subdirectory(viscoplastic_laws)
add_subdirectory(yield_criteria)

# Kernel state of plasticity (muscle_plastic_state) in every kind: FE array layout (pack/unpack;
# it replaces the golden of the retired Plastic_material_history, D-018.18f), sp against the
# sp-rounded dp state and a pure @:DEVICE() point procedure
add_muscle_fypp_test(plastic_state_kinds test_plastic_state_kinds.fypp)
