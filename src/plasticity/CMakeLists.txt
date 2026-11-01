target_sources(muscle PRIVATE
                 muscle_plastic_history.f90
                 )

# Kernel state of a plastic material point (support module, every kind)
muscle_add_fypp_sources(muscle
                 muscle_plastic_state.fypp
                 )

add_subdirectory(hardening_laws)
add_subdirectory(viscoplastic_laws)
add_subdirectory(yield_criteria)