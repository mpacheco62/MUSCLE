target_sources(muscle PRIVATE
                 muscle_plastic_history.f90
                 )

add_subdirectory(hardening_laws)
add_subdirectory(viscoplastic_laws)
add_subdirectory(yield_criteria)