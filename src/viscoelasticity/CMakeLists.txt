# Historical OO API (host only, removed in Phase 7, D-017.1): the transactional history container
# and the polymorphic base
target_sources(muscle PRIVATE
                 muscle_visco_state.f90
                 muscle_visco_base.f90
                 )

# Kernel types (state, coefficients), the frozen elasticity (an elastic model), the models written
# with the model contract macros (include/muscle_models_common.fypp; OO adapters from
# muscle_visco_adapters.fypp) and the family module with the generics visco_* (D-017.13)
muscle_add_fypp_sources(muscle
                 muscle_visco_types.fypp
                 muscle_visco_frozen.fypp
                 muscle_visco_maxwell.fypp
                 muscle_visco_kelvin.fypp
                 muscle_visco_laws.fypp
                 DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/muscle_visco_adapters.fypp
                 )
