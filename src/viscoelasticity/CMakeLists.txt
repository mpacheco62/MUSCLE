# Kernel types (state, coefficients), the frozen elasticity (an elastic model), the models written
# with the model contract macros (include/muscle_models_common.fypp) and the family module with the
# generics visco_* (D-017.13)
muscle_add_fypp_sources(muscle
                 muscle_visco_types.fypp
                 muscle_visco_frozen.fypp
                 muscle_visco_maxwell.fypp
                 muscle_visco_kelvin.fypp
                 muscle_visco_laws.fypp
                 )
