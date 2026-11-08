# Kernel types (state, coefficients), the frozen elasticity (an elastic model), the models and the
# family module with the generics visco_*
muscle_add_fypp_sources(muscle
                 muscle_visco_types.fypp
                 muscle_visco_frozen.fypp
                 muscle_visco_maxwell.fypp
                 muscle_visco_kelvin.fypp
                 muscle_visco_laws.fypp
                 )
