# Polymorphic base of the historical OO API (host only, removed in Phase 7, D-017.1)
target_sources(muscle PRIVATE
                 muscle_elasticity_base.f90
                 )

# Models written with the model contract macros (include/muscle_models_common.fypp) and the
# family module with the generics elas_* (D-017)
muscle_add_fypp_sources(muscle
                 muscle_elasticity_linear.fypp
                 muscle_elasticity_laws.fypp
                 )
