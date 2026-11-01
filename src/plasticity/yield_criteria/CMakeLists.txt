# Polymorphic base of the historical OO API (host only, removed in Phase 7, D-017.1)
target_sources(muscle PRIVATE
                 muscle_yield_base.f90
                 )

# Criteria written with the model contract macros (include/muscle_models_common.fypp) and the
# family module with the generics yield_* (D-017); muscle_anisotropic_yield: linear transformations
# and principal-value terms shared by the anisotropic criteria
muscle_add_fypp_sources(muscle
                 muscle_anisotropic_yield.fypp
                 muscle_yield_vonmises.fypp
                 muscle_yield_hill48.fypp
                 muscle_yield_druckerprager.fypp
                 muscle_yield_cpb06.fypp
                 muscle_yield_cpb06ex2.fypp
                 muscle_yield_yld2004.fypp
                 muscle_yield_criteria.fypp
                 )
