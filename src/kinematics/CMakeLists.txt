# Historical OO API (host only, removed in Phase 7, D-017.1): the polymorphic bases (their bindings
# call the kernel procedures) and the umbrella module
target_sources(muscle PRIVATE
                 muscle_kinematics_base.f90
                 muscle_kin_finite_base.f90
                 muscle_kinematics.f90
                 )

# Kernel state and common functions, the kinematics written with the model contract macros
# (include/muscle_models_common.fypp; OO adapters from muscle_kin_adapters.fypp) and the family
# module with the generics kin_* (D-017.14)
muscle_add_fypp_sources(muscle
                 muscle_kin_finite.fypp
                 muscle_kin_small_strain.fypp
                 muscle_kin_total_lagrangian.fypp
                 muscle_kin_updated_lagrangian.fypp
                 muscle_kin_material_logarithmic.fypp
                 muscle_kin_spatial_logarithmic.fypp
                 muscle_kin_corotational.fypp
                 muscle_kin_laws.fypp
                 DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/muscle_kin_adapters.fypp
                 )
