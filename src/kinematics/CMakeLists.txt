# Kernel state and common functions, the kinematics and the family module
# with the generics kin_*
muscle_add_fypp_sources(muscle
                 muscle_kin_finite.fypp
                 muscle_kin_small_strain.fypp
                 muscle_kin_total_lagrangian.fypp
                 muscle_kin_updated_lagrangian.fypp
                 muscle_kin_material_logarithmic.fypp
                 muscle_kin_spatial_logarithmic.fypp
                 muscle_kin_corotational.fypp
                 muscle_kin_laws.fypp
                 )
