target_sources(muscle PRIVATE
                 muscle_hard_base.f90
                 )

# Laws written with the hardening_model contract macro (include/muscle_models_common.fypp) and
# the family module with the generics stress, dstress_dep, ddstress_ddep (D-017)
muscle_add_fypp_sources(muscle
                 muscle_hard_bilinear.fypp
                 muscle_hard_ludwik.fypp
                 muscle_hard_swift.fypp
                 muscle_hard_vocemod.fypp
                 muscle_hard_laws.fypp
                 )
