# Polymorphic bases of the historical OO API (host only, removed in Phase 7, D-017.1)
target_sources(muscle PRIVATE
                 muscle_vp_base.f90
                 muscle_vp_overstress_base.f90
                 )

# Laws written with the model contract macros (include/muscle_models_common.fypp) and the family
# modules with the generics vp_* (flow stress laws) and os_* (overstress laws) (D-017)
muscle_add_fypp_sources(muscle
                 muscle_vp_jc.fypp
                 muscle_vp_jcfull.fypp
                 muscle_vp_mrk.fypp
                 muscle_vp_nnl.fypp
                 muscle_vp_rk.fypp
                 muscle_vp_va.fypp
                 muscle_vp_laws.fypp
                 muscle_vp_chaboche.fypp
                 muscle_vp_perzyna.fypp
                 muscle_vp_overstress_laws.fypp
                 )
