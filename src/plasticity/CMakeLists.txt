set(HARDENING_LAWS   plasticity/hardening_laws/basis_hardening_law.f90
                     plasticity/hardening_laws/swift_hardening.f90
                     plasticity/hardening_laws/bilinear_hardening.f90
                     CACHE INTERNAL "")


set(RETURN_MAPPING   plasticity/return_mapping/cutting_plane.f90
                     plasticity/return_mapping/closest_point.f90
                     plasticity/return_mapping/closest_point_2.f90
                     CACHE INTERNAL "")


# viscoplastic models				  
set(VISCOPLASTIC_LAWS   plasticity/viscoplastic_laws/basis_viscoplastic_law.f90
                    	plasticity/viscoplastic_laws/JC_model.f90
                     	plasticity/viscoplastic_laws/JohnsonCook_full_visco_hardening.f90
                     	plasticity/viscoplastic_laws/MRK_model.f90
                     	plasticity/viscoplastic_laws/NNL_model.f90
                     	plasticity/viscoplastic_laws/RK_model.f90
                     	plasticity/viscoplastic_laws/VA_model.f90
                     	plasticity/viscoplastic_laws/voce_m_hardening.f90
                     	CACHE INTERNAL "")
