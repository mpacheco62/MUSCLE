set(HARDENING_LAWS   plasticity/hardening_laws/basis_hardening_law.f90
                     plasticity/hardening_laws/swift_hardening.f90
                     plasticity/hardening_laws/bilinear_hardening.f90
                     CACHE INTERNAL "")


set(RETURN_MAPPING   plasticity/return_mapping/cutting_plane.f90
                     plasticity/return_mapping/closest_point.f90
                     plasticity/return_mapping/closest_point_2.f90
                     CACHE INTERNAL "")


set(VISCOHARDENING_LAWS   plasticity/viscohardening_laws/basis_visco_hardening_law.f90
                          plasticity/viscohardening_laws/JohnsonCook_full_visco_hardening.f90
                          plasticity/viscohardening_laws/JohnsonCook_visco_hardening.f90
                          CACHE INTERNAL "")