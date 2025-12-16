# set(OPERATORS tensors_types/mixed_operators/operator_I2O_3D2Osym.f90
#               CACHE INTERNAL "")

set(HARDENING_LAWS   hardening_laws/basis_hardening_law.f90
                     hardening_laws/swift_hardening.f90
                     hardening_laws/bilinear_hardening.f90
                     hardening_laws/basis_viscoplastic_law.f90
                     hardening_laws/JC_model.f90
                     hardening_laws/MRK_model.f90
                     hardening_laws/NNL_model.f90
                     hardening_laws/RK_model.f90
                     hardening_laws/VA_model.f90
                     hardening_laws/voce_m_hardening.f90
                     CACHE INTERNAL "")

# add_library(umatlib ${HARDENING_LAWS})

