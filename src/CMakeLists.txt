# set(OPERATORS tensors_types/mixed_operators/operator_I3D2O_3D2Osym.f90
#               CACHE INTERNAL "")

set(TENSORS_TYPES tensors_types/basic_types/ten_3D2O.f90
                  tensors_types/basic_types/ten_3D2Osym.f90
                  tensors_types/basic_types/ten_3D4O3sym.f90
                  tensors_types/basic_types/ten_3D4O2sym.f90
                  tensors_types/basic_types/iden_3D2OMod.f90
                  tensors_types/basic_types/iden_3D2O.f90
                  tensors_types/basic_types/iden_3D4O3TMod.f90
                  tensors_types/basic_types/iden_3D4O3T.f90
                  tensors_types/basic_types/iden_3D4O4TMod.f90
                  tensors_types/basic_types/iden_3D4O4T.f90
                  tensors_types/mixed_operators/operator_I3D2O_3D2Osym.f90
                  tensors_types/mixed_operators/operator_I3D2O_3D2O.f90
                  tensors_types/mixed_operators/operator_I3D2OMod_3D2Osym.f90
                  tensors_types/mixed_operators/operator_I3D2OMod_3D2O.f90
                  tensors_types/mixed_operators/operator_I3D4O3T_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_I3D4O3T_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_I3D4O3TMod_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_I3D4O3TMod_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_I3D4O4T_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_I3D4O4T_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_I3D4O4TMod_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_I3D4O4TMod_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_3D2Osym_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_3D2Osym_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_3D4O2sym_3D4O3sym.f90
                  tensors_types/tensors_types.f90
                  CACHE INTERNAL "")

set(MY_MODS   basic_operations.f90
              basis_yield_criteria.f90
              basis_hardening_law.f90
              swift_hardening.f90
              vonMises.f90
              CPB06.f90
              cutting_plane.f90
              closest_point.f90
              basis_elasticity.f90
              inverses_mat.f90
              elasticity_linear.f90
              CACHE INTERNAL "")

      #set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(umatlib ${TENSORS_TYPES} ${MY_MODS})
set_target_properties(umatlib PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/umatlib)
install(TARGETS umatlib LIBRARY)
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/umatlib/ TYPE INCLUDE)


add_executable(my_exe main.f90)
target_link_libraries(my_exe PRIVATE umatlib)


