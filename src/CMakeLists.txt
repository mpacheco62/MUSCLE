# set(OPERATORS tensors_types/mixed_operators/operator_I2O_3D2Osym.f90
#               CACHE INTERNAL "")

set(TENSORS_TYPES tensors_types/basic_types/ten_3D2O.f90
                  tensors_types/basic_types/ten_3D2Osym.f90
                  tensors_types/basic_types/ten_2D2Osym.f90
                  tensors_types/basic_types/ten_3D4O3sym.f90
                  tensors_types/basic_types/ten_2D4O3sym.f90
                  tensors_types/basic_types/ten_3D4O2sym.f90
                  tensors_types/basic_types/iden_2OS.f90
                  tensors_types/basic_types/iden_2O.f90
                  tensors_types/basic_types/iden_4O3TS.f90
                  tensors_types/basic_types/iden_4O3T.f90
                  tensors_types/basic_types/iden_4O4TS.f90
                  tensors_types/basic_types/iden_4O4T.f90
                  tensors_types/mixed_operators/operator_I2O_3D2Osym.f90
                  tensors_types/mixed_operators/operator_I2O_2D2Osym.f90
                  tensors_types/mixed_operators/operator_I2O_3D2O.f90
                  tensors_types/mixed_operators/operator_I2OS_3D2Osym.f90
                  tensors_types/mixed_operators/operator_I2OS_2D2Osym.f90
                  tensors_types/mixed_operators/operator_I2OS_3D2O.f90
                  tensors_types/mixed_operators/operator_I4O3T_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_I4O3T_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_I4O3TS_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_I4O3TS_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_I4O4T_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_I4O4T_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_I4O4TS_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_I4O4TS_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_3D2Osym_3D2O.f90
                  tensors_types/mixed_operators/operator_3D2Osym_3D4O2sym.f90
                  tensors_types/mixed_operators/operator_3D2Osym_3D4O3sym.f90
                  tensors_types/mixed_operators/operator_2D2Osym_2D4O3sym.f90
                  tensors_types/mixed_operators/operator_3D4O2sym_3D4O3sym.f90
                  tensors_types/tensors_types.f90
                  CACHE INTERNAL "")

set(HARDENING_LAWS   hardening_laws/basis_hardening_law.f90
                     hardening_laws/swift_hardening.f90
                     hardening_laws/bilinear_hardening.f90
                     CACHE INTERNAL "")

set(YIELD_CRITERIA   yield_criteria/basis_yield_criteria.f90
                     yield_criteria/vonMises.f90
                     yield_criteria/CPB06.f90
                     CACHE INTERNAL "")

set(ELASTICITY       elasticity/basis_elasticity.f90
                     elasticity/elasticity_linear.f90
                     CACHE INTERNAL "")

set(RETURN_MAPPING   return_mapping/cutting_plane.f90
                     return_mapping/closest_point.f90
                     CACHE INTERNAL "")

set(OPERATORS   operators/basic_operations.f90
                operators/inverses_mat.f90
                operators/derivatives.f90
               CACHE INTERNAL "")

      #set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(umatlib 
                     ${TENSORS_TYPES}
                     ${HARDENING_LAWS}
                     ${YIELD_CRITERIA}
                     ${ELASTICITY}
                     ${RETURN_MAPPING}
                     ${OPERATORS}
                     )
target_compile_options(umatlib PRIVATE 
                                -c -cpp -fopenmp -fPIC -mcmodel=medium
                                $<$<CONFIG:Debug>:-g>
                                $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                        )
set_target_properties(umatlib PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/umatlib)
install(TARGETS umatlib EXPORT umatlib-targets LIBRARY)
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/umatlib/ TYPE INCLUDE)

# Create and install package configuration and version files.
configure_file(
  umatlib-config.cmake.in
  ${CMAKE_BINARY_DIR}/umatlib-config.cmake @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/umatlib-config.cmake
	TYPE LIB
        )

install(EXPORT umatlib-targets
        DESTINATION lib)


add_executable(my_exe main.f90)
target_include_directories(my_exe PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
target_compile_options(my_exe PRIVATE 
                        $<$<CONFIG:Debug>:-g>
                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                        )
target_link_options(my_exe PRIVATE 
                        $<$<CONFIG:Debug>:-g>
                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                        )
target_link_libraries(my_exe PRIVATE umatlib)


