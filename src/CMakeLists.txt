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

<<<<<<< HEAD
add_subdirectory(plasticity)
=======
#set(HARDENING_LAWS   hardening_laws/basis_hardening_law.f90
#                     hardening_laws/swift_hardening.f90
#                     hardening_laws/bilinear_hardening.f90
#                     CACHE INTERNAL "")

set(HARDENING_LAWS   plasticity/hardening_laws/basis_hardening_law.f90
                     plasticity/hardening_laws/swift_hardening.f90
                     plasticity/hardening_laws/bilinear_hardening.f90
                     CACHE INTERNAL "")

set(VISCO_HARDENING_LAWS   plasticity/viscohardening_laws/basis_visco_hardening_law.f90
                     	   plasticity/viscohardening_laws/JohnsonCook_full_visco_hardening.f90
                     	   plasticity/viscohardening_laws/JohnsonCook_visco_hardening.f90
                     CACHE INTERNAL "")
>>>>>>> df8a4a6 (Reorganizacion de carpetas)

set(YIELD_CRITERIA   yield_criteria/basis_yield_criteria.f90
                     yield_criteria/vonMises.f90
                     yield_criteria/CPB06.f90
                     CACHE INTERNAL "")

set(ELASTICITY       elasticity/basis_elasticity.f90
                     elasticity/elasticity_linear.f90
                     CACHE INTERNAL "")

<<<<<<< HEAD
set(OPERATORS   operators/basic_operations.f90
                operators/inverses_mat.f90
                operators/derivatives.f90
                operators/lin_solvers.f90
=======
#set(RETURN_MAPPING   return_mapping/cutting_plane.f90
#                     return_mapping/closest_point.f90
#                     CACHE INTERNAL "")

set(RETURN_MAPPING   plasticity/return_mapping/cutting_plane.f90
                     plasticity/return_mapping/closest_point.f90
		     plasticity/return_mapping/closest_point2.f90
                     CACHE INTERNAL "")

set(OPERATORS   operators/basic_operations.f90
                operators/inverses_mat.f90
                operators/derivatives.f90
		operators/lin_solvers.f90
>>>>>>> df8a4a6 (Reorganizacion de carpetas)
               CACHE INTERNAL "")

      #set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(umatlib 
                     ${TENSORS_TYPES}
                     ${HARDENING_LAWS}
                     ${YIELD_CRITERIA}
                     ${ELASTICITY}
                     ${RETURN_MAPPING}
                     ${OPERATORS}
                     ${VISCOHARDENING_LAWS}
                     )

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        message(STATUS "INTEL Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        target_compile_options(umatlib PRIVATE 
                                        /c /Qopenmp
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )

        # target_compile_options(umatlib PRIVATE 
        #                                 -c -Qopenmp -fpic
        #                         )
elseif  (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(STATUS "GNU Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        target_compile_options(umatlib PRIVATE 
                                        -c -cpp -fopenmp -fPIC -mcmodel=medium
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                )
else()
    # Default options or a warning if an unsupported compiler is found
    error(STATUS "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
endif()

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

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        target_compile_options(my_exe PRIVATE
                                        /c /Qopenmp
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )
        target_link_options(my_exe PRIVATE 
                                        /c /Qopenmp
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
                                )
elseif  (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(my_exe PRIVATE 
                                        -c -cpp -fopenmp -fPIC -mcmodel=medium
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                )
        target_link_options(my_exe PRIVATE 
                                        -c -cpp -fopenmp -fPIC -mcmodel=medium
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
                                )
else()
    # Default options or a warning if an unsupported compiler is found
    error(ERROR "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
endif()

target_link_libraries(my_exe PRIVATE umatlib)


