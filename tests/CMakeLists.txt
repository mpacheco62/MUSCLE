add_library(compile_options INTERFACE)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        message(STATUS "INTEL Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        target_compile_options(compile_options INTERFACE "/O3" "/warn:all")
        target_compile_options(compile_options INTERFACE
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )
        target_link_options(compile_options INTERFACE
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )

        # target_compile_options(umatlib PRIVATE 
        #                                 -c -Qopenmp -fpic
        #                         )
elseif  (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(STATUS "GNU Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        target_compile_options(compile_options INTERFACE "-O3" "-Wall")
        target_compile_options(compile_options INTERFACE
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile_unittests>:-g -O0 -fprofile-arcs -ftest-coverage>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                )
        target_link_options(compile_options INTERFACE 
                        $<$<CONFIG:Debug>:-g>
                        $<$<CONFIG:Profile_unittests>:--coverage>
                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                        )
else()
    message(FATAL_ERROR "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
endif()


add_subdirectory(tensors_types)
add_subdirectory(operators)

add_executable(test_vonMises test_vonMises.f90)
target_link_libraries(test_vonMises PRIVATE umatlib compile_options)
target_include_directories(test_vonMises PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_vonMises COMMAND test_vonMises)

add_executable(test_CPB06 test_CPB06.f90)
target_link_libraries(test_CPB06 PRIVATE umatlib compile_options)
target_include_directories(test_CPB06 PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_CPB06 COMMAND test_CPB06)

add_executable(test_swift_law test_swift_law.f90)
target_link_libraries(test_swift_law PRIVATE umatlib compile_options)
target_include_directories(test_swift_law PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_swift_law COMMAND test_swift_law)

add_executable(test_bilinear_hardening test_bilinear_hardening.f90)
target_link_libraries(test_bilinear_hardening PRIVATE umatlib compile_options)
target_include_directories(test_bilinear_hardening PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_bilinear_hardening COMMAND test_bilinear_hardening)

add_executable(test_elasticity_linear test_elasticity_linear.f90)
target_link_libraries(test_elasticity_linear PRIVATE umatlib compile_options)
target_include_directories(test_elasticity_linear PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_elasticity_linear COMMAND test_elasticity_linear)

add_executable(test_cutting_plane test_cutting_plane.f90)
target_link_libraries(test_cutting_plane PRIVATE umatlib compile_options)
target_include_directories(test_cutting_plane PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_cutting_plane COMMAND test_cutting_plane)

add_executable(test_closest_point test_closest_point.f90)
target_link_libraries(test_closest_point PRIVATE umatlib compile_options)
target_include_directories(test_closest_point PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_closest_point COMMAND test_closest_point)

add_executable(test_closest_point_time test_closest_point_time.f90)
target_link_libraries(test_closest_point_time PRIVATE umatlib compile_options)
target_include_directories(test_closest_point_time PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_closest_point_time COMMAND test_closest_point_time)

add_executable(test_closest_point_2 test_closest_point_2.f90)
target_link_libraries(test_closest_point_2 PRIVATE umatlib compile_options ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(test_closest_point_2 PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_closest_point_2 COMMAND test_closest_point_2)

add_executable(test_basic_operations test_basic_operations.f90)
target_link_libraries(test_basic_operations PRIVATE umatlib compile_options)
target_include_directories(test_basic_operations PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_basic_operations COMMAND test_basic_operations)

add_executable(test_viscoplastic_models test_viscoplastic_models.f90)
target_link_libraries(test_viscoplastic_models PRIVATE umatlib compile_options)
target_include_directories(test_viscoplastic_models PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_viscoplastic_models COMMAND test_viscoplastic_models)