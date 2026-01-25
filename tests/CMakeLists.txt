add_library(compile_options INTERFACE)
add_custom_target(build_all_tests)

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


function(add_umat_test NAME)
    add_executable(${NAME} ${ARGN})
    target_link_libraries(${NAME} PRIVATE umatlib compile_options)
    target_include_directories(${NAME} PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
    add_dependencies(build_all_tests ${NAME})
    add_test(NAME ${NAME} COMMAND ${NAME})
endfunction()


add_subdirectory(tensors_types)
add_subdirectory(operators)

add_umat_test(test_vonMises test_vonMises.f90)
add_umat_test(test_CPB06 test_CPB06.f90)
add_umat_test(test_swift_law test_swift_law.f90)
add_umat_test(test_bilinear_hardening test_bilinear_hardening.f90)
add_umat_test(test_elasticity_linear test_elasticity_linear.f90)
add_umat_test(test_cutting_plane test_cutting_plane.f90)
add_umat_test(test_closest_point test_closest_point.f90)
add_umat_test(test_closest_point_time test_closest_point_time.f90)
add_umat_test(test_closest_point_2 test_closest_point_2.f90)
target_link_libraries(test_closest_point_2 PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_umat_test(test_basic_operations test_basic_operations.f90)
add_umat_test(test_viscoplastic_models test_viscoplastic_models.f90)