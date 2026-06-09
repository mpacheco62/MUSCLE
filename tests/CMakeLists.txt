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


# Store the root tests directory path
set(BASE_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Base tests directory")

function(add_umat_test NAME)
    # Define executable name with 'test_' prefix
    set(EXE_NAME "test_${NAME}")

    # Create executable target
    add_executable(${EXE_NAME} ${ARGN})
    target_link_libraries(${EXE_NAME} PRIVATE umatlib compile_options)
    target_include_directories(${EXE_NAME} PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
    add_dependencies(build_all_tests ${EXE_NAME})


    # Get relative path from tests root to current folder
    file(RELATIVE_PATH REL_PATH "${BASE_TESTS_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    
    # Set delimiter matching VS Code 'cmake.ctest.testSuiteDelimiter'
    set(DELIMITER "-") 
    
    if (REL_PATH STREQUAL "")
        # No prefix if defined in the root tests folder
        set(TEST_DISPLAY_NAME "${NAME}")
    else()
        # Replace '/' with the delimiter for nested folders
        string(REPLACE "/" "${DELIMITER}" CATEGORY_PATH "${REL_PATH}")
        set(TEST_DISPLAY_NAME "${CATEGORY_PATH}${DELIMITER}${NAME}")
    endif()
    
    # Add test using prefix-free name for CTest and prefixed name for executable
    add_test(NAME ${TEST_DISPLAY_NAME} COMMAND ${EXE_NAME})
    
    target_link_libraries(${EXE_NAME} PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endfunction()


add_subdirectory(tensors_types)
add_subdirectory(math)

add_umat_test(vonMises test_vonMises.f90)
add_umat_test(hill48 test_hill48.f90)
add_umat_test(CPB06 test_CPB06.f90)
add_umat_test(swift_law test_swift_law.f90)
add_umat_test(bilinear_hardening test_bilinear_hardening.f90)
add_umat_test(elasticity_linear test_elasticity_linear.f90)
add_umat_test(cutting_plane test_cutting_plane.f90)
add_umat_test(closest_point test_closest_point.f90)
add_umat_test(closest_point_time test_closest_point_time.f90)
add_umat_test(closest_point_2 test_closest_point_2.f90)
add_umat_test(viscoplastic_models test_viscoplastic_models.f90)
add_umat_test(ludwik_law test_ludwik_law.f90)
