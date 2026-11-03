# Performance benchmarks of MUSCLE.
#
# They are NOT part of the test suite (not registered in CTest): they report timings instead of
# a pass/fail result. Enable them with -DMUSCLE_BUILD_BENCHMARKS=ON and build the target
# build_all_benchmarks. Use a Release build to obtain meaningful timings.

find_package(OpenMP REQUIRED COMPONENTS Fortran)
add_custom_target(build_all_benchmarks)

function(add_muscle_benchmark NAME)
    set(EXE_NAME "bench_${NAME}")
    add_executable(${EXE_NAME} ${ARGN})
    target_link_libraries(${EXE_NAME} PRIVATE muscle OpenMP::OpenMP_Fortran
                          ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    target_include_directories(${EXE_NAME} PRIVATE ${CMAKE_BINARY_DIR}/modules/muscle)
    # The benchmark prints the build type of the library it is linked against
    set_target_properties(${EXE_NAME} PROPERTIES Fortran_PREPROCESS ON)
    target_compile_definitions(${EXE_NAME} PRIVATE MUSCLE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${EXE_NAME} PRIVATE -O3 -Wall)
    endif()
    add_dependencies(build_all_benchmarks ${EXE_NAME})
endfunction()

# Batch API of materials of the test registry (tools/fypp/muscle_materials.py): not built with a user
# registry (-DMUSCLE_MATERIALS)
if(NOT MUSCLE_MATERIALS)
    add_muscle_benchmark(closest_point_scalability solvers/bench_closest_point_scalability.f90)
endif()
