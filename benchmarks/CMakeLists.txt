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

    # Batch API GPU vs CPU (Phase 8, D-019): one executable bench_batch_gpu_<material>_<kind> per case of
    # gpu/bench_cases.py, so that each device image stays small. The template takes the sample inputs of
    # tests/materials/props_samples.py and the device tolerance of muscle_test_utils. They are programs
    # without modules of their own (no module directory per case, unlike the target tests).
    set(_bench_dir "${CMAKE_CURRENT_SOURCE_DIR}/gpu")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1
                ${Python3_EXECUTABLE} "${_bench_dir}/bench_cases.py" "${MUSCLE_REAL_KINDS}"
        OUTPUT_VARIABLE _bench_cases
        RESULT_VARIABLE _bench_result
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT _bench_result EQUAL 0)
        message(FATAL_ERROR "MUSCLE: ${_bench_dir}/bench_cases.py failed")
    endif()
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_bench_dir}/bench_cases.py")
    _muscle_fypp_source("${_bench_dir}/bench_batch_gpu.fypp" _bench_abs _bench_rel)
    foreach(_case IN LISTS _bench_cases)
        set(_out "${MUSCLE_BINARY_DIR}/fypp/benchmarks/gpu/bench_batch_gpu_${_case}.f90")
        _muscle_fypp_command("${_bench_abs}" "${_out}" "FYPP ${_bench_rel} (${_case})"
                             FLAGS "-DBENCH_CASE='${_case}'"
                             DEPENDS "${_bench_dir}/bench_cases.py"
                                     "${BASE_TESTS_DIR}/materials/props_samples.py"
                                     "${BASE_TESTS_DIR}/materials/case_params.py")
        set_source_files_properties("${_out}" PROPERTIES Fortran_PREPROCESS ON)
        add_muscle_benchmark(batch_gpu_${_case} "${_out}")
        target_link_libraries(bench_batch_gpu_${_case} PRIVATE muscle_test_utils)
    endforeach()
endif()
