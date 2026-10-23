add_executable(muscle_reference_cpu "${PROJECT_SOURCE_DIR}/validation/reference/cpu_reference.f90")
target_link_libraries(muscle_reference_cpu PRIVATE muscle compile_options)
add_executable(muscle_reference_concrete)
set(_reference "${PROJECT_SOURCE_DIR}/validation/reference/cpu_reference.f90")
configure_file(reference_concrete.fypp.in "${CMAKE_CURRENT_BINARY_DIR}/reference_concrete.fypp" @ONLY)
muscle_generate_fypp(muscle_reference_concrete "${CMAKE_CURRENT_BINARY_DIR}/reference_concrete.fypp" reference_concrete)
target_link_libraries(muscle_reference_concrete PRIVATE muscle compile_options)
add_test(NAME generated_cpu_reference
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/validation/scripts/compare_reference.py"
    --reference "$<TARGET_FILE:muscle_reference_cpu>"
    --concrete "$<TARGET_FILE:muscle_reference_concrete>"
    --baseline "${PROJECT_SOURCE_DIR}/validation/reference/reference-release.csv"
    --output "${CMAKE_CURRENT_BINARY_DIR}/reference-comparison")
set_tests_properties(generated_cpu_reference PROPERTIES
  ENVIRONMENT "OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;OMP_TARGET_OFFLOAD=DISABLED" TIMEOUT 60)

add_executable(muscle_test_generated_kernels test_generated_kernels.f90)
target_link_libraries(muscle_test_generated_kernels PRIVATE muscle muscle_target compile_options)
set(_policy DISABLED)
set(_argument --host)
if(MUSCLE_TARGET_BACKEND STREQUAL "nvidia_gcc12")
  set(_policy MANDATORY)
  set(_argument --require-device)
endif()
add_test(NAME generated_kernels COMMAND muscle_test_generated_kernels "${_argument}")
set_tests_properties(generated_kernels PROPERTIES
  ENVIRONMENT "OMP_NUM_THREADS=4;OMP_DYNAMIC=FALSE;OPENBLAS_NUM_THREADS=1;OMP_TARGET_OFFLOAD=${_policy}" TIMEOUT 60)
add_dependencies(build_all_tests muscle_reference_cpu muscle_reference_concrete muscle_test_generated_kernels)
target_compile_definitions(muscle_reference_concrete PRIVATE MUSCLE_CONCRETE_REFERENCE)
set_source_files_properties("${PROJECT_SOURCE_DIR}/validation/reference/cpu_reference.f90" PROPERTIES Fortran_PREPROCESS ON)
set_target_properties(muscle_reference_cpu muscle_reference_concrete PROPERTIES Fortran_PREPROCESS ON)
