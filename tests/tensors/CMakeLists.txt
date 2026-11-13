add_subdirectory(operators)
add_subdirectory(types)

add_muscle_fypp_test(functions test_muscle_tensor_functions.fypp DEPENDS ${BASE_TESTS_DIR}/gpu/target_cases.py)
# The free functions of muscle_tensor_functions in an OpenMP target region against the CPU, one test
# per family, dimension and kind (add_muscle_fypp_target_tests, tests/CMakeLists.txt)
add_muscle_fypp_target_tests(tensor_target test_muscle_tensor_target.fypp tensor)
