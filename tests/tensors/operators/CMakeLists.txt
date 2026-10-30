add_muscle_test(test_ddot_4O3_4O3 test_ddot_4O3_4O3.f90)
add_muscle_test(test_operator_contraction_single test_operator_contraction_single.f90)
add_muscle_test(test_operator_dyadic test_operator_dyadic.f90)
add_muscle_test(test_operator_transform test_operator_transform.f90)

add_muscle_test(test_operator_2d_vs_3d test_operator_2d_vs_3d.f90)

# Golden (characterization) tests of every operator between tensor types, double precision,
# one per module of src/tensors/operators/ (ops_golden_cases.fypp; references captured from
# the original code in tests/reference/tensors/operators/, regenerate.sh there).
foreach(_module IN ITEMS addition_subtraction assignment contraction_double contraction_single
                         dyadic transform)
    add_muscle_fypp_test(ops_golden_${_module} test_ops_golden_${_module}.fypp
                         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/ops_golden_cases.fypp)
    set_tests_properties(tensors-operators-ops_golden_${_module} PROPERTIES
                         ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")
endforeach()

# Every entry of the operator table (op_entries) against the same operation on the full tensors,
# every dimension and kind (D-015.8), one test per generic (.ddot.) or operator module
# (ops_full_common.fypp).
foreach(_name IN ITEMS ddot addition_subtraction assignment contraction_single dyadic transform)
    add_muscle_fypp_test(ops_${_name}_full test_ops_${_name}_full.fypp
                         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/ops_full_common.fypp
                                 ${MUSCLE_FYPP_INCLUDE_DIR}/muscle_tensor_ops_common.fypp)
endforeach()
