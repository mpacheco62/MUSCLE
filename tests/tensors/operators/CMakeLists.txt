# Every entry of the operator table (op_entries) against the same operation on the full tensors,
# every dimension and kind (D-015.8), one test per generic (.ddot.) or operator module
# (ops_full_common.fypp).
foreach(_name IN ITEMS ddot addition_subtraction assignment contraction_single dyadic transform)
    # The test file is named after the operator module (`.ddot.` is muscle_tensor_ops_contraction_double)
    set(_module ${_name})
    if(_name STREQUAL "ddot")
        set(_module contraction_double)
    endif()
    add_muscle_fypp_test(ops_${_name}_full test_muscle_tensor_ops_${_module}_full.fypp
                         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/ops_full_common.fypp
                                 ${MUSCLE_FYPP_INCLUDE_DIR}/muscle_tensor_ops_common.fypp)
endforeach()
