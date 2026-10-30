target_sources(muscle PRIVATE
                   muscle_tensor_functions.f90
                   )

muscle_add_fypp_sources(muscle
                   types/muscle_tensor_2osym.fypp
                   types/muscle_tensor_2o.fypp
                   types/muscle_tensor_4o3sym.fypp
                   types/muscle_tensor_4o2sym.fypp
                   types/muscle_tensor_4o.fypp
                   types/muscle_tensor_iden.fypp
                   muscle_tensor_io.fypp
                   muscle_tensor_minus.fypp
                   muscle_tensors.fypp
                   )

# Operators between types (D-015): generated from the pair table of tools/fypp/muscle_meta.py
muscle_add_fypp_sources(muscle
                   operators/muscle_tensor_ops_addition_subtraction.fypp
                   operators/muscle_tensor_ops_assignment.fypp
                   operators/muscle_tensor_ops_dyadic.fypp
                   operators/muscle_tensor_ops_contraction_double.fypp
                   operators/muscle_tensor_ops_contraction_single.fypp
                   operators/muscle_tensor_ops_transform.fypp
                   DEPENDS "${MUSCLE_FYPP_INCLUDE_DIR}/muscle_tensor_ops_common.fypp"
                   )
