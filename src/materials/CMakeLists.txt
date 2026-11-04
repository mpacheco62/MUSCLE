# Declared materials (D-008, D-018): one module muscle_mat_<name> per entry of the material registry
# (tools/fypp/muscle_materials.py, or -DMUSCLE_MATERIALS=<path>) from the same template, which
# includes the body of the return mapping of the material (src/solvers/muscle_<solver>_body.fypp,
# one per muscle_meta.MATERIAL_SOLVERS: muscle_add_materials adds them as dependencies) and the
# point kernel and history layout (muscle_material_kernel.fypp) and the batch API
# (muscle_material_batch.fypp) and the point API for FE codes (muscle_material_point.fypp)
muscle_add_materials(muscle muscle_material.fypp
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../solvers/muscle_closest_point_tangent.fypp
                             ${CMAKE_CURRENT_SOURCE_DIR}/muscle_material_kernel.fypp
                             ${CMAKE_CURRENT_SOURCE_DIR}/muscle_material_batch.fypp
                             ${CMAKE_CURRENT_SOURCE_DIR}/muscle_material_point.fypp
                     )
