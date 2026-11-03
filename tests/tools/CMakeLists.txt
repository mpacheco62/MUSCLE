find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Unit tests of the FYPP metadata (tools/fypp/muscle_meta.py)
add_test(NAME tools-muscle_meta
         COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_muscle_meta
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Coverage of the golden operator tests (tests/tensors/operators/ops_golden_cases.fypp) against
# the inventory of historical double precision operator signatures (tools/fypp/ref/)
add_test(NAME tools-ops_golden_coverage
         COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_ops_golden_coverage
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Registry of declared materials (tools/fypp/muscle_materials.py) and its validation in muscle_meta
add_test(NAME tools-muscle_materials
         COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_muscle_materials
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Parameters and initializers of the registered models (MODEL_PARAMS, MODEL_INITS of muscle_meta,
# D-018.17e)
add_test(NAME tools-muscle_model_params
         COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_muscle_model_params
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Parameter array of the declared materials (material_props, `muscle_meta.py props NAME`, D-018.17a)
add_test(NAME tools-muscle_material_props
         COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_muscle_material_props
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
