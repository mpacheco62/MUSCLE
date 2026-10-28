find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Unit tests of the FYPP metadata (tools/fypp/muscle_meta.py)
add_test(NAME tools-muscle_meta
         COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_muscle_meta
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
