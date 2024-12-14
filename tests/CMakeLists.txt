# add_executable(test_basic_operations test_basic_operations.f90)
# target_link_libraries(test_basic_operations PRIVATE mylib)

add_executable(test_vonMises test_vonMises.f90)
target_link_libraries(test_vonMises PRIVATE umatlib)
add_test(Test_vonMises "${CMAKE_SOURCE_DIR}/bin/test_vonMises")


add_subdirectory(tensors_types)

# add_executable(test_tensors_types_ten_3d2Osym tensors_types/test_tensors_types_ten_3d2Osym.f90)
# target_link_libraries(test_tensors_types_ten_3d2Osym PRIVATE mylib)

# add_test(Test_basic_operations "${CMAKE_SOURCE_DIR}/bin/test_basic_operations")
#add_test(Test_tensors_types "${CMAKE_SOURCE_DIR}/bin/test_tensors_types")
# add_test(Test_tensors_types_ten_3d2Osym "${CMAKE_SOURCE_DIR}/bin/test_tensors_types_ten_3d2Osym")