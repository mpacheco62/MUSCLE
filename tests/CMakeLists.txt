# add_executable(test_basic_operations test_basic_operations.f90)
# target_link_libraries(test_basic_operations PRIVATE mylib)
add_subdirectory(tensors_types)

add_executable(test_vonMises test_vonMises.f90)
target_link_libraries(test_vonMises PRIVATE umatlib)
add_test(Test_vonMises "${CMAKE_SOURCE_DIR}/bin/test_vonMises")


add_executable(test_swift_law test_swift_law.f90)
target_link_libraries(test_swift_law PRIVATE umatlib)
add_test(Test_swift_law "${CMAKE_SOURCE_DIR}/bin/test_swift_law")

add_executable(test_elasticity_linear test_elasticity_linear.f90)
target_link_libraries(test_elasticity_linear PRIVATE umatlib)
add_test(Test_elasticity_linear "${CMAKE_SOURCE_DIR}/bin/test_elasticity_linear")

add_executable(test_cutting_plane test_cutting_plane.f90)
target_link_libraries(test_cutting_plane PRIVATE umatlib)
add_test(test_cutting_plane "${CMAKE_SOURCE_DIR}/bin/test_cutting_plane")

add_executable(test_closest_point test_closest_point.f90)
target_link_libraries(test_closest_point PRIVATE umatlib)
add_test(test_closest_point "${CMAKE_SOURCE_DIR}/bin/test_closest_point")
# add_executable(test_tensors_types_ten_3d2Osym tensors_types/test_tensors_types_ten_3d2Osym.f90)
# target_link_libraries(test_tensors_types_ten_3d2Osym PRIVATE mylib)

# add_test(Test_basic_operations "${CMAKE_SOURCE_DIR}/bin/test_basic_operations")
#add_test(Test_tensors_types "${CMAKE_SOURCE_DIR}/bin/test_tensors_types")
# add_test(Test_tensors_types_ten_3d2Osym "${CMAKE_SOURCE_DIR}/bin/test_tensors_types_ten_3d2Osym")