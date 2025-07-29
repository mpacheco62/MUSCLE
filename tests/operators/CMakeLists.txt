add_executable(test_derivatives test_derivatives.f90)
target_link_libraries(test_derivatives PRIVATE umatlib)
target_include_directories(test_derivatives PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(test_derivatives "${CMAKE_SOURCE_DIR}/bin/test_derivatives")

