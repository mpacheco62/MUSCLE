add_executable(test_iden_3d2Omod test_iden_3d2Omod.f90)
target_link_libraries(test_iden_3d2Omod PRIVATE umatlib)
target_include_directories(test_iden_3d2Omod PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(Test_iden_3d2Omod "${CMAKE_SOURCE_DIR}/bin/test_iden_3d2Omod")

add_executable(test_ten_3d2Osym test_ten_3d2Osym.f90)
target_link_libraries(test_ten_3d2Osym PRIVATE umatlib)
target_include_directories(test_ten_3d2Osym PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(Test_ten_3d2Osym "${CMAKE_SOURCE_DIR}/bin/test_ten_3d2Osym")

add_executable(test_ten_3d2O test_ten_3d2O.f90)
target_link_libraries(test_ten_3d2O PRIVATE umatlib)
target_include_directories(test_ten_3d2O PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(Test_ten_3d2O "${CMAKE_SOURCE_DIR}/bin/test_ten_3d2O")

# set(TEST_LIST
#   test_iden_3D2OMod_sum
# )

# create_test_sourcelist(_ tests_iden_3D2OMod_main.c ${TEST_LIST})
# add_executable(my_tests tests_iden_3D2OMod_main.c test_iden_3d2Omod2.f90)
# target_link_libraries(my_tests PRIVATE mylib)

# foreach (test ${TEST_LIST})
#   add_test(NAME ${test} COMMAND my_tests ${test})
# endforeach ()


# set(TEST_LIST
#   my_test_success
#   my_test_fail
# )
# create_test_sourcelist(_ my_tests_main.c ${TEST_LIST})
# add_executable(my_tests my_tests_main.c my_test.f90)

# foreach (test ${TEST_LIST})
#   add_test(NAME ${test} COMMAND my_tests ${test})
# endforeach ()
