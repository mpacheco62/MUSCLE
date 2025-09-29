add_executable(test_iden_2OS test_iden_2OS.f90)
target_link_libraries(test_iden_2OS PRIVATE umatlib)
target_include_directories(test_iden_2OS PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_iden_2OS COMMAND test_iden_2OS)

add_executable(test_ten_3d2Osym test_ten_3d2Osym.f90)
target_link_libraries(test_ten_3d2Osym PRIVATE umatlib)
target_include_directories(test_ten_3d2Osym PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_ten_3d2Osym COMMAND test_ten_3d2Osym)

add_executable(test_ten_2d4O3sym test_ten_2d4O3sym.f90)
target_link_libraries(test_ten_2d4O3sym PRIVATE umatlib)
target_include_directories(test_ten_2d4O3sym PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_ten_2d4O3sym COMMAND test_ten_2d4O3sym)


add_executable(test_ten_3d2O test_ten_3d2O.f90)
target_link_libraries(test_ten_3d2O PRIVATE umatlib)
target_include_directories(test_ten_3d2O PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)
add_test(NAME test_ten_3d2O COMMAND test_ten_3d2O)

# set(TEST_LIST
#   test_iden_2OS_sum
# )

# create_test_sourcelist(_ tests_iden_2OS_main.c ${TEST_LIST})
# add_executable(my_tests tests_iden_2OS_main.c test_iden_2OS2.f90)
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
