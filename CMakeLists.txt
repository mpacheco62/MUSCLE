cmake_minimum_required(VERSION 3.20)

project(Lib_Umats VERSION 0.1
        DESCRIPTION "Library for umats"
        LANGUAGES Fortran C) # C for the tests
include(CTest)
enable_testing()

enable_language(Fortran)


SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/output")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")


# Source code
add_subdirectory(src)
add_subdirectory(tests)

