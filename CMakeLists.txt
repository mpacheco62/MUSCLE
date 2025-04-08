cmake_minimum_required(VERSION 3.20)

project(Lib_Umats VERSION 0.1
        DESCRIPTION "Library for umats"
        LANGUAGES Fortran C) # C for the tests
include(CTest)
enable_testing()

enable_language(Fortran)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${CMAKE_SOURCE_DIR}/output")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")


# Source code
add_subdirectory(src)
add_subdirectory(tests)

