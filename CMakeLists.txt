cmake_minimum_required(VERSION 3.20)
# set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CXX_FLAGS_RELEASE};profiling" CACHE STRING "Available build types.")
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;Profile" CACHE STRING "Available build types.")

project(Lib_Umats VERSION 0.1
        DESCRIPTION "Library for umats"
        LANGUAGES Fortran C) # C for the tests
include(CTest)
enable_testing()

enable_language(Fortran)


# This generate a file compile_commands.json with commands use on compile
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${CMAKE_SOURCE_DIR}/output")
endif()

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")

# Source code
add_subdirectory(src)
add_subdirectory(tests)

