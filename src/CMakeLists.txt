
#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(muscle )

# Detect whether the compiler supports UDTIO (write(formatted)) safely inside OpenMP regions.
# MUSCLE_UDTIO=ON/OFF forces the choice (tools/check.sh builds that path with ON).
set(MUSCLE_UDTIO "AUTO" CACHE STRING "Formatted output of tensors (UDTIO): AUTO, ON or OFF")
set_property(CACHE MUSCLE_UDTIO PROPERTY STRINGS AUTO ON OFF)
if(MUSCLE_UDTIO STREQUAL "AUTO")
    if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "16.3")
        set(_muscle_udtio ON)
    else()
        set(_muscle_udtio OFF)
    endif()
else()
    set(_muscle_udtio ${MUSCLE_UDTIO})
endif()
if(_muscle_udtio)
    message(STATUS "MUSCLE: Habilitando soporte UDTIO write(formatted) (MUSCLE_UDTIO=${MUSCLE_UDTIO})")
    target_compile_definitions(muscle PUBLIC ENABLE_UDTIO)
    list(APPEND MUSCLE_FYPP_DEFINES ENABLE_UDTIO)
else()
    message(STATUS "MUSCLE: Desactivando UDTIO write(formatted) por bug conocido de OpenMP en GFortran < 16.3")
endif()

target_link_libraries(muscle PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
# The batch API calls the OpenMP runtime: a program that links the library also links OpenMP (the library itself
# compiles with the OpenMP flag of MUSCLE_OPT_LIBRARY, cmake/MuscleCompiler.cmake)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
target_link_libraries(muscle INTERFACE OpenMP::OpenMP_Fortran)

if(BLAS_LIBRARIES MATCHES "mkl" OR LAPACK_LIBRARIES MATCHES "mkl")
    message(STATUS "MUSCLE: Intel MKL detectado. Activando directiva USE_MKL.")
    list(APPEND MUSCLE_FYPP_DEFINES USE_MKL)
endif()
# GPU offload (MUSCLE_GPU, cmake/MuscleGpu.cmake): fypp define MUSCLE_GPU and the offload flags, PUBLIC
# so that the tests and benchmarks compile and link with the same flags
include(MuscleGpu)
target_compile_options(muscle PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${MUSCLE_GPU_COMPILE_OPTIONS}>)
target_link_options(muscle PUBLIC ${MUSCLE_GPU_LINK_OPTIONS})
# Device image of the whole library for tools/audit_kernels.sh (gfortran: the offload compiler builds
# the device code at link time; -save-temps keeps its assembly next to the probe, in audit/)
if(MUSCLE_GPU_ENABLED AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_executable(muscle_audit_probe EXCLUDE_FROM_ALL ${PROJECT_SOURCE_DIR}/tools/audit_kernels_probe.f90)
    target_link_libraries(muscle_audit_probe PRIVATE muscle)
    target_link_options(muscle_audit_probe PRIVATE
                        "LINKER:--whole-archive" "$<TARGET_FILE:muscle>" "LINKER:--no-whole-archive"
                        # absolute: the AMDGCN mkoffload embeds its image by a path relative to the temporaries
                        -save-temps=obj -dumpdir ${CMAKE_BINARY_DIR}/audit/muscle_audit_probe.)
    set_target_properties(muscle_audit_probe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/audit)
endif()
# The tests (sibling directory) preprocess their templates with the same switches
set(MUSCLE_FYPP_DEFINES "${MUSCLE_FYPP_DEFINES}" PARENT_SCOPE)

set_target_properties(muscle PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/muscle)
target_include_directories(muscle PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules/muscle>
    $<INSTALL_INTERFACE:include>
)

add_subdirectory(common)
add_subdirectory(plasticity)
add_subdirectory(math)
add_subdirectory(tensors)
add_subdirectory(kinematics)
add_subdirectory(elasticity)
add_subdirectory(viscoelasticity)
add_subdirectory(materials)

# Flags of the library by compiler and OS: cmake/MuscleCompiler.cmake
target_compile_options(muscle PRIVATE ${MUSCLE_OPT_LIBRARY})
target_link_options(muscle PRIVATE ${MUSCLE_LINK_LIBRARY})


# Installation: the library in lib/, the Fortran modules in include/ and the CMake package in lib/cmake/muscle/,
# where find_package(muscle) looks for it (CMAKE_PREFIX_PATH=<prefix>)
install(TARGETS muscle EXPORT muscle-targets)
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/muscle/ TYPE INCLUDE)

include(CMakePackageConfigHelpers)
configure_file(muscle-config.cmake.in ${CMAKE_BINARY_DIR}/muscle-config.cmake @ONLY)
write_basic_package_version_file(${CMAKE_BINARY_DIR}/muscle-config-version.cmake
                                 VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
install(FILES ${CMAKE_BINARY_DIR}/muscle-config.cmake ${CMAKE_BINARY_DIR}/muscle-config-version.cmake
        DESTINATION lib/cmake/muscle)
install(EXPORT muscle-targets DESTINATION lib/cmake/muscle)
