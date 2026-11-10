
#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(muscle )

# Detectar si el compilador soporta UDTIO (write(formatted)) de forma thread-safe en OpenMP.
# MUSCLE_UDTIO=ON/OFF fuerza la decisión (tools/check.sh compila ese camino con ON).
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

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        message(STATUS "INTEL Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        if (WIN32)
            target_compile_options(muscle PRIVATE
                                        /c /Qopenmp
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )
        else()
            # Linux: the Windows-style options above are read as file names
            target_compile_options(muscle PRIVATE
                                        -fpp -qopenmp -fPIC
                                        $<$<CONFIG:Debug>:-g -O0>
                                        $<$<CONFIG:Profile>:-Ofast -ipo -fp-model=fast -DNDEBUG>
                                )
        endif()

        # target_compile_options(muscle PRIVATE 
        #                                 -c -Qopenmp -fpic
        #                         )
elseif  (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(STATUS "GNU Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")

        target_compile_options(muscle PRIVATE 
                                        -c -cpp -fopenmp -fPIC -mcmodel=medium
                                        $<$<CONFIG:Debug>:-g -O0 -Wall -Wextra -Wimplicit-interface -Wsurprising -fcheck=all -fbacktrace -ffpe-trap=zero,overflow,invalid>
                                        # Unused dummy arguments are the uniform contract of a family (a kinematics without
                                        # parameters still takes prm); unused local variables still warn. The exact real
                                        # comparisons are deliberate guards (exact zero before a division, ...): silenced for
                                        # now, to be removed (docs/fypp/PROBLEMAS.md)
                                        $<$<CONFIG:Debug>:-Wno-unused-dummy-argument -Wno-compare-reals>
                                        $<$<CONFIG:Profile_unittests>:-g -O0 -fprofile-arcs -ftest-coverage>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                )

        target_link_options(muscle PRIVATE
                $<$<CONFIG:Profile_unittests>:--coverage>
                $<$<CONFIG:Debug>:-g -fbacktrace>
                )
else()
    # Default options or a warning if an unsupported compiler is found
    message(FATAL_ERROR "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
endif()


install(TARGETS muscle EXPORT muscle-targets LIBRARY)
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/muscle/ TYPE INCLUDE)

# Create and install package configuration and version files.
configure_file(
  muscle-config.cmake.in
  ${CMAKE_BINARY_DIR}/muscle-config.cmake @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/muscle-config.cmake
	TYPE LIB
        )

install(EXPORT muscle-targets
        DESTINATION lib)



