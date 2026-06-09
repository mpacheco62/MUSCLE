
#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(muscle )

set_target_properties(muscle PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/muscle)

add_subdirectory(plasticity)
add_subdirectory(math)
add_subdirectory(tensors)
add_subdirectory(elasticity)
add_subdirectory(solvers)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        message(STATUS "INTEL Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        target_compile_options(muscle PRIVATE 
                                        /c /Qopenmp
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )

        # target_compile_options(muscle PRIVATE 
        #                                 -c -Qopenmp -fpic
        #                         )
elseif  (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(STATUS "GNU Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")

        target_compile_options(muscle PRIVATE 
                                        -c -cpp -fopenmp -fPIC -mcmodel=medium
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile_unittests>:-g -O0 -fprofile-arcs -ftest-coverage>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                )

        target_link_options(muscle PRIVATE
                $<$<CONFIG:Profile_unittests>:--coverage>
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



