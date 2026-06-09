
#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_library(umatlib )

set_target_properties(umatlib PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/umatlib)

add_subdirectory(yield_criteria)
add_subdirectory(plasticity)
add_subdirectory(math)
add_subdirectory(tensors)
add_subdirectory(elasticity)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        message(STATUS "INTEL Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        target_compile_options(umatlib PRIVATE 
                                        /c /Qopenmp
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )

        # target_compile_options(umatlib PRIVATE 
        #                                 -c -Qopenmp -fpic
        #                         )
elseif  (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(STATUS "GNU Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")

        target_compile_options(umatlib PRIVATE 
                                        -c -cpp -fopenmp -fPIC -mcmodel=medium
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile_unittests>:-g -O0 -fprofile-arcs -ftest-coverage>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                )

        target_link_options(umatlib PRIVATE
                $<$<CONFIG:Profile_unittests>:--coverage>
                )
else()
    # Default options or a warning if an unsupported compiler is found
    message(FATAL_ERROR "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
endif()


install(TARGETS umatlib EXPORT umatlib-targets LIBRARY)
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/umatlib/ TYPE INCLUDE)

# Create and install package configuration and version files.
configure_file(
  umatlib-config.cmake.in
  ${CMAKE_BINARY_DIR}/umatlib-config.cmake @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/umatlib-config.cmake
	TYPE LIB
        )

install(EXPORT umatlib-targets
        DESTINATION lib)


add_executable(my_exe main.f90)
target_include_directories(my_exe PUBLIC ${CMAKE_BINARY_DIR}/modules/umatlib)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        target_compile_options(my_exe PRIVATE
                                        /c /Qopenmp
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )
        target_link_options(my_exe PRIVATE 
                                        /c /Qopenmp
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
                                )
elseif  (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(my_exe PRIVATE 
                                        -c -cpp -fopenmp -fPIC -mcmodel=medium
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                )
        target_link_options(my_exe PRIVATE 
                                        -c -cpp -fopenmp -fPIC -mcmodel=medium
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                        ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
                                )
else()
    # Default options or a warning if an unsupported compiler is found
    error(ERROR "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
endif()

target_link_libraries(my_exe PRIVATE umatlib)


