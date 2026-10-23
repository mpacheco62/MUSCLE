muscle_generate_fypp(muscle closest_point_lsvm.fypp muscle_solver_closest_point_lsvm)
muscle_generate_fypp(muscle material_kernels.fypp muscle_kernels_cpu -SVARIANT=cpu)

# Target routines are isolated from the original polymorphic CPU library.
add_library(muscle_target STATIC)
set_target_properties(muscle_target PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/target"
  POSITION_INDEPENDENT_CODE ON)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/target")
target_include_directories(muscle_target PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules/target>
  $<INSTALL_INTERFACE:include>)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
target_link_libraries(muscle_target PUBLIC OpenMP::OpenMP_Fortran)
if(MUSCLE_TARGET_BACKEND STREQUAL "nvidia_gcc12")
  get_filename_component(_driver "${CMAKE_Fortran_COMPILER}" REALPATH)
  get_filename_component(_expected "${PROJECT_SOURCE_DIR}/validation/openmp-fypp/cmake/gfortran12-nvidia.sh" REALPATH)
  if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR
      CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 12 OR
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 13 OR
      NOT _driver STREQUAL _expected)
    message(FATAL_ERROR "nvidia_gcc12 requires the pinned GCC12 wrapper; use the fypp-nvidia-gcc12 preset in a fresh build directory")
  endif()
  target_compile_options(muscle PRIVATE -foffload=disable)
  target_compile_options(muscle_target PUBLIC -foffload=nvptx-none -foffload-options=nvptx-none=-march=sm_80)
  target_link_options(muscle_target PUBLIC -fopenmp -foffload=nvptx-none
    -foffload-options=nvptx-none=-march=sm_80
    -foffload-options=nvptx-none=-lm -foffload-options=nvptx-none=-latomic)
elseif(MUSCLE_TARGET_BACKEND STREQUAL "host")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(muscle_target PUBLIC -foffload=disable)
    target_link_options(muscle_target PUBLIC -foffload=disable)
  endif()
else()
  message(FATAL_ERROR "MUSCLE_TARGET_BACKEND must be host or nvidia_gcc12; AMD is not yet validated")
endif()
muscle_generate_fypp(muscle_target material_kernels.fypp muscle_kernels_target -SVARIANT=target)
install(TARGETS muscle_target EXPORT muscle-targets ARCHIVE DESTINATION lib)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/target/" TYPE INCLUDE)
