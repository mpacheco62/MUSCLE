# GCC/NVPTX requires the same driver wrapper for compiler checks, object
# compilation and final links. CMake 3.22 has no Fortran linker launcher.
set(CMAKE_Fortran_COMPILER "${CMAKE_CURRENT_LIST_DIR}/gfortran12-nvidia.sh")
