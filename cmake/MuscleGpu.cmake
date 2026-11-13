# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# GPU build of MUSCLE (OpenMP offload, D-005, D-019).
#
#   MUSCLE_GPU       OFF (default), NVPTX, AMDGCN or SPIR64
#   MUSCLE_GPU_ARCH  device architecture; empty: the default of the target (NVPTX: the compiler
#                    default, portable PTX compiled at load time; AMDGCN: gfx90a; SPIR64: JIT)
#
# With a target other than OFF:
# - the fypp define MUSCLE_GPU is appended to MUSCLE_FYPP_DEFINES: @:DEVICE() expands to
#   `!$omp declare target` and BATCH_LOOP to a target region (include/muscle_common.fypp);
# - MUSCLE_GPU_COMPILE_OPTIONS / MUSCLE_GPU_LINK_OPTIONS hold the offload flags. The library adds
#   them as PUBLIC options (src/CMakeLists.txt), so the tests and benchmarks that link `muscle`
#   compile and link with the same flags;
# - MUSCLE_GPU_ENABLED (ON/OFF, cache) tells every directory whether offload is on.
#
#   target  compiler          flags
#   NVPTX   gfortran (GNU)    -fopenmp -foffload=nvptx-none [-foffload-options=nvptx-none=-march=<arch>]
#   AMDGCN  gfortran (GNU)    -fopenmp -foffload=amdgcn-amdhsa -foffload-options=amdgcn-amdhsa=-march=<arch>
#   SPIR64  ifx (IntelLLVM)   -fiopenmp -fopenmp-targets=spir64
#   (gfortran links libgomp with --no-as-needed too, see below)
#
# The combination is validated at configure time (compiler, installed offload backend) by linking a
# small program with a target region and a `declare target` procedure.

include_guard(GLOBAL)

set(MUSCLE_GPU "OFF" CACHE STRING "OpenMP offload target of the kernels: OFF, NVPTX, AMDGCN or SPIR64")
set_property(CACHE MUSCLE_GPU PROPERTY STRINGS OFF NVPTX AMDGCN SPIR64)
set(MUSCLE_GPU_ARCH "" CACHE STRING
    "Device architecture for MUSCLE_GPU (e.g. sm_80, gfx90a); empty: the default of the target")
string(TOUPPER "${MUSCLE_GPU}" _muscle_gpu)

set(MUSCLE_GPU_COMPILE_OPTIONS)
set(MUSCLE_GPU_LINK_OPTIONS)
# MUSCLE_GPU_ENABLED: ON/OFF for every directory (tests/ too), whatever the case of MUSCLE_GPU
set(MUSCLE_GPU_ENABLED OFF CACHE INTERNAL "MUSCLE_GPU is a device target (not OFF)")

if(_muscle_gpu STREQUAL "OFF" OR _muscle_gpu STREQUAL "")
    message(STATUS "MUSCLE: GPU offload disabled (MUSCLE_GPU=OFF)")
    return()
endif()

set(_muscle_gpu_arch "${MUSCLE_GPU_ARCH}")
if(_muscle_gpu STREQUAL "NVPTX" OR _muscle_gpu STREQUAL "AMDGCN")
    if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "MUSCLE: MUSCLE_GPU=${MUSCLE_GPU} needs gfortran (GNU); the Fortran compiler is "
                            "${CMAKE_Fortran_COMPILER_ID} (${CMAKE_Fortran_COMPILER}). Use MUSCLE_GPU=SPIR64 with ifx.")
    endif()
    if(_muscle_gpu STREQUAL "NVPTX")
        set(_muscle_gpu_triple "nvptx-none")
    else()
        set(_muscle_gpu_triple "amdgcn-amdhsa")
        if(_muscle_gpu_arch STREQUAL "")
            set(_muscle_gpu_arch "gfx90a")
        endif()
    endif()
    # The offload targets gfortran was configured with (`gfortran -v`: OFFLOAD_TARGET_NAMES=a:b)
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -v ERROR_VARIABLE _muscle_gpu_v OUTPUT_QUIET)
    string(REGEX MATCH "OFFLOAD_TARGET_NAMES=([^ \n]*)" _muscle_gpu_names "${_muscle_gpu_v}")
    string(REPLACE ":" ";" _muscle_gpu_names "${CMAKE_MATCH_1}")
    if(NOT _muscle_gpu_triple IN_LIST _muscle_gpu_names)
        message(FATAL_ERROR "MUSCLE: MUSCLE_GPU=${MUSCLE_GPU}: ${CMAKE_Fortran_COMPILER} was not configured for "
                            "offloading to ${_muscle_gpu_triple} (offload targets: '${CMAKE_MATCH_1}').")
    endif()
    set(MUSCLE_GPU_COMPILE_OPTIONS -fopenmp -foffload=${_muscle_gpu_triple})
    if(NOT _muscle_gpu_arch STREQUAL "")
        list(APPEND MUSCLE_GPU_COMPILE_OPTIONS -foffload-options=${_muscle_gpu_triple}=-march=${_muscle_gpu_arch})
    endif()
    # The device code is compiled at link time (LTO): the link needs the same flags. The offload
    # registration object (GOMP_offload_register_ver) comes after every library on the link line:
    # with --as-needed (default of Debian/Ubuntu) libgomp is dropped from a program that calls no
    # OpenMP routine itself, and the link fails. Keep it.
    set(MUSCLE_GPU_LINK_OPTIONS ${MUSCLE_GPU_COMPILE_OPTIONS}
        "LINKER:--push-state,--no-as-needed,-lgomp,--pop-state")
    set(_muscle_gpu_hint "install the offload backend (Debian/Ubuntu: gcc-<version>-offload-nvptx or "
                         "gcc-<version>-offload-amdgcn) or check MUSCLE_GPU_ARCH")
elseif(_muscle_gpu STREQUAL "SPIR64")
    if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        message(FATAL_ERROR "MUSCLE: MUSCLE_GPU=SPIR64 needs ifx (IntelLLVM); the Fortran compiler is "
                            "${CMAKE_Fortran_COMPILER_ID} (${CMAKE_Fortran_COMPILER}). Use NVPTX or AMDGCN with gfortran.")
    endif()
    if(WIN32)
        message(FATAL_ERROR "MUSCLE: MUSCLE_GPU=SPIR64 is only set up for ifx on Linux.")
    endif()
    if(NOT _muscle_gpu_arch STREQUAL "")
        message(FATAL_ERROR "MUSCLE: MUSCLE_GPU_ARCH is not used with SPIR64 (JIT compilation); leave it empty.")
    endif()
    set(MUSCLE_GPU_COMPILE_OPTIONS -fiopenmp -fopenmp-targets=spir64)
    set(MUSCLE_GPU_LINK_OPTIONS ${MUSCLE_GPU_COMPILE_OPTIONS})
    set(_muscle_gpu_hint "check the oneAPI environment (source setvars.sh)")
else()
    message(FATAL_ERROR "MUSCLE: unknown MUSCLE_GPU='${MUSCLE_GPU}' (allowed: OFF, NVPTX, AMDGCN, SPIR64).")
endif()

# Link check of the combination: a declare target procedure called in a target region
set(_muscle_gpu_dir "${CMAKE_BINARY_DIR}/CMakeFiles/muscle_gpu_check")
file(WRITE "${_muscle_gpu_dir}/check.f90" [=[
module muscle_gpu_check_mod
contains
    pure function twice(x) result(y)
        !$omp declare target
        real(kind(1.0d0)), intent(in) :: x
        real(kind(1.0d0)) :: y
        y = 2*x
    end function twice
end module muscle_gpu_check_mod
program muscle_gpu_check
    use muscle_gpu_check_mod, only: twice
    real(kind(1.0d0)) :: a(4)
    integer :: i
    a = 1
    !$omp target teams distribute parallel do map(tofrom: a)
    do i = 1, 4
        a(i) = twice(a(i))
    end do
    print *, a
end program muscle_gpu_check
]=])
try_compile(_muscle_gpu_ok "${_muscle_gpu_dir}/build"
            SOURCES "${_muscle_gpu_dir}/check.f90"
            COMPILE_DEFINITIONS ${MUSCLE_GPU_COMPILE_OPTIONS}
            LINK_OPTIONS ${MUSCLE_GPU_LINK_OPTIONS}
            OUTPUT_VARIABLE _muscle_gpu_out)
if(NOT _muscle_gpu_ok)
    message(FATAL_ERROR "MUSCLE: MUSCLE_GPU=${MUSCLE_GPU}: cannot link an offload program with "
                        "'${MUSCLE_GPU_COMPILE_OPTIONS}' (${_muscle_gpu_hint}).\n${_muscle_gpu_out}")
endif()

list(APPEND MUSCLE_FYPP_DEFINES MUSCLE_GPU)
set(MUSCLE_GPU_ENABLED ON CACHE INTERNAL "MUSCLE_GPU is a device target (not OFF)")
message(STATUS "MUSCLE: GPU offload MUSCLE_GPU=${MUSCLE_GPU}: ${MUSCLE_GPU_COMPILE_OPTIONS}")
