# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Compiler and OS dependent flags of MUSCLE, in one place (D-040).
#
# The library (src/), the tests (tests/) and the benchmarks (benchmarks/) take every flag that depends on the
# Fortran compiler from the variables below; no other CMakeLists.txt asks for CMAKE_Fortran_COMPILER_ID, except
# the GPU offload (cmake/MuscleGpu.cmake) and the build checks of src/CMakeLists.txt (UDTIO, audit probe).
#
# What every compiler has to give (the unit tests assume it):
#   - strict floating point evaluation (no reassociation or contraction) in the code of the tests: gfortran does it
#     by default; ifx asks for the precise model, but in the Profile build, which wants speed. The library compiled by
#     ifx uses its fast model (D-042; it may flush the subnormal numbers to zero: no branch of a model depends on
#     them since D-043);
#     gfortran's is strict
#   - the standard sign of a negative zero: gfortran does it by default; ifx asks for /assume:minus0
#   - the diagnostics of the Debug library: -Wall -Wextra and run time checks of bounds (gfortran: -fcheck=all;
#     ifx: bounds, pointers and shape; `-check all` of ifx is much slower and finds nothing more here)
#   - traps of invalid operations, divisions by zero and overflows in the main program of the tests that ask
#     for them (muscle_fpe_traps; with ifx only in Debug, D-042)
#   - optimization -O3 in the tests, in every configuration, and in the Release library (with ifx on Windows CMake
#     would give /O2)
# A build with GPU offload (MUSCLE_GPU) keeps the ifx flags it had before D-040 (see _muscle_gpu_build below).
#
# Variables (lists; an element may be a generator expression with several options):
#   MUSCLE_OPT_LIBRARY / MUSCLE_LINK_LIBRARY   compile and link options of the library (PRIVATE)
#   MUSCLE_OPT_TESTS / MUSCLE_LINK_TESTS       compile and link options of the tests (tests/CMakeLists.txt)
#   MUSCLE_OPT_BENCH                           compile options of the benchmarks
#   MUSCLE_FLAG_NO_INLINE                      the option that keeps one copy of each procedure
#   MUSCLE_FLAGS_HOST_OFFLOAD                  the options that make omp_target_* work on the host (ifx)
#   MUSCLE_LINK_STACK                          the stack of the main thread of the tests (Windows)
#   MUSCLE_FLAG_FPE_TRAPS                      the options of muscle_fpe_traps
#   MUSCLE_FLAGS_COMPILE_FAIL, MUSCLE_FLAG_INCLUDE, MUSCLE_COMPILE_FAIL_PASS / _FAIL
#                                              the compiler, the include option and the messages of the
#                                              compile-fail tests (muscle_add_compile_fail_test)
#
# Measured on Windows (gfortran 16.2, ifx 2026.2) and on Linux (WSL: gfortran 15.2, ifx 2026.2); see
# docs/fypp/DECISIONES.md D-037 to D-042 and PROBLEMAS.md.

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(MUSCLE_COMPILER "gnu")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    set(MUSCLE_COMPILER "intel")
else()
    message(FATAL_ERROR "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
endif()
message(STATUS "MUSCLE: Fortran compiler ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (flags: cmake/MuscleCompiler.cmake)")

set(MUSCLE_OPT_BENCH "")
set(MUSCLE_LINK_STACK "")

# A build with GPU offload (MUSCLE_GPU, cmake/MuscleGpu.cmake) keeps the ifx flags it had before D-040: the precise
# floating point model, the signed zero and the Debug diagnostics were measured on the CPU only, and the tests of a GPU
# build compare the results of the device with those of the host bit for bit
if(DEFINED MUSCLE_GPU AND NOT MUSCLE_GPU STREQUAL "" AND NOT MUSCLE_GPU STREQUAL "OFF")
    set(_muscle_gpu_build ON)
else()
    set(_muscle_gpu_build OFF)
endif()

if(MUSCLE_COMPILER STREQUAL "gnu")
    set(MUSCLE_OPT_LIBRARY
        -c -cpp -fopenmp -fPIC -mcmodel=medium
        $<$<CONFIG:Debug>:-g -O0 -Wall -Wextra -Wimplicit-interface -Wsurprising -fcheck=all -fbacktrace -ffpe-trap=zero,overflow,invalid>
        # Unused dummy arguments are the uniform contract of a family (a kinematics without parameters still takes
        # prm); unused local variables still warn. The exact real comparisons are deliberate guards (exact zero before
        # a division, ...): silenced for now, to be removed (docs/fypp/PROBLEMAS.md)
        $<$<CONFIG:Debug>:-Wno-unused-dummy-argument -Wno-compare-reals>
        $<$<CONFIG:Profile_unittests>:-g -O0 -fprofile-arcs -ftest-coverage>
        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>)
    set(MUSCLE_LINK_LIBRARY
        $<$<CONFIG:Profile_unittests>:--coverage>
        $<$<CONFIG:Debug>:-g -fbacktrace>)
    # Unused dummy arguments are part of uniform interfaces; unused local variables still warn
    set(MUSCLE_OPT_TESTS
        -O3 -Wall -Wno-unused-dummy-argument
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Profile_unittests>:-g -O0 -fprofile-arcs -ftest-coverage>
        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>)
    set(MUSCLE_LINK_TESTS
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Profile_unittests>:--coverage>
        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>)
    set(MUSCLE_OPT_BENCH -O3 -Wall)
    set(MUSCLE_FLAG_NO_INLINE -fno-inline)
    set(MUSCLE_FLAGS_HOST_OFFLOAD "")
    set(MUSCLE_FLAG_FPE_TRAPS -ffpe-trap=zero,overflow,invalid)
    set(MUSCLE_FLAGS_COMPILE_FAIL -fsyntax-only)
    set(MUSCLE_FLAG_INCLUDE -I)
    set(MUSCLE_COMPILE_FAIL_PASS "PRIVATE component")
    set(MUSCLE_COMPILE_FAIL_FAIL "Cannot open module|Fatal Error")
    if(WIN32)
        # The main thread of a Windows program has a 1 MB stack (8 MB on Linux) and mat_batch_openmp needs more
        set(MUSCLE_LINK_STACK -Wl,--stack,16777216)
    endif()
elseif(WIN32)
    # ifx on Windows: the options are /option, not -option. Release optimizes at /O3, as gfortran and ifx on Linux
    # (the default of CMake is /O2 here)
    string(REPLACE "/O2" "/O3" CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")
    if(_muscle_gpu_build)
        set(_muscle_fp_library "")
        set(_muscle_fp_tests "")
        set(_muscle_diagnostics "")
    else()
        # The library uses the fast floating point model of ifx (D-042; the subnormal numbers may be flushed to zero,
        # no branch of a model depends on them since D-043) and, in Debug, safe speculation (the traps of the Debug tests must not
        # stop on an operation the fast model evaluates speculatively). The tests keep the precise model: the fast
        # one assumes there are no infinities and folds a test against +Inf or NaN to false.
        set(_muscle_fp_library
            /assume:minus0
            $<$<CONFIG:Debug>:/Qfp-speculation:safe>)
        set(_muscle_fp_tests
            /assume:minus0
            $<$<NOT:$<CONFIG:Profile>>:/fp:precise>)
        set(_muscle_diagnostics $<$<CONFIG:Debug>:/warn:all /check:bounds /check:pointers /check:shape /traceback>)
    endif()
    set(MUSCLE_OPT_LIBRARY
        /c /Qopenmp ${_muscle_fp_library}
        $<$<CONFIG:Debug>:/Z7>
        ${_muscle_diagnostics}
        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>)
    set(MUSCLE_LINK_LIBRARY "")
    set(MUSCLE_OPT_TESTS
        /O3 /warn:all
        ${_muscle_fp_tests}
        $<$<CONFIG:Debug>:/Z7>
        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>)
    set(MUSCLE_LINK_TESTS
        $<$<CONFIG:Debug>:/Z7>
        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>)
    set(MUSCLE_FLAG_NO_INLINE /Ob0)
    set(MUSCLE_FLAGS_HOST_OFFLOAD -Qopenmp-targets:spir64)
    # Traps only in Debug: optimized, the fast model of the library evaluates operations speculatively (D-042)
    set(MUSCLE_FLAG_FPE_TRAPS $<$<CONFIG:Debug>:/fpe:0>)
    set(MUSCLE_FLAGS_COMPILE_FAIL /nologo /syntax-only)
    set(MUSCLE_FLAG_INCLUDE /I)
    set(MUSCLE_COMPILE_FAIL_PASS "error #(8289|6053)")
    set(MUSCLE_COMPILE_FAIL_FAIL "error #7002|catastrophic error")
    set(MUSCLE_LINK_STACK /F16777216)
else()
    # ifx on Linux
    if(_muscle_gpu_build)
        set(_muscle_fp_library "")
        set(_muscle_fp_tests "")
        set(_muscle_diagnostics "")
    else()
        # As on Windows (D-042, D-043): fast model in the library (safe speculation in Debug), precise model in the
        # tests
        set(_muscle_fp_library
            -assume minus0
            $<$<CONFIG:Debug>:-fp-speculation=safe>)
        set(_muscle_fp_tests
            -assume minus0
            $<$<NOT:$<CONFIG:Profile>>:-fp-model=precise>)
        # one `-check` only: CMake removes the repeated options of a target, and `-check bounds -check pointers`
        # would become `-check bounds pointers`
        set(_muscle_diagnostics $<$<CONFIG:Debug>:-warn all -check bounds,pointers,shape -traceback>)
    endif()
    set(MUSCLE_OPT_LIBRARY
        -fpp -qopenmp -fPIC ${_muscle_fp_library}
        $<$<CONFIG:Debug>:-g -O0>
        ${_muscle_diagnostics}
        $<$<CONFIG:Profile>:-Ofast -ipo -fp-model=fast -DNDEBUG>)
    set(MUSCLE_LINK_LIBRARY "")
    set(MUSCLE_OPT_TESTS
        -O3 -warn all
        ${_muscle_fp_tests}
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Profile>:-Ofast -ipo -fp-model=fast -DNDEBUG>)
    set(MUSCLE_LINK_TESTS
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Profile>:-Ofast -ipo -fp-model=fast -DNDEBUG>)
    set(MUSCLE_FLAG_NO_INLINE -fno-inline)
    set(MUSCLE_FLAGS_HOST_OFFLOAD -fopenmp-targets=spir64)
    set(MUSCLE_FLAG_FPE_TRAPS $<$<CONFIG:Debug>:-fpe0>)
    set(MUSCLE_FLAGS_COMPILE_FAIL -syntax-only)
    set(MUSCLE_FLAG_INCLUDE -I)
    set(MUSCLE_COMPILE_FAIL_PASS "error #(8289|6053)")
    set(MUSCLE_COMPILE_FAIL_FAIL "error #7002|catastrophic error")
endif()

function(muscle_fpe_traps)
    # Aborts a division by zero, an overflow or an invalid operation in the library, for the tests (targets) that
    # ask for it: the traps are set by the options of the main program
    foreach(_target IN LISTS ARGN)
        target_compile_options(${_target} PRIVATE ${MUSCLE_FLAG_FPE_TRAPS})
    endforeach()
endfunction()


function(muscle_add_compile_fail_test NAME SOURCE)
    # Compile-fail check (D-017.10b): passes when the compiler rejects SOURCE because it uses a PRIVATE component
    # of a kernel type; any other error (a missing module, a fatal error) fails it. The messages are those of each
    # compiler (the file name itself contains "private", so a loose pattern would always pass: ifx is matched by
    # its error numbers, #8289 and #6053).
    add_test(NAME ${NAME}
             COMMAND ${CMAKE_Fortran_COMPILER} ${MUSCLE_FLAGS_COMPILE_FAIL}
                     ${MUSCLE_FLAG_INCLUDE}$<TARGET_PROPERTY:muscle,Fortran_MODULE_DIRECTORY>
                     ${SOURCE})
    set_tests_properties(${NAME} PROPERTIES
                         PASS_REGULAR_EXPRESSION "${MUSCLE_COMPILE_FAIL_PASS}"
                         FAIL_REGULAR_EXPRESSION "${MUSCLE_COMPILE_FAIL_FAIL}")
endfunction()
