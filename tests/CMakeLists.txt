find_package(OpenMP REQUIRED COMPONENTS Fortran)
add_library(compile_options INTERFACE)
target_link_libraries(compile_options INTERFACE OpenMP::OpenMP_Fortran)
add_custom_target(build_all_tests)
# Test-only macros (tests/common/muscle_fd.fypp: central finite differences without closures): the templates
# of the tests include them by name and every test is regenerated when they change.
list(APPEND MUSCLE_FYPP_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/common")
list(APPEND MUSCLE_FYPP_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/common/muscle_fd.fypp")
# The tests that run kernels in an OpenMP target region (muscle_target_test_options, ctest label `target`)
add_custom_target(build_target_tests)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        message(STATUS "INTEL Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        if (WIN32)
            target_compile_options(compile_options INTERFACE "/O3" "/warn:all")
            target_compile_options(compile_options INTERFACE
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )
            target_link_options(compile_options INTERFACE
                                        $<$<CONFIG:Debug>:/Z7>
                                        $<$<CONFIG:Profile>:/Ofast /qipo /fp:fast /DNDEBUG>
                                )
        else()
            # Linux: the Windows-style options above are read as file names
            target_compile_options(compile_options INTERFACE "-O3" "-warn" "all")
            target_compile_options(compile_options INTERFACE
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile>:-Ofast -ipo -fp-model=fast -DNDEBUG>
                                )
            target_link_options(compile_options INTERFACE
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile>:-Ofast -ipo -fp-model=fast -DNDEBUG>
                                )
        endif()

        # target_compile_options(muscle PRIVATE 
        #                                 -c -Qopenmp -fpic
        #                         )
elseif  (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        message(STATUS "GNU Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
        # Unused dummy arguments are part of uniform interfaces; unused local variables still warn
        target_compile_options(compile_options INTERFACE "-O3" "-Wall" "-Wno-unused-dummy-argument")
        target_compile_options(compile_options INTERFACE
                                        $<$<CONFIG:Debug>:-g>
                                        $<$<CONFIG:Profile_unittests>:-g -O0 -fprofile-arcs -ftest-coverage>
                                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                                )
        target_link_options(compile_options INTERFACE 
                        $<$<CONFIG:Debug>:-g>
                        $<$<CONFIG:Profile_unittests>:--coverage>
                        $<$<CONFIG:Profile>:-Ofast -pg -flto -ffast-math -faggressive-loop-optimizations -DNDEBUG>
                        )
else()
    message(FATAL_ERROR "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}.")
endif()


# Shared assertions (muscle_test_utils: check_approx, check_true, check_scaled) and helpers
# (muscle_test_helpers: random generator, reference matrix operations, record) of the tests
muscle_fypp_preprocess(MUSCLE_TEST_UTILS_SRC common/muscle_test_utils.fypp common/muscle_test_helpers.fypp)
# muscle_test_reference: reference checks (tolerance by category, D-021) against tests/reference/ (register the
# test with muscle_reference_test); muscle_test_solver_sequence: strain
# sequence and record layout of the reference tests of the return mapping solvers
add_library(muscle_test_utils STATIC ${MUSCLE_TEST_UTILS_SRC} common/muscle_test_reference.f90
            materials/muscle_test_solver_sequence.f90)
target_link_libraries(muscle_test_utils PUBLIC muscle PRIVATE compile_options)
set_target_properties(muscle_test_utils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/tests)
target_include_directories(muscle_test_utils PUBLIC ${CMAKE_BINARY_DIR}/modules/tests)

# Store the root tests directory path
set(BASE_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Base tests directory")

function(add_muscle_test NAME)
    # Define executable name with 'test_' prefix
    set(EXE_NAME "test_${NAME}")

    # Create executable target
    add_executable(${EXE_NAME} ${ARGN})
    target_link_libraries(${EXE_NAME} PRIVATE muscle muscle_test_utils compile_options)
    target_include_directories(${EXE_NAME} PUBLIC ${CMAKE_BINARY_DIR}/modules/muscle)
    add_dependencies(build_all_tests ${EXE_NAME})


    # Get relative path from tests root to current folder
    file(RELATIVE_PATH REL_PATH "${BASE_TESTS_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    
    # Set delimiter matching VS Code 'cmake.ctest.testSuiteDelimiter'
    set(DELIMITER "-") 
    
    if (REL_PATH STREQUAL "")
        # No prefix if defined in the root tests folder
        set(TEST_DISPLAY_NAME "${NAME}")
    else()
        # Replace '/' with the delimiter for nested folders
        string(REPLACE "/" "${DELIMITER}" CATEGORY_PATH "${REL_PATH}")
        set(TEST_DISPLAY_NAME "${CATEGORY_PATH}${DELIMITER}${NAME}")
    endif()
    
    # Add test using prefix-free name for CTest and prefixed name for executable
    add_test(NAME ${TEST_DISPLAY_NAME} COMMAND ${EXE_NAME})
    set_target_properties(${EXE_NAME} PROPERTIES MUSCLE_TEST_NAME ${TEST_DISPLAY_NAME})
    target_link_libraries(${EXE_NAME} PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endfunction()


function(add_muscle_fypp_test NAME SOURCE)
    # Same as add_muscle_test, for a single .fypp template (see cmake/MuscleFypp.cmake).
    # DEPENDS lists the files the template includes (e.g. tests/common/*.fypp macros) so that
    # it is regenerated when they change.
    cmake_parse_arguments(PARSE_ARGV 2 arg "" "" "DEPENDS")
    muscle_fypp_preprocess(GENERATED ${SOURCE} DEPENDS ${arg_DEPENDS})
    add_muscle_test(${NAME} ${GENERATED})
    get_target_property(test_name test_${NAME} MUSCLE_TEST_NAME)
    set_tests_properties(${test_name} PROPERTIES LABELS kinds)
endfunction()


# Reference test (plan R): compares results with the records of tests/reference/ (muscle_test_reference);
# sets the directory of the references and the label `reference`. Extra ENVIRONMENT entries (e.g.
# MUSCLE_REFERENCE_STRICT=1) are passed after ENV.
function(muscle_reference_test TEST_NAME)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "" "ENV")
    set(env "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference" ${arg_ENV})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS reference ENVIRONMENT "${env}")
endfunction()


function(muscle_target_test_options TARGET)
    # Options of a test that runs kernels in an OpenMP target region and compares them with the CPU
    # (testing.md, GPU). Without MUSCLE_GPU the library has no `declare target`: an offload-enabled
    # gfortran (e.g. Ubuntu's, offload-defaulted to nvptx/amdgcn) would fail to link the device image
    # (unresolved library symbols), so the region is compiled for host fallback only. -fno-inline keeps
    # one copy of each point procedure for both calls (inlined copies may be vectorized differently,
    # and the host comparison is bit for bit).
    # The test gets the ctest label `target` and the target joins build_target_tests (tools/check.sh
    # builds and runs those in the GPU builds).
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND NOT MUSCLE_GPU_ENABLED)
        target_compile_options(${TARGET} PRIVATE -foffload=disable -fno-inline)
        target_link_options(${TARGET} PRIVATE -foffload=disable)
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        target_compile_options(${TARGET} PRIVATE -fno-inline)
    endif()
    get_target_property(test_name ${TARGET} MUSCLE_TEST_NAME)
    set_tests_properties(${test_name} PROPERTIES LABELS target)
    add_dependencies(build_target_tests ${TARGET})
endfunction()


function(add_muscle_fypp_target_tests NAME SOURCE GROUP)
    # One target test per case of GROUP of tests/gpu/target_cases.py (`tensor`, `family`, `material`, `batch_dev`):
    # SOURCE preprocessed with the fypp define TARGET_CASE='<case>' into its own executable
    # test_<NAME>_<case> (ctest <dir>-<NAME>_<case>), so that every device image stays small. Each one
    # keeps its module files in a directory of its own (the templates define the same modules).
    cmake_parse_arguments(PARSE_ARGV 3 arg "" "" "DEPENDS")
    set(cases_py "${BASE_TESTS_DIR}/gpu/target_cases.py")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=${MUSCLE_FYPP_META_DIR}
                ${Python3_EXECUTABLE} "${cases_py}" ${GROUP} "${MUSCLE_REAL_KINDS}"
        OUTPUT_VARIABLE cases
        RESULT_VARIABLE cases_result
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT cases_result EQUAL 0)
        message(FATAL_ERROR "MUSCLE: ${cases_py} ${GROUP} failed")
    endif()
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${cases_py}")
    _muscle_fypp_source("${SOURCE}" source_abs source_rel)
    string(REGEX REPLACE "\\.fypp$" "" output_base "${source_rel}")
    foreach(case IN LISTS cases)
        set(output "${MUSCLE_BINARY_DIR}/fypp/${output_base}_${case}.f90")
        _muscle_fypp_command("${source_abs}" "${output}" "FYPP ${source_rel} (${case})"
                             FLAGS "-DTARGET_CASE='${case}'" DEPENDS "${cases_py}" ${arg_DEPENDS})
        set_source_files_properties("${output}" PROPERTIES Fortran_PREPROCESS ON)
        add_muscle_test(${NAME}_${case} "${output}")
        set_target_properties(test_${NAME}_${case} PROPERTIES
                              Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules_${NAME}_${case}")
        muscle_target_test_options(test_${NAME}_${case})
    endforeach()
endfunction()


add_subdirectory(common)
add_subdirectory(tensors)
add_subdirectory(kinematics)
add_subdirectory(math)
add_subdirectory(elasticity)
add_subdirectory(viscoelasticity)
add_subdirectory(plasticity)
add_subdirectory(materials)
add_subdirectory(gpu)
add_subdirectory(tools)


