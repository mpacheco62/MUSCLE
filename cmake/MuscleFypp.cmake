# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# FYPP integration for MUSCLE.
#
# Templates (*.fypp) are preprocessed into ${MUSCLE_BINARY_DIR}/fypp/<same relative path>.f90.
# The generated files are build artifacts and are never committed.
#
#   muscle_fypp_preprocess(<out_var> <file.fypp>... [DEPENDS <file>...])
#       Adds the generation rules and returns the generated .f90 paths in <out_var>.
#   muscle_add_fypp_sources(<target> <file.fypp>... [DEPENDS <file>...])
#       Generates the files and adds them to <target>, which may be defined in another directory.
#   muscle_add_materials(<target> <template.fypp> [DEPENDS <file>...])
#       Generates one module per declared material (D-008, D-018) from the same template:
#       ${MUSCLE_BINARY_DIR}/fypp/<template dir>/muscle_mat_<name>.f90, with -DMUSCLE_MATERIAL='<name>'
#       (and -DMUSCLE_MATERIALS_FILE='<path>' for a user registry). The registry is
#       tools/fypp/muscle_materials.py or the file given by -DMUSCLE_MATERIALS=<path>; muscle_meta
#       validates it at configure time. The bodies of the solvers (muscle_meta.MATERIAL_SOLVERS,
#       `muscle_meta.py solver-bodies`) are dependencies; DEPENDS: the other files the template includes.
#
# DEPENDS lists extra files the templates include (e.g. a body shared by an area); the files
# in MUSCLE_FYPP_DEPENDS are always dependencies. Build-time switches that a template needs
# (USE_MKL, ENABLE_UDTIO, ...) are appended to the MUSCLE_FYPP_DEFINES variable before the
# call and reach fypp as -D<name>; templates test them with `#:if defined('<name>')`.

include_guard(GLOBAL)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

find_program(FYPP_EXECUTABLE fypp)
if(NOT FYPP_EXECUTABLE)
    message(FATAL_ERROR
        "MUSCLE: fypp not found. Install it with `pip install fypp` "
        "or `conda env update -f environment.yml`.")
endif()

execute_process(COMMAND ${FYPP_EXECUTABLE} --version
                OUTPUT_VARIABLE _muscle_fypp_version_out
                OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REGEX MATCH "[0-9]+\\.[0-9]+(\\.[0-9]+)?" FYPP_VERSION "${_muscle_fypp_version_out}")
if(FYPP_VERSION VERSION_LESS 3.1)
    message(FATAL_ERROR "MUSCLE: fypp >= 3.1 is required, found '${_muscle_fypp_version_out}'.")
endif()
message(STATUS "MUSCLE: fypp ${FYPP_VERSION} (${FYPP_EXECUTABLE})")

# The allowed kinds are validated again inside the templates (muscle_common.fypp,
# muscle_meta.KINDS); this check only gives an early, readable error at configure time.
set(MUSCLE_REAL_KINDS "sp;dp" CACHE STRING
    "Real kinds instantiated by the FYPP templates (dp is always required)")
list(REMOVE_DUPLICATES MUSCLE_REAL_KINDS)
foreach(_kind IN LISTS MUSCLE_REAL_KINDS)
    if(NOT _kind MATCHES "^(sp|dp)$")
        message(FATAL_ERROR "MUSCLE: unknown kind '${_kind}' in MUSCLE_REAL_KINDS (allowed: sp, dp).")
    endif()
endforeach()
# The historical type names are aliases of the dp instances (D-003), so dp cannot be dropped.
if(NOT "dp" IN_LIST MUSCLE_REAL_KINDS)
    message(FATAL_ERROR "MUSCLE: MUSCLE_REAL_KINDS must contain dp.")
endif()
message(STATUS "MUSCLE: real kinds generated: ${MUSCLE_REAL_KINDS}")

set(MUSCLE_FYPP_INCLUDE_DIR "${MUSCLE_SOURCE_DIR}/include")
set(MUSCLE_FYPP_META_DIR "${MUSCLE_SOURCE_DIR}/tools/fypp")
# Every template is regenerated when one of these files changes.
set(MUSCLE_FYPP_DEPENDS
    "${MUSCLE_FYPP_INCLUDE_DIR}/muscle_common.fypp"
    "${MUSCLE_FYPP_INCLUDE_DIR}/muscle_tensor_common.fypp"
    "${MUSCLE_FYPP_META_DIR}/muscle_meta.py")

# Python list literal for fypp, e.g. "sp;dp" -> "['sp','dp']"
string(REPLACE ";" "','" _muscle_kinds_py "${MUSCLE_REAL_KINDS}")
# Line markers (all configurations) make compiler errors and coverage point to the .fypp
# template, and keep the command identical across configurations so multi-config generators
# share one generated file. Only the `#line` (std) format survives CMake's separate Fortran
# preprocessing step without "linemarker ignored due to incorrect nesting" warnings (D-009).
# No markers inside a continued statement (nocontlines): fypp folds long lines and a `#line`
# between the two halves of a folded character literal is an "unterminated character constant"
# for ifx (gfortran accepts it).
set(MUSCLE_FYPP_FLAGS
    "-I${MUSCLE_FYPP_INCLUDE_DIR}"
    "-M${MUSCLE_FYPP_META_DIR}" "-mmuscle_meta"
    "-DMUSCLE_REAL_KINDS=['${_muscle_kinds_py}']"
    "--line-numbering" "--line-marker-format=std" "--line-numbering-mode=nocontlines")


# Material registry (D-018.9): empty = the test registry tools/fypp/muscle_materials.py
set(MUSCLE_MATERIALS "" CACHE FILEPATH
    "Python file with the declared materials (MATERIALS); empty: the test registry tools/fypp/muscle_materials.py")


# _muscle_fypp_source(<source> <source_abs_var> <source_rel_var>): absolute and repository-relative
# path of a template, checked (a .fypp inside the source tree).
function(_muscle_fypp_source SOURCE ABS_VAR REL_VAR)
    get_filename_component(source_abs "${SOURCE}" ABSOLUTE)
    if(NOT source_abs MATCHES "\\.fypp$")
        message(FATAL_ERROR "MUSCLE: '${SOURCE}' is not a .fypp template.")
    endif()
    file(RELATIVE_PATH source_rel "${MUSCLE_SOURCE_DIR}" "${source_abs}")
    if(source_rel MATCHES "^\\.\\.")
        message(FATAL_ERROR "MUSCLE: template '${source_abs}' is outside ${MUSCLE_SOURCE_DIR}.")
    endif()
    set(${ABS_VAR} "${source_abs}" PARENT_SCOPE)
    set(${REL_VAR} "${source_rel}" PARENT_SCOPE)
endfunction()


# _muscle_fypp_command(<source_abs> <output> <comment> FLAGS <flag>... DEPENDS <file>...): rule that
# preprocesses one template into <output> with the common flags, the build switches and FLAGS.
function(_muscle_fypp_command SOURCE_ABS OUTPUT COMMENT)
    cmake_parse_arguments(PARSE_ARGV 3 arg "" "" "FLAGS;DEPENDS")
    set(defines)
    foreach(define IN LISTS MUSCLE_FYPP_DEFINES)
        list(APPEND defines "-D${define}")
    endforeach()
    get_filename_component(output_dir "${OUTPUT}" DIRECTORY)
    add_custom_command(
        OUTPUT "${OUTPUT}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${output_dir}"
        # No bytecode: fypp imports muscle_meta from the source tree
        COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1
                ${FYPP_EXECUTABLE} ${MUSCLE_FYPP_FLAGS} ${defines} ${arg_FLAGS} "${SOURCE_ABS}" "${OUTPUT}"
        DEPENDS "${SOURCE_ABS}" ${MUSCLE_FYPP_DEPENDS} ${arg_DEPENDS}
        COMMENT "${COMMENT}"
        VERBATIM)
endfunction()


function(muscle_fypp_preprocess OUT_VAR)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "" "DEPENDS")
    set(generated)
    foreach(source IN LISTS arg_UNPARSED_ARGUMENTS)
        _muscle_fypp_source("${source}" source_abs source_rel)
        string(REGEX REPLACE "\\.fypp$" ".f90" output_rel "${source_rel}")
        set(output "${MUSCLE_BINARY_DIR}/fypp/${output_rel}")
        _muscle_fypp_command("${source_abs}" "${output}" "FYPP ${source_rel}" DEPENDS ${arg_DEPENDS})
        list(APPEND generated "${output}")
    endforeach()
    # Needed by the `#line` markers (see MUSCLE_FYPP_FLAGS)
    set_source_files_properties(${generated} PROPERTIES Fortran_PREPROCESS ON)
    set(${OUT_VAR} ${generated} PARENT_SCOPE)
endfunction()


function(muscle_add_fypp_sources TARGET)
    muscle_fypp_preprocess(generated ${ARGN})
    _muscle_add_generated(${TARGET} ${generated})
endfunction()


# _muscle_add_generated(<target> <file.f90>...): adds generated sources to <target>, which may be
# defined in another directory.
function(_muscle_add_generated TARGET)
    set(generated ${ARGN})
    target_sources(${TARGET} PRIVATE ${generated})
    set_source_files_properties(${generated} TARGET_DIRECTORY ${TARGET}
                                PROPERTIES Fortran_PREPROCESS ON)

    # A custom command only gets a build rule in the directory where it is defined. When the
    # target lives elsewhere (e.g. `muscle` in src/), a local target drives the generation.
    get_target_property(target_dir ${TARGET} SOURCE_DIR)
    if(NOT target_dir STREQUAL CMAKE_CURRENT_SOURCE_DIR)
        file(RELATIVE_PATH current_rel "${MUSCLE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
        string(MAKE_C_IDENTIFIER "fypp_${TARGET}_${current_rel}" driver_base)
        set(driver "${driver_base}")
        set(index 1)
        while(TARGET ${driver})
            math(EXPR index "${index} + 1")
            set(driver "${driver_base}_${index}")
        endwhile()
        add_custom_target(${driver} DEPENDS ${generated})
        add_dependencies(${TARGET} ${driver})
    endif()
endfunction()


function(muscle_add_materials TARGET TEMPLATE)
    cmake_parse_arguments(PARSE_ARGV 2 arg "" "" "DEPENDS")
    _muscle_fypp_source("${TEMPLATE}" template_abs template_rel)
    get_filename_component(template_dir "${template_rel}" DIRECTORY)

    set(registry_flags)
    set(registry "${MUSCLE_FYPP_META_DIR}/muscle_materials.py")
    if(MUSCLE_MATERIALS)
        get_filename_component(registry "${MUSCLE_MATERIALS}" ABSOLUTE)
        set(registry_flags "-DMUSCLE_MATERIALS_FILE='${registry}'")
    endif()
    if(NOT EXISTS "${registry}")
        message(FATAL_ERROR "MUSCLE: material registry '${registry}' not found (MUSCLE_MATERIALS).")
    endif()
    # Names of the validated registry (muscle_meta.load_materials); reconfigure when it changes
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1
                ${Python3_EXECUTABLE} "${MUSCLE_FYPP_META_DIR}/muscle_meta.py" materials "${registry}"
        OUTPUT_VARIABLE materials
        ERROR_VARIABLE materials_error
        RESULT_VARIABLE materials_result
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT materials_result EQUAL 0)
        message(FATAL_ERROR "MUSCLE: invalid material registry '${registry}':\n${materials_error}")
    endif()
    # (the validation lives in muscle_meta.py: a change there may change the list too)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                 "${registry}" "${MUSCLE_FYPP_META_DIR}/muscle_meta.py")
    message(STATUS "MUSCLE: materials (${registry}): ${materials}")
    # For the tests generated for every material of the registry (tests/materials/test_mat_smoke.fypp)
    set_property(GLOBAL PROPERTY MUSCLE_MATERIALS_REGISTRY "${registry}")
    set_property(GLOBAL PROPERTY MUSCLE_MATERIALS_FLAGS "${registry_flags}")
    # Bodies of the solvers of muscle_meta.MATERIAL_SOLVERS (src/solvers/), included by the template
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1
                ${Python3_EXECUTABLE} "${MUSCLE_FYPP_META_DIR}/muscle_meta.py" solver-bodies
        OUTPUT_VARIABLE bodies
        RESULT_VARIABLE bodies_result
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT bodies_result EQUAL 0)
        message(FATAL_ERROR "MUSCLE: muscle_meta.py solver-bodies failed")
    endif()
    list(TRANSFORM bodies PREPEND "${MUSCLE_SOURCE_DIR}/src/solvers/")

    set(generated)
    foreach(name IN LISTS materials)
        set(output "${MUSCLE_BINARY_DIR}/fypp/${template_dir}/muscle_mat_${name}.f90")
        _muscle_fypp_command("${template_abs}" "${output}" "FYPP ${template_rel} (material ${name})"
                             FLAGS "-DMUSCLE_MATERIAL='${name}'" ${registry_flags}
                             DEPENDS "${registry}" ${bodies} ${arg_DEPENDS})
        list(APPEND generated "${output}")
    endforeach()
    _muscle_add_generated(${TARGET} ${generated})
endfunction()

