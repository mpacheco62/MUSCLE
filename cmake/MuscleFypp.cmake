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
#
# DEPENDS lists extra files the templates include (e.g. a body shared by an area); the files
# in MUSCLE_FYPP_DEPENDS are always dependencies. Build-time switches that a template needs
# (USE_MKL, ENABLE_UDTIO, ...) are appended to the MUSCLE_FYPP_DEFINES variable before the
# call and reach fypp as -D<name>; templates test them with `#:if defined('<name>')`.

include_guard(GLOBAL)

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
    "${MUSCLE_FYPP_META_DIR}/muscle_meta.py")

# Python list literal for fypp, e.g. "sp;dp" -> "['sp','dp']"
string(REPLACE ";" "','" _muscle_kinds_py "${MUSCLE_REAL_KINDS}")
# Line markers (all configurations) make compiler errors and coverage point to the .fypp
# template, and keep the command identical across configurations so multi-config generators
# share one generated file. Only the `#line` (std) format survives CMake's separate Fortran
# preprocessing step without "linemarker ignored due to incorrect nesting" warnings (D-009).
set(MUSCLE_FYPP_FLAGS
    "-I${MUSCLE_FYPP_INCLUDE_DIR}"
    "-M${MUSCLE_FYPP_META_DIR}" "-mmuscle_meta"
    "-DMUSCLE_REAL_KINDS=['${_muscle_kinds_py}']"
    "--line-numbering" "--line-marker-format=std")


function(muscle_fypp_preprocess OUT_VAR)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "" "DEPENDS")
    set(defines)
    foreach(define IN LISTS MUSCLE_FYPP_DEFINES)
        list(APPEND defines "-D${define}")
    endforeach()

    set(generated)
    foreach(source IN LISTS arg_UNPARSED_ARGUMENTS)
        get_filename_component(source_abs "${source}" ABSOLUTE)
        if(NOT source_abs MATCHES "\\.fypp$")
            message(FATAL_ERROR "MUSCLE: '${source}' is not a .fypp template.")
        endif()
        file(RELATIVE_PATH source_rel "${MUSCLE_SOURCE_DIR}" "${source_abs}")
        if(source_rel MATCHES "^\\.\\.")
            message(FATAL_ERROR "MUSCLE: template '${source_abs}' is outside ${MUSCLE_SOURCE_DIR}.")
        endif()
        string(REGEX REPLACE "\\.fypp$" ".f90" output_rel "${source_rel}")
        set(output "${MUSCLE_BINARY_DIR}/fypp/${output_rel}")
        get_filename_component(output_dir "${output}" DIRECTORY)

        add_custom_command(
            OUTPUT "${output}"
            COMMAND ${CMAKE_COMMAND} -E make_directory "${output_dir}"
            COMMAND ${FYPP_EXECUTABLE} ${MUSCLE_FYPP_FLAGS} ${defines} "${source_abs}" "${output}"
            DEPENDS "${source_abs}" ${MUSCLE_FYPP_DEPENDS} ${arg_DEPENDS}
            COMMENT "FYPP ${source_rel}"
            VERBATIM)
        list(APPEND generated "${output}")
    endforeach()
    # Needed by the `#line` markers (see MUSCLE_FYPP_FLAGS)
    set_source_files_properties(${generated} PROPERTIES Fortran_PREPROCESS ON)
    set(${OUT_VAR} ${generated} PARENT_SCOPE)
endfunction()


function(muscle_add_fypp_sources TARGET)
    muscle_fypp_preprocess(generated ${ARGN})
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
