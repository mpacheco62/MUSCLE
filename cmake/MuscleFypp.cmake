include_guard(GLOBAL)
find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
set(MUSCLE_FYPP_SCRIPT "${PROJECT_SOURCE_DIR}/tools/fypp/fypp.py")
file(READ "${PROJECT_SOURCE_DIR}/tools/fypp/release.lock.json" _fypp_lock)
string(JSON _source_hash GET "${_fypp_lock}" source_sha256)
string(JSON _license_hash GET "${_fypp_lock}" license_sha256)
file(SHA256 "${MUSCLE_FYPP_SCRIPT}" _actual_source)
file(SHA256 "${PROJECT_SOURCE_DIR}/tools/fypp/LICENSE" _actual_license)
if(NOT _actual_source STREQUAL _source_hash OR NOT _actual_license STREQUAL _license_hash)
  message(FATAL_ERROR "Vendored Fypp release does not match its locked SHA256 hashes")
endif()
execute_process(COMMAND "${Python3_EXECUTABLE}" "${MUSCLE_FYPP_SCRIPT}" --version
  OUTPUT_VARIABLE _fypp_version RESULT_VARIABLE _fypp_status)
if(NOT _fypp_status EQUAL 0 OR NOT _fypp_version MATCHES "fypp 3[.]3[\r\n]*$")
  message(FATAL_ERROR "Fypp release 3.3 is required: ${_fypp_version}")
endif()
message(STATUS "MUSCLE: generating library sources with pinned Fypp 3.3")

function(muscle_generate_fypp TARGET TEMPLATE NAME)
  get_filename_component(_template "${TEMPLATE}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  set(_directory "${CMAKE_CURRENT_BINARY_DIR}/generated/${TARGET}")
  file(MAKE_DIRECTORY "${_directory}")
  set(_output "${_directory}/${NAME}.f90")
  add_custom_command(OUTPUT "${_output}"
    COMMAND "${Python3_EXECUTABLE}" "${MUSCLE_FYPP_SCRIPT}"
      -I "${PROJECT_SOURCE_DIR}/src" ${ARGN} --depfile "${_output}.d" "${_template}" "${_output}"
    DEPENDS "${_template}" "${MUSCLE_FYPP_SCRIPT}" "${PROJECT_SOURCE_DIR}/tools/fypp/release.lock.json"
    DEPFILE "${_output}.d" VERBATIM)
  # Some sources are registered from subdirectories of the target definition.
  # An explicit dependency connects their generation before Fortran scanning.
  add_custom_target(${TARGET}_${NAME}_fypp DEPENDS "${_output}")
  add_dependencies(${TARGET} ${TARGET}_${NAME}_fypp)
  target_sources(${TARGET} PRIVATE "${_output}")
endfunction()

# Existing public modules are emitted by Fypp from their canonical source.
# This preserves CPU module/type names and leaves other library sources alone.
function(muscle_add_source SOURCE)
  if(MUSCLE_USE_FYPP)
    get_filename_component(_name "${SOURCE}" NAME_WE)
    set(_wrapper "${CMAKE_CURRENT_BINARY_DIR}/${_name}.fypp")
    configure_file("${PROJECT_SOURCE_DIR}/cmake/cpu_source.fypp.in" "${_wrapper}" @ONLY)
    muscle_generate_fypp(muscle "${_wrapper}" "${_name}")
    # Fortran INCLUDEs inside the canonical module are resolved here.
    target_include_directories(muscle PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  else()
    target_sources(muscle PRIVATE "${SOURCE}")
  endif()
endfunction()
