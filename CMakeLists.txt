cmake_minimum_required(VERSION 3.20)
# set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CXX_FLAGS_RELEASE};profiling" CACHE STRING "Available build types.")
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;Profile;Profile_unittests" CACHE STRING "Available build types.")

project(MUSCLE VERSION 0.1
        DESCRIPTION "Mechanics User-defined Solids Constitutive Library & Extensions"
        LANGUAGES Fortran C) # C for the tests
include(CTest)
enable_testing()

enable_language(Fortran)

message(STATUS "********MUSCLE - CMAKE_BUILD_TYPE is: ${CMAKE_BUILD_TYPE}********")

# This generate a file compile_commands.json with commands use on compile
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set_property(CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE "${CMAKE_SOURCE_DIR}/output")
endif()

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(MuscleFypp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Source code
add_subdirectory(src)
add_subdirectory(tests)

# Benchmarks (optional, not part of the test suite)
option(MUSCLE_BUILD_BENCHMARKS "Build the performance benchmarks in benchmarks/" OFF)
if(MUSCLE_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# Documentation (FORD, optional): `cmake --build <dir> --target doc` writes <dir>/doc. FORD runs
# from the source root because the fypp module path (-M tools/fypp) in my_project.md is relative.
find_program(FORD_EXECUTABLE ford)
if(FORD_EXECUTABLE)
    add_custom_target(doc
        COMMAND ${CMAKE_COMMAND} -E env PYTHONDONTWRITEBYTECODE=1
                ${FORD_EXECUTABLE} --output_dir "${PROJECT_BINARY_DIR}/doc" my_project.md
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        COMMENT "FORD documentation in ${PROJECT_BINARY_DIR}/doc"
        VERBATIM)
endif()


# Solo configurar el reporte si estamos en el modo de cobertura
if(CMAKE_BUILD_TYPE STREQUAL "Profile_unittests")

    find_program(LCOV_PATH lcov)
    find_program(GENHTML_PATH genhtml)

    if(LCOV_PATH AND GENHTML_PATH)
        message(STATUS "** Se encontro lcov y genhtml: configurando target de cobertura **")
        # lcov >= 2 reports the line markers of the templates (#line, D-009) as inconsistent function
        # end lines and stops: tolerate that category (lcov 1.x does not know it)
        execute_process(COMMAND ${LCOV_PATH} --version OUTPUT_VARIABLE _lcov_version_out ERROR_QUIET)
        string(REGEX MATCH "[0-9]+\\.[0-9]+" _lcov_version "${_lcov_version_out}")
        set(LCOV_IGNORE)
        if(_lcov_version VERSION_GREATER_EQUAL 2.0)
            set(LCOV_IGNORE --ignore-errors inconsistent,inconsistent,unused,unused)
        endif()
        
        add_custom_target(coverage_report
            # # PASO 0: Limpiar contadores antiguos para evitar datos basura
            COMMAND ${CMAKE_COMMAND} -E echo " *** PASO 0: Limpiando contadores previos..."
            COMMAND ${LCOV_PATH} --directory . --zerocounters

            # Re-ejecutar tests (aunque fallen, continuará con el siguiente PASO)
            COMMAND ${CMAKE_COMMAND} -E echo " *** Ejecutando CTest (se ignorarán fallos para generar el reporte)..."
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure || ${CMAKE_COMMAND} -E echo "Algún test falló, pero continuamos..."

            # PASO 1: Capturar datos brutos
            COMMAND ${CMAKE_COMMAND} -E echo " *** PASO 1: Capturando datos de cobertura (.gcda)..."
            COMMAND ${LCOV_PATH} --directory . --capture --output-file coverage.info --rc lcov_branch_coverage=1 ${LCOV_IGNORE}
            
            # PASO 2: Filtrar para dejar solo la librería (src/ e include/; sin los tests)
            COMMAND ${CMAKE_COMMAND} -E echo " *** PASO 2: Filtrando archivos externos y tests..."
            COMMAND ${LCOV_PATH} --extract coverage.info "${CMAKE_SOURCE_DIR}/src/*" "${CMAKE_SOURCE_DIR}/include/*"
                    --output-file coverage_cleaned.info --rc lcov_branch_coverage=1 ${LCOV_IGNORE}
            
            # PASO 3: Generar el reporte HTML
            COMMAND ${CMAKE_COMMAND} -E echo " *** PASO 3: Generando reporte HTML final..."
            COMMAND ${GENHTML_PATH} coverage_cleaned.info --output-directory html_report --legend ${LCOV_IGNORE}
            
            COMMAND ${CMAKE_COMMAND} -E echo " *** EXITO: Reporte generado en ${CMAKE_BINARY_DIR}/html_report/index.html"
            
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Procesando cobertura de codigo..."
            VERBATIM
        )
        add_dependencies(coverage_report muscle)
        add_dependencies(coverage_report build_all_tests)
    else()
        message(WARNING "lcov o genhtml no encontrados. No se podra generar el reporte visual.")
    endif()

endif()