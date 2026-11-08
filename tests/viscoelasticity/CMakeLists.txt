# Steps of a material point for the physical checks below (visco_test_steps)
muscle_fypp_preprocess(VISCO_TEST_STEPS_SRC visco_test_steps.fypp)
add_library(muscle_visco_test_steps STATIC ${VISCO_TEST_STEPS_SRC})
target_link_libraries(muscle_visco_test_steps PUBLIC muscle PRIVATE compile_options)
set_target_properties(muscle_visco_test_steps PROPERTIES
                      Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/tests)
target_include_directories(muscle_visco_test_steps PUBLIC ${CMAKE_BINARY_DIR}/modules/tests)

# Physical checks of Maxwell and Kelvin-Voigt through the kernel API (analytic relaxation, constant
# rate, time step limits, small-x series of gamma, initial elastic strain, frozen tangents in 2D)
foreach(_muscle_visco maxwell kelvin)
    add_muscle_fypp_test(${_muscle_visco} test_muscle_visco_${_muscle_visco}.fypp)
    target_link_libraries(test_${_muscle_visco} PRIVATE muscle_visco_test_steps)
endforeach()

# Golden (characterization) test of the viscoelasticity family (Maxwell with 0, 1, 2, 3 and
# MAX_VISCO_BRANCHES branches, Kelvin-Voigt, the frozen elasticity in 3D and 2D), multi-step load
# histories through the kernel API, double precision, within the tolerance of its category (D-021) against
# tests/reference/viscoelasticity/ (regenerate.sh there), captured from the original OO code.
add_muscle_fypp_test(visco_golden test_visco_golden.fypp)
set_tests_properties(viscoelasticity-visco_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the viscoelastic family (D-017.13) in every kind: load histories of Maxwell and
# Kelvin-Voigt through the generics of muscle_visco_laws, sp against dp, the status of invalid time
# steps and parameters (D-017.13b), the frozen elasticity 2D vs 3D, and a pure @:DEVICE() procedure
# that advances several steps at a point
add_muscle_fypp_test(visco_kinds test_visco_kinds.fypp)
