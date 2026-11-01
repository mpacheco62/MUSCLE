add_muscle_test(maxwell test_muscle_visco_maxwell.f90)
add_muscle_test(kelvin test_muscle_visco_kelvin.f90)

# Golden (characterization) test of the viscoelasticity family (Generalized_maxwell with 0, 1, 2, 3
# and MAX_VISCO_BRANCHES branches, Kelvin_voigt, Elasticity_frozen through Base_elasticity in 3D and
# 2D, Viscoelastic_history), multi-step load histories through the polymorphic base, double
# precision, bit for bit against tests/reference/viscoelasticity/ (regenerate.sh there), captured
# from the original OO code.
add_muscle_test(visco_golden test_visco_golden.f90)
set_tests_properties(viscoelasticity-visco_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the viscoelastic family (D-017.13) in every kind: load histories of Maxwell and
# Kelvin-Voigt through the generics of muscle_visco_laws against the historical OO types (dp, bit for
# bit) and sp against dp, the status of invalid time steps and parameters (D-017.13b), the frozen
# elasticity 2D vs 3D, and a pure @:DEVICE() procedure that advances several steps at a point
add_muscle_fypp_test(visco_kinds test_visco_kinds.fypp)
