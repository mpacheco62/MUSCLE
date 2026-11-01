add_muscle_test(bilinear test_muscle_hard_bilinear.f90)
add_muscle_test(ludwik test_muscle_hard_ludwik.f90)
add_muscle_test(swift test_muscle_hard_swift.f90)
add_muscle_test(vocemod test_muscle_hard_vocemod.f90)

# Golden (characterization) test of every hardening law (stress, dstress_dep, ddstress_ddep
# through the polymorphic base, constructors), double precision, bit for bit against
# tests/reference/plasticity/hardening_laws/ (regenerate.sh there), captured from the original OO code.
add_muscle_test(hard_golden test_hard_golden.f90)
set_tests_properties(plasticity-hardening_laws-hard_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the hardening laws (D-017) in every kind: free functions of the hardening_model
# contract macro through the generics of muscle_hard_laws, finite-difference defaults of the macro,
# historical OO type against the dp functions, and a pure @:DEVICE() point procedure
add_muscle_fypp_test(hard_laws_kinds test_hard_laws_kinds.fypp)
