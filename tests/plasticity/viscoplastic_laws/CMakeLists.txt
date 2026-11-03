# Physical checks of the laws through the kernel API: fits against the Python reference values, and
# the Chaboche law (independent reference value, chain rule, equivalence with Perzyna, threshold)
add_muscle_test(vp_models test_muscle_vp_models.f90)
add_muscle_fypp_test(vp_chaboche test_muscle_vp_chaboche.fypp)

# Golden (characterization) test of every viscoplastic law (flow_stress/dstress_dep and
# flow_rate/dflow_dsigma/dflow_dep through the kernel API, parameter types), composed with
# bilinear/Swift/modified Voce hardening, double precision, bit for bit against
# tests/reference/plasticity/viscoplastic_laws/ (regenerate.sh there), captured from the original
# OO code. Includes ep = 0 and epd = 0 (first plastic increment).
add_muscle_fypp_test(vp_golden test_vp_golden.fypp)
set_tests_properties(plasticity-viscoplastic_laws-vp_golden PROPERTIES
                     ENVIRONMENT "MUSCLE_TEST_REFERENCE_DIR=${BASE_TESTS_DIR}/reference")

# Kernel API of the viscoplastic and overstress laws (D-017) in every kind: free functions through
# the generics of muscle_vp_laws/muscle_vp_overstress_laws (composed laws with the hardening values
# evaluated outside, D-017.6b), default derivatives of the contract, and a pure @:DEVICE() point
# procedure that composes hard_* with vp_* and os_*
add_muscle_fypp_test(vp_laws_kinds test_vp_laws_kinds.fypp)
