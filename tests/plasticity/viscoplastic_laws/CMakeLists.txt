# The Chaboche law (independent reference value, chain rule, equivalence with Perzyna, threshold)
add_muscle_fypp_test(vp_chaboche test_muscle_vp_chaboche.fypp)

# Kernel API of the viscoplastic and overstress laws (D-017) in every kind: free functions through
# the generics of muscle_vp_laws/muscle_vp_overstress_laws (composed laws with the hardening values
# evaluated outside, D-017.6b), default derivatives of the contract, and a pure @:DEVICE() point
# procedure that composes hard_* with vp_* and os_*
add_muscle_fypp_test(vp_laws_kinds test_muscle_vp_laws_kinds.fypp)

# Limits, clamps and branches of the laws in every kind, at the documented boundaries (value on the
# boundary and on both sides), against closed forms: thermal factor of the full Johnson-Cook, rate
# floors of Johnson-Cook, Perzyna at psi = 0, VA with m = 0 and 1, the clamps of the derivatives of
# Rusinek-Klepaczko and MRK, and the clamp of X in NNL
add_muscle_fypp_test(vp_limits test_muscle_vp_laws_limits.fypp)
