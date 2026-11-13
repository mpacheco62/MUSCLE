# Kernel API of the kinematics family (D-017.14) in every kind: sp against dp, the status of invalid
# deformation gradients (D-017.14a), objectivity, and pure @:DEVICE() point procedures in an OpenMP
# target region (muscle_target_test_options, tests/CMakeLists.txt)
add_muscle_fypp_test(kin_kinds test_muscle_kin_laws_kinds.fypp)
muscle_target_test_options(test_kin_kinds)

# Unit tests of muscle_kin_finite and of the corotational update in every kind against an independent truth
# in quad precision (plan R): strain measures, J and R of F = R0 Q diag(s) Q^T, the ten stress conversions,
# the geometric terms of the spatial tangent of every kinematics, the material parts of the logarithmic and
# Euler-Almansi spatial tangents (Daleckii-Krein; the truth checked against finite differences), on distinct,
# nearly repeated (log series) and dropped eigenvalues of C, and Hencky strains near the identity (D-032). It replaces the former reference test
# kin_golden (plan R, tests/mutations/kin.toml).
add_muscle_fypp_test(kin_finite test_muscle_kin_finite.fypp)
