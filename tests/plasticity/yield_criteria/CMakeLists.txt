# Kernel API of the yield criteria (D-017) in every kind: free functions through the generics of
# muscle_yield_laws against closed forms and central differences (dp) and against dp (sp), zero
# derivatives at a null equivalent stress (D-017.9), default joint evaluation, finite-difference
# finite-difference derivatives (FD_GRAD/FD_HESS of the stress), a pure @:DEVICE() point procedure, and
# CPB06 with |k| = 1 at pure shears (closed forms; Hessian across a null principal value). It replaced
# the reference test of the criteria (plan R, tests/mutations/yield.toml)
add_muscle_fypp_test(yield_kinds test_muscle_yield_laws_kinds.fypp)

# Compile-fail checks (D-017.10b): the kernel types with derived parameters (CPB06: B, L;
# Drucker-Prager: f_factor, I1_factor; CPB06ex2: L, Lp, B; Yld2004: Lp, Lpp) have them private, so
# they cannot be set outside the model module (only by yield_<m>_init). Each test passes when the
# compiler rejects <m>_private_params.f90 (muscle_add_compile_fail_test, cmake/MuscleCompiler.cmake)
foreach(_muscle_model cpb06 druckerprager cpb06ex2 yld2004)
    muscle_add_compile_fail_test(plasticity-yield_criteria-${_muscle_model}_private_params
                                 ${CMAKE_CURRENT_SOURCE_DIR}/${_muscle_model}_private_params.f90)
endforeach()
