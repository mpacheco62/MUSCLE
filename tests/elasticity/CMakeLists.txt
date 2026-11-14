# Kernel API of the elastic family (D-017.12) in every kind and dimension: the linear model through
# the generics of muscle_elas_laws against Hooke law (dp reference), 2D vs 3D embedding of
# every function, the generics against the dp specific functions (bit for bit), and pure @:DEVICE()
# point procedures
add_muscle_fypp_test(elas_kinds test_muscle_elas_laws_kinds.fypp)

# Compile-fail check (D-017.10b): every component of the kernel type elas_linear_<k> is derived
# (private), so it can only be built with elas_linear_init (muscle_add_compile_fail_test,
# cmake/MuscleCompiler.cmake)
muscle_add_compile_fail_test(elasticity-linear_private_params
                             ${CMAKE_CURRENT_SOURCE_DIR}/linear_private_params.f90)
