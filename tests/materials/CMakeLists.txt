# Declared materials (Phase 7, D-018).
#
# Smoke test of every material of the registry the library was built with (the test registry or the
# user registry of -DMUSCLE_MATERIALS): init, elastic step, history round trip, batch == point. The
# template reads the registry passed to the material modules (MUSCLE_MATERIALS_FILE).
get_property(_registry GLOBAL PROPERTY MUSCLE_MATERIALS_REGISTRY)
get_property(_registry_flags GLOBAL PROPERTY MUSCLE_MATERIALS_FLAGS)
set(_defines "${MUSCLE_FYPP_DEFINES}")
if(_registry_flags)
    list(APPEND MUSCLE_FYPP_DEFINES "MUSCLE_MATERIALS_FILE='${_registry}'")
endif()
add_muscle_fypp_test(mat_smoke test_muscle_material_smoke.fypp
                     DEPENDS "${_registry}" ${CMAKE_CURRENT_SOURCE_DIR}/props_samples.py
                             ${CMAKE_CURRENT_SOURCE_DIR}/case_params.py)
set(MUSCLE_FYPP_DEFINES "${_defines}")

# Point kernels of the materials of the test registry (tools/fypp/muscle_materials.py). With a user
# registry (-DMUSCLE_MATERIALS) these modules do not exist and the tests are not built.
if(NOT MUSCLE_MATERIALS)
    # Strain sequence tests of the return mapping solvers, one per solver (the cases of the historical
    # reference tests): sp against dp, analytic vs FD tangent, coverage of the registry; the common part is
    # mat_sequence_common.fypp, the parameter sets case_params.py.
    foreach(solver IN ITEMS closest_point closest_point_perzyna closest_point_time closest_point_2 cutting_plane)
        add_muscle_fypp_test(mat_${solver} test_muscle_material_${solver}.fypp
                             DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mat_sequence_common.fypp
                                     ${CMAKE_CURRENT_SOURCE_DIR}/case_params.py)
    endforeach()
    # Reduced reference (plan R, R6): the solvers whose mutations the unit tests cannot all detect keep the
    # dp reference records of a few cases (tests/reference/solvers/<solver>_dp.txt, within the tolerance of its
    # category, D-021; `--write` regenerates it, regenerate.sh there). Each `mat_<solver>_reference_mutation`
    # test runs it with `--mutate rel` (one parameter perturbed by 1e-6 relative, above every reference
    # tolerance) and passes only if the reference comparison reports it; `*_reference_mutation_strict` runs
    # `--mutate ulp` (one ulp) under MUSCLE_REFERENCE_STRICT=1 and must report the bit-for-bit difference.
    foreach(solver IN ITEMS closest_point_perzyna closest_point_time)
        muscle_reference_test(materials-mat_${solver})
        add_test(NAME materials-mat_${solver}_reference_mutation COMMAND test_mat_${solver} --mutate rel)
        muscle_reference_test(materials-mat_${solver}_reference_mutation)
        set_tests_properties(materials-mat_${solver}_reference_mutation PROPERTIES
                             PASS_REGULAR_EXPRESSION "outside the tolerance of its category"
                             FAIL_REGULAR_EXPRESSION "no reference for|not produced by the test|cannot open|malformed")
        add_test(NAME materials-mat_${solver}_reference_mutation_strict COMMAND test_mat_${solver} --mutate ulp)
        muscle_reference_test(materials-mat_${solver}_reference_mutation_strict ENV MUSCLE_REFERENCE_STRICT=1)
        set_tests_properties(materials-mat_${solver}_reference_mutation_strict PROPERTIES
                             PASS_REGULAR_EXPRESSION "differs from the reference value \\(bit-for-bit"
                             FAIL_REGULAR_EXPRESSION "no reference for|not produced by the test|cannot open|malformed")
    endforeach()
    # Physical checks of the anisotropic yield criteria with closest_point (Lankford coefficients,
    # tension/compression, shears, large increments), ported from the OO solver tests (D-018.18g)
    add_muscle_fypp_test(mat_cp_anisotropic test_muscle_material_cp_anisotropic.fypp
                         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mat_sequence_common.fypp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/case_params.py
                                 ${CMAKE_CURRENT_SOURCE_DIR}/anisotropic_cases.py)
    # Physical checks of the return mapping kernels (closed forms, iteration limits, FE history
    # boundary), ported from the OO solver tests (D-018.18g)
    add_muscle_fypp_test(mat_solver_checks test_muscle_material_solver_checks.fypp
                         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mat_sequence_common.fypp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/case_params.py)
    # Physical checks of the Perzyna kernel (flow rule, exact relaxation, eta -> 0, static limit of the
    # anisotropic criteria), ported from the OO solver tests (D-018.18g)
    add_muscle_fypp_test(mat_perzyna_checks test_muscle_material_perzyna_checks.fypp
                         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mat_sequence_common.fypp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/case_params.py)
    # Return mapping kernels of every solver against independent truths (plan R, R6): closed-form radial
    # return and tangents, elastic threshold, homogeneity in the stress units, the equations of a converged
    # step, quadratic convergence of the Newton solvers
    add_muscle_fypp_test(mat_return_mapping test_muscle_material_return_mapping.fypp
                         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mat_sequence_common.fypp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/case_params.py)
    # Singular Newton system of closest_point_2 (own LU, D-018.1)
    add_muscle_fypp_test(mat_cp2_singular test_muscle_material_cp2_singular.fypp)
    # Error paths of the time step of a viscoelastic material and a model with fewer branches than the
    # history array keeps (Phase 7 step 6)
    add_muscle_fypp_test(mat_visco_errors test_muscle_material_visco_errors.fypp)
    # Finite strain materials (D-018.15): physical oracles (tangent vs finite differences of the
    # Kirchhoff stress, objectivity, F -> I, uniaxial closed form)
    add_muscle_fypp_test(mat_kinematics test_muscle_material_kinematics.fypp)
    # Batch API muscle_update_<name> (D-008.2, D-018.7): bit for bit the point kernel point by point
    # (every solver, frozen viscoelasticity, every kinematics; tangent, index list, 1 and 4 threads,
    # rejected time step) and its thread safety (the historical test of the OO solvers, test_muscle_openmp)
    add_muscle_fypp_test(mat_batch test_muscle_material_batch.fypp)
    add_muscle_fypp_test(mat_batch_openmp test_muscle_material_batch_openmp.fypp)
    # Batch API on device memory muscle_update_<name>_dev (D-019.7a): arrays of omp_target_alloc against
    # muscle_update_<name> on mapped arrays, one test per material and kind (target_cases.py)
    add_muscle_fypp_target_tests(mat_batch_dev test_muscle_material_batch_dev.fypp batch_dev)
    # Point kernels of materials that use every model and every solver in an OpenMP target region
    # against the CPU (testing.md, GPU), one test per material (add_muscle_fypp_target_tests,
    # tests/CMakeLists.txt; the materials are in tests/gpu/target_cases.py)
    add_muscle_fypp_target_tests(mat_kernels_target test_muscle_material_kernels_target.fypp material
                                 DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mat_sequence_common.fypp
                                         ${CMAKE_CURRENT_SOURCE_DIR}/case_params.py)
    # Generic API for FE codes (D-018.17): point API over flat arrays bit for bit the point kernel,
    # rotation of the history array
    add_muscle_fypp_test(mat_point test_muscle_material_point.fypp)
    # Material initializers from a flat parameter array (D-018.17a, e): every material of the registry,
    # positions of the parameter table, same behaviour as the material built by hand, error statuses
    add_muscle_fypp_test(mat_init test_muscle_material_init.fypp
                         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/props_samples.py
                                 ${CMAKE_CURRENT_SOURCE_DIR}/case_params.py)
endif()
