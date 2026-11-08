title: Migration from the OO API

# Migration from the object-oriented API

Up to the FYPP migration MUSCLE had an object-oriented API: abstract base types (`Base_elasticity`,
`Base_hardening_laws`, `Base_yield_critera`, `Base_kinematics`, ...), models as type extensions with
type-bound procedures, solver objects holding polymorphic (`class`, `allocatable`) copies of the models, and
a transactional history container. That layer cannot be compiled for GPUs (dynamic dispatch, allocatable
polymorphic components, procedure arguments) and has been removed (decisions D-007, D-008 and D-018 in
`docs/fypp/DECISIONES.md`). There is no compatibility layer: code that used it is rewritten against the
declared materials (see [FEA Integration](integration.html)) or against the free functions of the model
families. The tensor API (`muscle_tensors`, the type-bound `%init()`, `%xx()`, ... of the tensors and the
formatted output of `muscle_tensor_io`) is unchanged.

## Equivalences

| Old OO API | New API |
|------------|---------|
| Combination of models chosen at run time (`class(Base_...)` components of a solver) | Declared material: one `Material(...)` entry in a registry (`tools/fypp/muscle_materials.py` or `-DMUSCLE_MATERIALS=<path>`), generated module `muscle_mat_<name>` |
| `type(Closest_point) :: solver` + `call solver%init(elasticity, hardening, yield, iter_nw, tol)` | `type(mat_<name>_dp) :: mat` (or `_sp`) + `call mat_init_<name>(mat, props, status)`, or set its components (`mat%elas`, `mat%yield`, `mat%hard`, ..., `mat%tol`, `mat%iter_nw`) |
| Solver modules `muscle_solver_closest_point`, `_closest_point_perzyna`, `_closest_point_time`, `_closest_point_2`, `_cutting_plane` (types `Closest_point`, `Closest_point_2`, ...) | Solver of the material entry: `'closest_point'`, `'closest_point_perzyna'`, `'closest_point_time'`, `'closest_point_2'`, `'cutting_plane'` |
| `call solver%solve(strain, history, status, iters)` | `call muscle_point_<name>(mat, strain, dt, hsv_n, hsv_np1, stress, status, iters)` (flat arrays), or the kernel `call mat_update(mat, strain, dt, state_n, state_np1, stress, tangent, status, iters)` (tensor types) |
| Loop of `solver%solve` over the Gauss points (OpenMP `private` solver copies) | `call muscle_update_<name>(mat, eps, dt, hsv_n, hsv_np1, sig, status, iters [, ddsdde] [, idx])` over SoA arrays `eps(npts, 6)` |
| `call solver%tangent(strain, history, tangent)` (`ten_3D4O2sym`) | Optional `ddsdde` of the point / batch API (`ddsdde(6, 6)`, engineering shear columns as UMAT DDSDDE), or `tangent` of `mat_update` |
| `solver%tangent_numerical(...)` | No kernel equivalent: finite differences of the point API in tests (macros of `tests/common/muscle_fd.fypp`) |
| Free wrapper `closest_point2(strain, elasticity, hardening, yield, stress, strain_pf, strain_p, status)` | `muscle_point_<name>` of a `'closest_point_2'` material |
| `type(Plastic_material_history)`: `%init`, `%unpack_from_fea`, `%pack_to_fea` | History array of `MAT_<NAME>_NSTATE` reals per point (`hsv_n`, `hsv_np1`), `mat_pack_<name>` / `mat_unpack_<name>` for the kernel states; the stress is not part of it |
| `history%commit()` / `history%rollback()` | None: the caller (FE code) keeps the converged `hsv_n`; the update writes only `hsv_np1` |
| `Viscoelastic_history`, `Visco_coefficients` | State inside the material history; `call mat_step(mat, dt, coefs, status)` once per step for materials with `elas='frozen'` |
| `type(Elasticity_linear)` + `%set_parameters(young=, poisson=)` | `type(elas_linear_dp)` + `call elas_linear_init(prm, young, poisson, status)` (`muscle_elasticity_linear`) |
| `elasticity%stress(strain)`, `%dstress_dstrain(...)`, `%stress_3D` / `%stress_2D` | `elas_stress(prm, strain)`, `elas_dstress_dstrain(prm, strain)`, ... (`muscle_elasticity_laws`; generic over 3D/2D and kind) |
| `Elasticity_frozen`, Maxwell / `Kelvin_voigt` objects | `elas_frozen_<k>`, `visco_maxwell_<k>`, `visco_kelvin_<k>`, functions `visco_step_coefficients`, `visco_freeze`, `visco_update` (`muscle_visco_laws`) |
| `type(Swift_hardening)` (and the other `Base_hardening_laws` extensions), `%stress(ep)`, `%dstress_dep(ep)` | `type(hard_swift_dp)` with public components (`k`, `e0`, `n`), `hard_stress(prm, ep)`, `hard_dstress_dep(prm, ep)`, `hard_ddstress_ddep(prm, ep)` (`muscle_hard_laws`) |
| `type(VonMises)`, `type(Hill48)`, ... (`Base_yield_critera`), `%stress_eq(sigma)` and derivatives | `yield_vonmises_dp`, `yield_hill48_dp`, ... (`yield_<m>_init` for the anisotropic ones and Drucker-Prager), `yield_stress_eq`, `yield_dstress_eq`, `yield_ddstress_eq`, `yield_stress_eq_and_derivatives` (`muscle_yield_criteria`) |
| `Perzyna_overstress`, `Chaboche_overstress`, `JohnsonCook_full_visco_hardening`, ... (with a `class` hardening component) | `os_perzyna_<k>`, `os_chaboche_<k>`, `vp_jc_<k>`, ... with the composed hardening law as a separate argument (`sy`, `dsy_dep`): `os_flow_rate`, `vp_flow_stress`, ... (`muscle_vp_overstress_laws`, `muscle_vp_laws`) |
| Facade `muscle_kinematics` (`Small_strain_kinematics`, `Total_lagrangian_kinematics`, ...) | `kin='<kinematics>'` (and `tangent_rate`) in the material entry; the point API takes `F_n(9)`, `F_np1(9)`; family functions `kin_update`, `kin_to_cauchy`, ... in `muscle_kin_laws` |
| `derivative(fun, x[, eps])`, `derivative2O(...)` (`procedure(...)` argument, `*_host` modules) | Macros `FD_DERIV`, `FD_DERIV2`, ... of `tests/common/muscle_fd.fypp` (FYPP, call expanded in place) |
| Solver status `STATUS_ELASTIC_CASE`, `STATUS_CONVERGED`, `STATUS_NONCONVERGED` | `STATUS_OK` (elastic if `iters = 0`) and `STATUS_NOT_CONVERGED` of `muscle_status` |
| Double precision only | Every type and procedure in `_sp` and `_dp` |

## What is lost

- **Transactional commit / rollback** of the history (`Plastic_material_history%commit`, `%rollback`):
  the FE code keeps the converged state `hsv_n` and passes it again to retry an increment; `hsv_n` and
  `hsv_np1` must be different arrays.
- **Run-time choice of the model**: the combination of models of a material is fixed at build time (only
  its parameters are read at run time). Several materials are several registry entries, and a multiphase
  model calls one batch per material with the index list `idx` of its phase. A tagged union with a
  generated `select case` remains an optional future extension.
- **User models by class extension**: a new model is added to the library sources as an FYPP template
  registered in `muscle_meta` (see the [Developer Guide](contributing.html)), not by extending a base type
  in user code.
- **`set_parameters` and `to_string`** of the models: parameters are public components of the parameter
  types, set directly or by `<prefix>_<model>_init` / `mat_init_<name>`; there is no textual description
  of a model (the parameter table of a material is printed by
  `python3 tools/fypp/muscle_meta.py props <name>`).
