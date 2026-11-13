title: Architecture

# Architecture

How the pieces of MUSCLE fit together, from the tensor types to the API that an FE code calls. See the
[glossary](glossary.html) for the abbreviations and [Conventions](conventions.html) for the rules.

## Layers

```
 FE code (UMAT, USERMAT, umat41, own solver)
        │  flat arrays: strain(6) or F_n(9), F_np1(9), props(:), hsv(:)
        ▼
 declared material  muscle_mat_<name>          ← generated, one module per entry of the registry
   mat_init_<name>, muscle_point_<name>, muscle_update_<name>[_dev], mat_update, mat_step, ...
        │  calls the solver body with the components of the material
        ▼
 solver body  src/solvers/muscle_<solver>_body.fypp   ← one text, included in every material module
        │  calls the family generics only: elas_stress, yield_stress_eq, hard_stress, kin_update, ...
        ▼
 model families  muscle_<family>_laws / muscle_yield_laws   ← generic interfaces, written by hand
        │  resolve at compile time by the type of the parameters (hard_swift_dp, yield_vonmises_dp, ...)
        ▼
 models  muscle_<family>_<model>                ← one file per model: type + pure functions
        │
        ▼
 tensors and math  muscle_tensors, muscle_math_*   ← types ten_<dim><order><sym>_<k>, operators, LU, spectral (math.html)
```

Every layer below the material is written once for both kinds (`sp`, `dp`) and, for the tensors, for 2D and
3D. There is no dynamic dispatch anywhere: a material module only differs from another one in the types of its
components, so each call resolves at compile time and the whole point kernel can run on a GPU.

## From the registry to a material module

1. **Registry.** `tools/fypp/muscle_materials.py` (or your own file, `-DMUSCLE_MATERIALS=<path>`) lists the
   materials: name, solver and one model per family the solver needs.
   `Material('cp_vm_swift', 'closest_point', elas='linear', yield_='vonmises', hard='swift')`.
   `muscle_meta.load_materials` validates it at configure time.
2. **Generation.** CMake (`muscle_add_materials` in `cmake/MuscleFypp.cmake`) preprocesses the template
   `src/materials/muscle_material.fypp` once per entry, with the define `MUSCLE_MATERIAL=<name>`, into
   `<build>/fypp/src/materials/muscle_mat_<name>.f90`.
3. **What the template writes**, per kind:
   - the material type `mat_<name>_<k>`: one component per model, plus the solver settings `tol` and `iter_nw`;
   - the solver body, by `#:include` of `src/solvers/muscle_<solver>_body.fypp` (and the shared consistent
     tangent of the closest point solvers);
   - the point kernel `mat_update` and the history layout (`mat_pack`/`mat_unpack`):
     `src/materials/muscle_material_kernel.fypp`;
   - the batch API `muscle_update_<name>`: `src/materials/muscle_material_batch.fypp`;
   - the FE API `mat_init_<name>` (from `props`) and `muscle_point_<name>`: `src/materials/muscle_material_point.fypp`.
4. **Where the variants are.** The material template branches on the few things that change the code of a
   material: a viscoelastic elasticity (the kernel freezes and updates the viscoelastic state), finite strain
   kinematics (the kernel takes deformation gradients) and the corotational kinematics. Everything else is the
   same text with other type names.

## What the generated code looks like

Excerpt of `<build>/fypp/src/materials/muscle_mat_cp_vm_swift.f90` (double precision, comments removed):

```fortran
type :: mat_cp_vm_swift_dp
    type(elas_linear_dp) :: elas
    type(yield_vonmises_dp) :: yield
    type(hard_swift_dp) :: hard
    integer :: iter_nw = 200
    real(dp) :: tol = NEWTON_TOL_dp
end type mat_cp_vm_swift_dp

pure subroutine mat_update_cp_vm_swift_dp(mat, strain, dt, state_n, state_np1, stress, tangent, status, iters)
    ...
    call mat_solve_cp_vm_swift_dp(mat, mat%elas, strain, dt, state_n, state_np1, stress, tangent, status, iters)
end subroutine mat_update_cp_vm_swift_dp

! inside mat_solve_cp_vm_swift_dp: the closest point body, the same for every material of this solver
stress_trial = elas_stress(elas, strain - state_n%strain_p)
hard_n = hard_stress(mat%hard, state_n%strain_pf)
f_trial = yield_stress_eq(mat%yield, stress_trial) - hard_n
if (f_trial/max(hard_n, SY_MIN_dp) <= TOL_ELASTIC_dp) then
    ...
```

`elas_stress`, `hard_stress` and `yield_stress_eq` are the family generics: here they resolve to
`elas_stress_linear_3D_dp`, `hard_stress_swift_dp` and `yield_stress_eq_vonmises_3D_dp`. To read the complete
module of any material, build once and open `<build>/fypp/src/materials/muscle_mat_<name>.f90`; compiler errors
and coverage point to the lines of the templates (`#line` markers).

## Where to change what

| To ... | Edit |
|--------|------|
| add a model | its template, the family module, `tools/fypp/muscle_meta.py` and the tests ([Developer Guide](contributing.html)) |
| add a material | the registry (and, for the test registry, a case of the sequence test of its solver) |
| add a solver | a body in `src/solvers/`, `MATERIAL_SOLVERS`, the solver chain of `muscle_material.fypp` |
| change the API of every material | `src/materials/muscle_material*.fypp` |
| change a tensor operation | `src/tensors/operators/`, whose bodies are unrolled from the layouts of `muscle_meta.py` |
