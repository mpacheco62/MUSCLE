title: Adding New Models

# Developer Guide

Models are not added by extending a base class: there are no abstract types nor type-bound model
procedures. A model is a **concrete parameter type** plus **free `pure` functions** with the signature of
its family, written once as an FYPP template for both kinds (`sp`, `dp`) and, for tensor families, for
both dimensions. Solvers call the family generics (`hard_stress`, `yield_stress_eq`, ...), which resolve at
compile time from the type of the parameters, so a new model is usable by every solver of its family
without touching the solvers. The rules every contribution follows are in [Conventions](conventions.html);
the [architecture](architecture.html) page shows how a material is put together and the
[glossary](glossary.html) explains the abbreviations of the names.

A model is written **by hand, in Fortran, in one file**: the parameter type, the signatures with their `intent`
and documentation, the bodies and the closings are in the template. FYPP only instantiates the text for each
kind (`#:for k in REAL_KINDS`, `${k}$`) and marks the kernels (`@:DEVICE()`), which is what Fortran cannot do
(rule «minimal FYPP», see [Conventions](conventions.html#minimal-fypp)). No model declaration is built from metadata. The
only generated modules are the materials, one per entry of the registry (type, parameter array and APIs of each
material), because Fortran cannot declare one module per entry; they are written once in `src/materials/`.

`tools/fypp/muscle_meta.py` holds **data**, used where a material is generated (the parameter array of
`mat_init_<name>`) and by the tests:

| Registry | Content |
|----------|---------|
| `MODEL_PARAMS[(prefix, model)]` | Names and order of the components of the parameter type, with their doc. The generated `mat_init_<name>` assigns them by name (`mat%yield%sigy = props(3)`), so the compiler rejects a name that the type does not have; `tests/tools/test_muscle_model_params.py` also checks the order, the types and the defaults against the type written in the template. |
| `MODEL_INITS[(prefix, model)]` | Inputs of the initializer `<prefix>_<model>_init`, when the parameters are derived (e.g. Lamé constants from E, ν) or validated, with their domains. Every initializer ends with `integer, intent(out) :: status` (`STATUS_OK` or an error code of `muscle_status`; always `STATUS_OK` if it validates nothing), registered as `status=`. |
| `HARDENING_MODELS`, `YIELD_MODELS`, `VISCOPLASTIC_MODELS`, ... | The registered models of each family (the materials look them up here). |
| `MODEL_LABELS` | Readable name of the model for the generated docs. |
| `MATERIAL_SOLVERS` | Solvers usable by a material, the families each one needs and its history size. |

The order of the props of `mat_init_<name>` comes from `MODEL_PARAMS` / `MODEL_INITS`
(`python3 tools/fypp/muscle_meta.py props <material>`).

## How to add a model to an existing family

Where each family lives (paths relative to `src/` and `tests/`; `<m>` is the model name):

| Family | Template | Family module | Lists in `tools/fypp/muscle_meta.py` (besides `MODEL_PARAMS`, `MODEL_LABELS`) | Material solvers that use it | Test of the model | Mutation area |
|--------|----------|---------------|------------------------------------------------------------------------------|------------------------------|-------------------|---------------|
| `hard` | `plasticity/hardening_laws/muscle_hard_<m>.fypp` | `muscle_hard_laws` | `HARDENING_MODELS` | `closest_point`, `closest_point_2`, `cutting_plane` (and composed by `vp`/`os` laws) | `test_muscle_hard_laws_kinds` | `hard` |
| `vp` | `plasticity/viscoplastic_laws/muscle_vp_<m>.fypp` | `muscle_vp_laws` | `VISCOPLASTIC_MODELS`; `VISCOPLASTIC_COMPOSED` if it composes a hardening law; `VISCOPLASTIC_FLOOR` if it has a rate floor | `closest_point_time` | `test_muscle_vp_laws_kinds` (boundaries: `test_muscle_vp_laws_limits`) | `vp` |
| `os` | `plasticity/viscoplastic_laws/muscle_os_<m>.fypp` | `muscle_os_laws` | `OVERSTRESS_MODELS`, `OVERSTRESS_COMPOSED` | `closest_point_perzyna` | `test_muscle_vp_laws_kinds` (Chaboche: `test_muscle_os_chaboche`) | `vp` |
| `yield` | `plasticity/yield_criteria/muscle_yield_<m>.fypp` | `muscle_yield_laws` | `YIELD_MODELS`; `YIELD_JOINT` if it writes its own joint evaluation; `YIELD_ISOTROPIC` if isotropic | all five | `test_muscle_yield_laws_kinds` | `yield` |
| `elas` | `elasticity/muscle_elas_<m>.fypp` (`frozen`: `viscoelasticity/muscle_elas_frozen.fypp`) | `muscle_elas_laws` | `ELASTICITY_MODELS`, `ELASTICITY_ISOTROPIC` | all five | `test_muscle_elas_laws_kinds` | `elas` |
| `visco` | `viscoelasticity/muscle_visco_<m>.fypp` | `muscle_visco_laws` | `VISCOELASTIC_MODELS` | all, through `elas='frozen'` | `test_muscle_visco_laws_kinds` or `test_muscle_visco_<m>` | `visco` |
| `kin` | `kinematics/muscle_kin_<m>.fypp` | `muscle_kin_laws` | `KINEMATICS_MODELS`; `KINEMATICS_SPATIAL`, `KINEMATICS_ROTATION`, `KINEMATICS_FE_ROTATED` as the case may be | all, through `kin=` | `test_muscle_kin_laws_kinds`, `test_muscle_kin_finite` | `kin` (`mat_kin` adds the kinematic materials) |

Besides: `MODEL_INITS` only if the model has an initializer of its own (derived parameters or validation), and
`MODEL_DOMAINS` if `mat_init` must check a condition between parameters (today only Johnson-Cook). The device
test of the family (`tests/gpu/`) takes every model of the family by itself, with the parameter sets of
`CASE_PARAMS` (step 6). Add the model to the list of models of the `README.md` (nothing checks it).

Example: a new hardening law `mylaw` with \( \sigma_y = \sigma_{y0} + K\epsilon_p^2 \).

1. **Write the template** `src/plasticity/hardening_laws/muscle_hard_mylaw.fypp` (copy
   `muscle_hard_bilinear.fypp`). The type, the three functions of the family contract, their `intent`s
   and their closings are written out:

   ```fortran
   #:include "muscle_common.fypp"

   module muscle_hard_mylaw
       !! Quadratic hardening law ...
       use muscle_kinds, only : ${', '.join(REAL_KINDS)}$
       implicit none
       private

   #:for k in REAL_KINDS
       public :: hard_mylaw_${k}$
       public :: hard_stress_mylaw_${k}$, hard_dstress_dep_mylaw_${k}$, hard_ddstress_ddep_mylaw_${k}$
   #:endfor

   #:for k in REAL_KINDS
       type :: hard_mylaw_${k}$
           !! Parameters of the quadratic model in `real(${k}$)`.
           real(${k}$) :: y0
               !! Initial yield stress \( \sigma_{y0} \).
           real(${k}$) :: K
               !! Hardening coefficient \( K \).
       end type hard_mylaw_${k}$

   #:endfor
   contains
   #:for k in REAL_KINDS
       pure function hard_stress_mylaw_${k}$(prm, ep) result(res)
           !! Flow stress \( \sigma_y(\epsilon_p) \) of the quadratic model in `real(${k}$)`.
           implicit none
           @:DEVICE()
           type(hard_mylaw_${k}$), intent(in) :: prm
               !! Parameters.
           real(${k}$), intent(in) :: ep
               !! Equivalent plastic strain \( \epsilon_p \ge 0 \).
           real(${k}$) :: res
               !! \( \sigma_y \).
           res = prm%y0 + prm%K*ep**2
       end function hard_stress_mylaw_${k}$

       ! ... hard_dstress_dep_mylaw_${k}$ and hard_ddstress_ddep_mylaw_${k}$, written the same way, with
       ! analytical derivatives (the library has no finite-difference code: the macros of
       ! tests/common/muscle_fd.fypp are for the tests that check these derivatives)
   #:endfor
   end module muscle_hard_mylaw
   ```

   Every kernel procedure is `pure` and device-safe (no `class`, `allocatable`, I/O, `stop`, LAPACK) and
   carries `@:DEVICE()` right after `implicit none`; literals and constants carry the kind (`1.0_${k}$`,
   constants of `muscle_kinds`).

2. **Add it to the family module** `muscle_hard_laws.fypp`: one `use` statement of its own (do not append the
   names to the `use` of another law) and one name in each of the three
   `interface` blocks (`hard_stress`, `hard_dstress_dep`, `hard_ddstress_ddep`) and in `public`. The
   compiler rejects a name that does not exist.

3. **Register it for the materials** in `tools/fypp/muscle_meta.py`: add `'mylaw'` to `HARDENING_MODELS`, its
   parameters to `MODEL_PARAMS` (same names and order as the type, a test checks it) and a label to
   `MODEL_LABELS`:

   ```python
   ('hard', 'mylaw'): (
       _P('y0', 'real', '', r'Initial yield stress \( \sigma_{y0} \).', None),
       _P('K', 'real', '', r'Hardening coefficient \( K \).', None),
   ),
   ```

4. **Add it to the build**: list the template in `muscle_add_fypp_sources` of the folder's
   `CMakeLists.txt` (before the family module; the order among the laws does not matter).

5. **Tests** (`tests/plasticity/hardening_laws/test_muscle_hard_laws_kinds.fypp`): add a function
   `test_mylaw_${k}$()` next to those of the other laws and register it with `call record(...)` in the main
   program. Compare the values with a `dp` reference computed in the test, and the derivatives with finite
   differences using `check_scaled` (relative to the size of the function: a plain relative check fails in
   `sp` where a derivative is close to zero). The law runs inside an existing ctest test, so the number of
   tests does not change. `tests/tools/test_muscle_model_params.py` checks the metadata by itself. Then run
   the [mutation check](#tests) of the family (its area is in the table above: `tools/mutation_check.py hard`).
   For a new file, list its mutations and run only those (about a minute, instead of 30–45 minutes for a whole
   area): `tools/mutation_check.py hard --samples all --list` gives the ids (`kind:file:line:column`), and
   `tools/mutation_check.py hard --samples all --only <id>,<id>,...` runs them. A mutation that changes a result
   but survives is a gap of your tests.

6. **Make it reachable from the materials** of the test registry `tools/fypp/muscle_materials.py`, which
   must cover every model with **every solver that uses its family** (`tests/tools/test_muscle_materials.py`
   checks it; the solvers of each family are in the table above; for a hardening law: `closest_point`,
   `closest_point_2` and `cutting_plane`):

   ```python
   Material('cp_vm_mylaw', 'closest_point', elas='linear', yield_='vonmises', hard='mylaw',
            doc='von Mises with quadratic hardening'),
   Material('cp2_vm_mylaw', 'closest_point_2', elas='linear', yield_='vonmises', hard='mylaw'),
   Material('cut_vm_mylaw', 'cutting_plane', elas='linear', yield_='vonmises', hard='mylaw'),
   ```

   Every material of the test registry runs in the strain sequence test of its solver, so also:
   - add parameter sets of the model to `CASE_PARAMS` in `tests/materials/case_params.py`
     (copy the entry of a neighbour law);
   - add a `solver_case(0, ...)` for each new material to `CASES` in `tests/materials/test_muscle_material_<solver>.fypp`
     (fields in `tests/materials/mat_sequence_common.fypp`: case number `0` (no reference records), the material,
     the elasticity `'steel'` (E = 210000) or `'lin1000'` (E = 1000), the parameter set names of the models,
     `iter_nw`, `dt` and the yield strain `a_y` of the strain sequence, about y0/E; copy a neighbour case). The
     test checks sp against dp and the tangent against finite differences along the sequence; no reference
     records are written for a new material.

   Your own materials go in a registry of your own, passed with `-DMUSCLE_MATERIALS=<path>`.

The same steps apply to the other families, with the names of the table above (copy the neighbour model of
the family): a yield criterion
writes `stress_eq` and its gradient and Hessian (and may write the joint evaluation by hand, as the
existing criteria do); a viscoplastic or overstress law that composes a hardening law receives its values
(`sy`, `dsy_dep`) as arguments.

## How to add a solver

A solver is a body `src/solvers/muscle_<solver>_body.fypp` included by the material template
`src/materials/muscle_material.fypp` with the kind fixed, plus a `MaterialSolver` entry in
`MATERIAL_SOLVERS` (families it calls, history size, whether it uses `dt`, default iteration limit). It
only calls family generics, so it works with every registered model of those families. Its body declares
the common dummy arguments of `SOLVE` by hand (copy them from another solver body). Add its
`#:include` to the solver chain of `muscle_material.fypp` (fypp only takes literal file names; the CMake
dependencies come from `MATERIAL_SOLVERS`, and `tests/tools/test_muscle_materials.py` checks that the chain
matches it), and tests in `tests/materials/`: a strain sequence test `test_muscle_material_<solver>.fypp` and its cases in
`test_muscle_material_return_mapping.fypp` (closed-form radial return, consistency at convergence, tangent against finite
differences). Add the body to `tests/mutations/solvers.toml` and run the mutation check; keep reference
records only for what no independent truth can detect (the Newton path, see `test_muscle_material_closest_point_time.fypp`).

## Checklists

**New model of an existing family**
- [ ] Template (the path of its family in the [table](#how-to-add-a-model-to-an-existing-family)) written by hand
      (type, functions, `@:DEVICE()`), listed in the
      folder's `CMakeLists.txt` before the family module.
- [ ] Family module: its own `use` statement and the specifics in every `interface` and in `public`.
- [ ] `tools/fypp/muscle_meta.py`: the lists of the family (table), `MODEL_PARAMS` (same names and order as the
      type), `MODEL_INITS` if it has an initializer, `MODEL_LABELS`; the model in the list of the `README.md`.
- [ ] Tests in `sp` and `dp` against independent truths (closed forms, derivatives against finite
      differences); `tools/mutation_check.py <area> --samples all` (area of the table; `--only` with the ids of
      the new file) with no surviving mutation that changes a
      result.
- [ ] One material per solver of the family in `tools/fypp/muscle_materials.py`, with `CASE_PARAMS` and a
      `solver_case(0, ...)` in the sequence test of each solver.
- [ ] `tools/check.sh --quick` while working; `tools/check.sh` before closing.

**New material of existing models** (test registry): the last two items above.

**New solver**: body, `MATERIAL_SOLVERS` entry, `#:include` in the solver chain of `muscle_material.fypp`,
strain sequence test `tests/materials/test_muscle_material_<solver>.fypp`, its cases in `test_muscle_material_return_mapping.fypp`, the
mutation check, and at least one material per model of its families in the test registry.

## Naming of the generated material API

Every material `<name>` of the registry gives a module `muscle_mat_<name>`. Its public names follow one
convention:

| Prefix | Meaning | Names |
|--------|---------|-------|
| `muscle_*` | User-facing entry points of one material, the calls an FE code or a driver makes per step | `muscle_update_<name>` (batch API), `muscle_point_<name>` (point API, with or without `ddsdde`) |
| `mat_*` | Per-material helpers and the generics shared by every material module | `mat_init_<name>`, `mat_pack_<name>`, `mat_unpack_<name>`, `mat_rotate_hsv_<name>`, the type `mat_<name>_<k>`, the generics `mat_update` and `mat_step` |
| `MAT_<NAME>_*` | Constants of one material | `MAT_<NAME>_NSTATE`, `MAT_<NAME>_NPROPS` |

The specific procedures behind these generics are `mat_<proc>_<name>_<k>` (`muscle_meta.material_proc`)
and are private.

## Tests

`tests/` mirrors `src/` (one folder per area), and so do the file names: the tests of
`src/<path>/muscle_<x>.fypp` are `tests/<path>/test_muscle_<x>[_<topic>].fypp` (e.g. `test_muscle_hard_laws_kinds`,
`test_muscle_material_point`; a test of no single module is named after its topic, `test_muscle_fypp_smoke`). The
ctest name is the folder path joined with `-` plus the name given to `add_muscle_fypp_test`
(`plasticity-hardening_laws-hard_laws_kinds`). `tests/common/` holds the shared assertions
(`muscle_test_utils`), helpers and finite-difference macros, `tests/tools/` the Python tests of the FYPP
metadata and `tests/gpu/` the kernel tests of every model family on the device (`test_muscle_family_target.fypp`,
one test per family and kind, cases in `target_cases.py`). Every ctest test of the Fortran code carries a label
(the Python tests of `tests/tools/` and the compile-fail checks `*_private_params` have none):

| Label | Tests | Run with |
|-------|-------|----------|
| `kinds` | templates instantiated for `sp` and `dp` against independent truths (closed forms, properties, finite differences, another computation path) | `ctest -L kinds` |
| `reference` | the few comparisons with stored values of `tests/reference/` (below) | `ctest -L reference` |
| `target` | kernels inside an OpenMP target region against the CPU | `ctest -L target` |

**Unit tests first.** A test needs an independent truth: a closed form, a property (symmetry, homogeneity,
`A A^-1 = I`), finite differences of the function, or the same result through another path. A *reference*
test only detects a change: it compares with values stored from a validated version
(`tests/reference/<area>/*_dp.txt`, written by `regenerate.sh` there, compared within the tolerance of their
category, `REFERENCE_TOL_*` of `tests/common/muscle_test_reference.f90`; `MUSCLE_REFERENCE_STRICT=1` compares bit
for bit and `MUSCLE_REFERENCE_REPORT=1` prints the largest error of each record). Three reference tests
remain: the Newton path of `closest_point_time` and `closest_point_perzyna` and three results of the math module
(`math_reference`); the header of each file says why. The other four tests labelled `reference`
(`*_reference_mutation[_strict]`) check that the comparison itself detects a perturbation.

**Mutation check.** `tools/mutation_check.py <area>` applies mutations to the sources of an area (an operator
`+`/`-`/`*`/`/` swapped, a literal changed, a subscript shifted, `<` relaxed; sampled, or all of them with
`--samples all`), rebuilds a `dp` Debug build and reports which ones the unit tests and the reference tests
detect. The areas are described in `tests/mutations/<area>.toml`. A reference test is deleted only if the unit
tests detect every mutation it detects; for new code, a mutation that changes a result but survives is a test
missing.

## Development cycle

Typical costs on a machine with 8 or more cores (gfortran, Ninja, `sp;dp`; orders of magnitude):

| When | What to run | Cost |
|------|-------------|------|
| While editing a model | `cmake --build build --target build_all_tests` and `ctest --test-dir build -R <area>` | ~1 s to rebuild after a change in the body of a function (nothing else recompiles); ~6 s after a change in the interface of a law (the material modules recompile) |
| After touching `tools/fypp/muscle_meta.py`, `include/muscle_common.fypp` or `muscle_kinds` | the same | a full rebuild, ~1 min (a new model always touches `muscle_meta.py`) |
| Every test | `ctest --test-dir build -j8` | ~10 s |
| A new or changed file | `tools/mutation_check.py <area> --samples all --only <ids>` | ~2 s per mutation, about a minute for a file (a whole area: 30–45 min) |
| Before every commit | `tools/check.sh --quick` | ~1–1.5 min |
| To close a piece of work | `tools/check.sh` | ~15 min, mostly the three GPU builds and the target tests on the Intel GPU |

The full check needs the offload toolchains (gfortran with NVPTX and AMDGCN offload, ifx) and an Intel GPU; without
them, `--quick` is the gate.

## Verification

`tools/check.sh --quick` (about a minute) lints the mechanical rules (kinds, literals, forbidden constructs in
kernels, SPDX) and builds and tests `sp;dp` in Debug: run it while developing. The build must have **no compiler
warnings** (library and tests, `-Wall -Wextra`; unused dummy arguments are not reported, since a family keeps a
uniform interface, but unused local variables are): `check.sh` fails on any warning and prints its template line
(with `--quick`, of the files it recompiled; the full check rebuilds from scratch). `tools/check.sh` without
`--quick` also builds the formatted-output (UDTIO) path, the FORD documentation and the three GPU targets,
runs the target tests on an Intel GPU and audits the device images: run it before closing a piece of work. See also the [migration guide](migration.html) for the mapping from the
former object-oriented API.

## Editor support

Language servers do not run FYPP: on a `.fypp` template, [fortls](https://fortls.fortran-lang.org) only sees
truncated names (`hard_stress_swift_` instead of `hard_stress_swift_${k}$`) and cannot resolve definitions. Point
it at the generated code instead, which is plain Fortran:

1. Build once (`tools/check.sh --quick`, or `cmake --build <build_dir> --target build_all_tests`).
2. Put a `.fortlsrc` at the repository root (not versioned: the build directory is yours), e.g.

   ```json
   {
     "source_dirs": ["build-check/fypp/src/**/"],
     "lowercase_intrinsics": true
   }
   ```

3. Open the generated file of the template you work on, `tools/generated.sh <template>` prints its path
   (the material template and the solver bodies give the module of every declared material). "Go to definition"
   and hover then work across modules; edit the template, not the generated file.

Compiler errors and coverage already point to the lines of the templates (`#line` markers).
