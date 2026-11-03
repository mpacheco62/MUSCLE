title: Adding New Models

# Developer Guide

Models are not added by extending a base class: there are no abstract types nor type-bound model
procedures. A model is a **concrete parameter type** plus **free `pure` functions** with the signature of
its family, written once as an FYPP template for both kinds (`sp`, `dp`) and, for tensor families, for
both dimensions. Solvers call the family generics (`hard_stress`, `yield_stress_eq`, ...), which resolve at
compile time from the type of the parameters, so a new model is usable by every solver of its family
without touching the solvers. Rules: `.claude/rules/` (templates, kinds, GPU-safe code) and the
decisions D-017 / D-018 in `docs/fypp/DECISIONES.md`.

The single source of truth is `tools/fypp/muscle_meta.py`:

| Registry | Content |
|----------|---------|
| `MODEL_FAMILIES` | Families (`hard`, `vp`, `os`, `yield`, `elas`, `visco`, `kin`): prefix, family module, arguments and functions of the contract (`ModelFunc`), finite-difference defaults, composition. |
| `HARDENING_MODELS`, `YIELD_MODELS`, `VISCOPLASTIC_MODELS`, ... | The registered models of each family (and subsets such as `YIELD_ISOTROPIC`, `VISCOPLASTIC_COMPOSED`). |
| `MODEL_PARAMS[(prefix, model)]` | Components of the parameter type (`ModelParam`: name, type, shape, doc, default, private). |
| `MODEL_INITS[(prefix, model)]` | Inputs of the initializer `<prefix>_<model>_init`, when the parameters are derived (e.g. Lamé constants from E, ν) or validated. |
| `MODEL_LABELS` | Readable name of the model for the generated docs. |
| `MATERIAL_SOLVERS` | Solvers usable by a material, the families each one needs and its history size. |

The order of the props of `mat_init_<name>` is generated from `MODEL_PARAMS` / `MODEL_INITS`: adding a
model automatically gives it a parameter table (`python3 tools/fypp/muscle_meta.py props <material>`).

## How to add a model to an existing family

Example: a new hardening law `mylaw` with \( \sigma_y = \sigma_{y0} + K\epsilon_p^2 \).

1. **Register it** in `tools/fypp/muscle_meta.py`: add `'mylaw'` to `HARDENING_MODELS`, its parameters to
   `MODEL_PARAMS` and a label to `MODEL_LABELS`:

   ```python
   ('hard', 'mylaw'): (
       _P('y0', 'real', '', r'Initial yield stress \( \sigma_{y0} \).', None),
       _P('K', 'real', '', r'Hardening coefficient \( K \).', None),
   ),
   ```

2. **Write the template** `src/plasticity/hardening_laws/muscle_hard_mylaw.fypp` with the contract macros
   of `include/muscle_models_common.fypp` (copy `muscle_hard_bilinear.fypp`):

   ```fortran
   #:include "muscle_models_common.fypp"
   #:set PARAMS = model_param_list('hard', 'mylaw')
   module muscle_hard_mylaw
       !! Quadratic hardening law ...
       use muscle_kinds, only : ${', '.join(REAL_KINDS)}$
       implicit none
       private
       $:MODEL_PUBLIC('hard', 'mylaw')
   #:for k in REAL_KINDS
       $:MODEL_PARAMS_TYPE('hard', 'mylaw', k, PARAMS)
   #:endfor
   contains
   #:for k in REAL_KINDS
       $:MODEL_FUNCTION('hard', 'mylaw', 'stress', k)
           res = prm%y0 + prm%K*ep**2
       end function hard_stress_mylaw_${k}$
       ! ... dstress_dep and ddstress_ddep (analytical; MODEL_FD gives finite differences if needed)
   #:endfor
   end module muscle_hard_mylaw
   ```

   Every kernel procedure is `pure` and device-safe (no `class`, `allocatable`, I/O, `stop`, LAPACK; the
   macros add `@:DEVICE()`); literals and constants carry the kind (`1.0_${k}$`, constants of
   `muscle_kinds`).

3. **Add it to the build**: list the template in `muscle_add_fypp_sources` of the folder's
   `CMakeLists.txt` (before the family module). The family module (`muscle_hard_laws`) re-exports the
   new specifics in the generics from the registry, so it needs no edit.

4. **Tests** (`tests/plasticity/hardening_laws/`): values in `sp` and `dp` (`*_kinds.fypp`), analytical
   derivatives against finite differences, and `tests/tools/test_muscle_meta.py` /
   `test_muscle_model_params.py` for the metadata.

5. **Make it reachable from a material**: add a `Material(...)` that uses it to the test registry
   `tools/fypp/muscle_materials.py` (so that CI compiles and runs it through a solver), or to your own
   registry passed with `-DMUSCLE_MATERIALS=<path>`:

   ```python
   Material('cp_vm_mylaw', 'closest_point', elas='linear', yield_='vonmises', hard='mylaw',
            doc='von Mises with quadratic hardening')
   ```

The same steps apply to the other families: a yield criterion writes `stress_eq` and its gradient and
Hessian (and may write the joint evaluation, see `ModelJoint`); a viscoplastic or overstress law that
composes a hardening law receives its values (`sy`, `dsy_dep`) as arguments (`composed`).

## How to add a solver

A solver is a body `src/solvers/muscle_<solver>_body.fypp` included by the material template
`src/materials/muscle_material.fypp` with the kind fixed, plus a `MaterialSolver` entry in
`MATERIAL_SOLVERS` (families it calls, history size, whether it uses `dt`, default iteration limit). It
only calls family generics, so it works with every registered model of those families. Its body declares
the common dummy arguments of `SOLVE` with the macro `SOLVE_DUMMIES` of the material template. Add its
`#:include` to the solver chain of `muscle_material.fypp` (fypp only takes literal file names; the CMake
dependencies come from `MATERIAL_SOLVERS`, and `tests/tools/test_muscle_materials.py` checks that the chain
matches it), and golden tests in `tests/materials/`.

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

## Verification

`tools/check.sh` lints the mechanical rules (kinds, literals, forbidden constructs in kernels, SPDX) and
builds and tests `sp;dp` in Debug. See also the [migration guide](migration.html) for the mapping from the
former object-oriented API.
