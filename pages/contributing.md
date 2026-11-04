title: Adding New Models

# Developer Guide

Models are not added by extending a base class: there are no abstract types nor type-bound model
procedures. A model is a **concrete parameter type** plus **free `pure` functions** with the signature of
its family, written once as an FYPP template for both kinds (`sp`, `dp`) and, for tensor families, for
both dimensions. Solvers call the family generics (`hard_stress`, `yield_stress_eq`, ...), which resolve at
compile time from the type of the parameters, so a new model is usable by every solver of its family
without touching the solvers. Rules: `.claude/rules/` (templates, kinds, GPU-safe code) and the
decisions D-017 / D-018 in `docs/fypp/DECISIONES.md`.

A model is written **by hand, in Fortran, in one file**: the parameter type, the signatures with their `intent`
and documentation, the bodies and the closings are in the template. FYPP only instantiates the text for each
kind (`#:for k in REAL_KINDS`, `${k}$`) and marks the kernels (`@:DEVICE()`), which is what Fortran cannot do
(rule «minimal FYPP», D-020: `.claude/rules/fypp-templates.md`). No model declaration is built from metadata. The
only generated modules are the materials, one per entry of the registry (type, parameter array and APIs of each
material), because Fortran cannot declare one module per entry; they are written once in `src/materials/`.

`tools/fypp/muscle_meta.py` holds **data**, used where a material is generated (the parameter array of
`mat_init_<name>`) and by the tests:

| Registry | Content |
|----------|---------|
| `MODEL_PARAMS[(prefix, model)]` | Names and order of the components of the parameter type, with their doc. The generated `mat_init_<name>` assigns them by name (`mat%yield%sigy = props(3)`), so the compiler rejects a name that the type does not have; `tests/tools/test_muscle_model_params.py` also checks the order, the types and the defaults against the type written in the template. |
| `MODEL_INITS[(prefix, model)]` | Inputs of the initializer `<prefix>_<model>_init`, when the parameters are derived (e.g. Lamé constants from E, ν) or validated, with their domains. |
| `HARDENING_MODELS`, `YIELD_MODELS`, `VISCOPLASTIC_MODELS`, ... | The registered models of each family (the materials look them up here). |
| `MODEL_LABELS` | Readable name of the model for the generated docs. |
| `MATERIAL_SOLVERS` | Solvers usable by a material, the families each one needs and its history size. |

The order of the props of `mat_init_<name>` comes from `MODEL_PARAMS` / `MODEL_INITS`
(`python3 tools/fypp/muscle_meta.py props <material>`).

## How to add a model to an existing family

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

       ! ... hard_dstress_dep_mylaw_${k}$ and hard_ddstress_ddep_mylaw_${k}$, written the same way
       ! (analytical; for a derivative by finite differences use the expression macros of
       ! tests/common/muscle_fd.fypp as the body: `@:FD_DERIV(res, hard_stress_mylaw_${k}$(prm, {x}), ep, ${k}$)`)
   #:endfor
   end module muscle_hard_mylaw
   ```

   Every kernel procedure is `pure` and device-safe (no `class`, `allocatable`, I/O, `stop`, LAPACK) and
   carries `@:DEVICE()` right after `implicit none`; literals and constants carry the kind (`1.0_${k}$`,
   constants of `muscle_kinds`).

2. **Add it to the family module** `muscle_hard_laws.fypp`: one `use` line and one name in each of the three
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
   `CMakeLists.txt` (before the family module).

5. **Tests** (`tests/plasticity/hardening_laws/`): values in `sp` and `dp` (`*_kinds.fypp`), analytical
   derivatives against finite differences, and `tests/tools/test_muscle_model_params.py` for the metadata.

6. **Make it reachable from a material**: add a `Material(...)` that uses it to the test registry
   `tools/fypp/muscle_materials.py` (so that CI compiles and runs it through a solver), or to your own
   registry passed with `-DMUSCLE_MATERIALS=<path>`:

   ```python
   Material('cp_vm_mylaw', 'closest_point', elas='linear', yield_='vonmises', hard='mylaw',
            doc='von Mises with quadratic hardening')
   ```

The same steps apply to the other families (copy the neighbour model of the family): a yield criterion
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
