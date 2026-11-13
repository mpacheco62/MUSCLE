# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Python helpers and tensor layout metadata for the MUSCLE FYPP templates.

fypp loads this module with ``-M tools/fypp -m muscle_meta``; ``include/muscle_common.fypp``
re-exports the helpers under short names. Conventions: ``pages/conventions.md``.

Contents (in file order; search for the section title)
------------------------------------------------------
1. Kinds, dimensions and model lists (``KINDS``, ``DIMS``, ``FAMILIES``, ``IDENTITIES``,
   ``*_MODELS``, ``MAX_VISCO_BRANCHES``, ``KINEMATICS_*``, ``MODEL_LABELS``):
   data read by every template (through ``muscle_common.fypp``), the materials and the tests.
2. Tensor layouts and term helpers (``Layout``, ``layout``, ``tname``, ``iname``, ``iden_*``,
   ``join_terms``, ``norm_terms``): ``src/tensors/``, ``include/muscle_tensor_common.fypp``; see
   «Tensor layouts» and «Identities» below.
3. «Operators between tensor types» (``OP_MODULES``, ``OP_TABLE``, ``op_*``):
   ``src/tensors/operators/``, ``include/muscle_tensor_ops_common.fypp``,
   ``tools/fypp/ops_signatures.py`` and ``tests/tensors/operators/``.
4. «Constitutive model families» (``ModelFamily``, ``ModelFunc``, ``MODEL_FAMILIES``,
   ``model_family``, ``model_type``, ``model_procs``, ...): names of the family modules, model types and
   procedures, which functions each model has, the composed values (``sy``, ``dsy_dep``) and the joint
   evaluation, for the material template and the solver bodies. The function signatures are written in
   Fortran; the generated target tests keep their own copy (``tests/gpu/target_cases.py``).
5. «Parameters and initializers of the registered models» (``MODEL_PARAMS``, ``MODEL_INITS``,
   ``model_param_list``, ``model_inputs``): the parameter array of ``mat_init_<name>``; checked against the
   hand-written types by ``tests/tools/test_muscle_model_params.py``.
6. «Declared materials» (``Material``, ``MATERIAL_SOLVERS``, ``load_materials``, ``material_*``): the
   material template (``src/materials/``), ``cmake/MuscleFypp.cmake`` and ``tests/materials/``.
7. Command line (``_main``): ``python3 tools/fypp/muscle_meta.py materials [REGISTRY]`` (names of the
   declared materials), ``solver-bodies`` (the solver body files) and ``props NAME [REGISTRY]`` (Markdown
   table of the parameter array of a material).

Tensor layouts
--------------
``layout(family, dim)`` describes how a tensor type stores its components in ``vals``. The
layouts reproduce the storage of the original ``src/tensors/types/*.f90`` exactly (D-004);
``tests/tools/test_muscle_meta.py`` and the Fortran test generated from this module check them.

Families: ``'2O'``, ``'2Osym'``, ``'4O'``, ``'4O2sym'``, ``'4O3sym'``; dims: ``'3D'``, ``'2D'``.

- 2nd order: ``vals(n)``. Unsymmetric tensors are stored column-major; symmetric ones in
  tensorial Voigt order.
- ``4O`` / ``4O2sym``: ``vals(n, n)`` with ``C_ijkl`` at ``(I(i,j), I(k,l))``, where ``I`` is the
  storage index of the matching 2nd-order layout (``2O`` / ``2Osym``).
- ``4O3sym`` (major and minor symmetry): ``vals(m)`` with the upper triangle of the ``2Osym``
  Voigt matrix stored diagonal by diagonal: (1,1), (2,2), ..., then (1,2), (2,3), ..., then
  (1,3), ... and finally (1,n).

Indices are 1-based, as in Fortran.

Identities
----------
``IDENTITIES`` lists the identity tensors (``iden_2O``, ``iden_4O3T``, ``iden_4O4T``): they do
not depend on the dimension. The unscaled type has no kind and no components; the scaled one
(``iden_<name>S_<k>``) stores the factor ``val``. ``iname`` builds their type names and
``iden_matrix`` the matrix shown by their formatted output.

Operators between types
-----------------------
``OP_TABLE`` lists, for each generic of ``src/tensors/operators/`` (``OP_MODULES``: one module
per operation, D-015), the operand pairs and result by family (``OpPair``), without dimension
or kind. ``op_entries(generic, dim, kind)`` expands it to the specific procedures of one
dimension and kind (``OpEntry``: concrete type names and the historical procedure name with the
kind appended); ``op_interface`` gives the Fortran generic spec. The dp expansion is exactly the
historical set of operators (``tools/fypp/ref/ops_signatures_dp.txt``).

Queries for the templates: ``Layout.index(i, j[, k, l])`` (storage position of any tensor
component, applying the symmetries of the family; ``None`` for a structural zero in 2D),
``Comp.mult`` (how many tensor components the stored value represents, e.g. 2 for a symmetric
shear) and ``Comp.voigt`` / ``Layout.base`` (the 2nd-order layout indexing a 4th-order one).
"""

import importlib.util
import itertools
import os
import re
import sys
from collections import namedtuple
from types import MappingProxyType

# ---------------------------------------------------------------------------------------------
# Kinds, dimensions and model lists
# ---------------------------------------------------------------------------------------------

#: Real kinds known to MUSCLE (Fortran constants ``sp`` and ``dp``).
KINDS = ('sp', 'dp')

#: Space dimensions.
DIMS = ('3D', '2D')

#: Tensor families.
FAMILIES = ('2O', '2Osym', '4O', '4O2sym', '4O3sym')

#: Identity tensors (unscaled names; the scaled ones append ``S``): ``2O`` is \delta_ij,
#: ``4O3T`` is \delta_ij \delta_kl and ``4O4T`` is (\delta_ik \delta_jl + \delta_il \delta_jk)/2.
IDENTITIES = ('2O', '4O3T', '4O4T')

#: Registered isotropic hardening laws (D-007.2): one template
#: ``src/plasticity/hardening_laws/muscle_hard_<name>.fypp`` each (D-017), exported by the family module ``muscle_hard_laws``.
HARDENING_MODELS = ('bilinear', 'ludwik', 'swift', 'vocemod')

#: Registered viscoplastic flow stress laws \( \sigma(\epsilon_p, \dot\epsilon_p) \) (D-017): one
#: template ``src/plasticity/viscoplastic_laws/muscle_vp_<name>.fypp`` each, exported by
#: ``muscle_vp_laws``. ``VISCOPLASTIC_COMPOSED``: those that compose a hardening law (D-017.6b).
VISCOPLASTIC_MODELS = ('jc', 'jcfull', 'mrk', 'nnl', 'rk', 'va')
VISCOPLASTIC_COMPOSED = ('jc', 'jcfull')
#: Viscoplastic laws with a rate floor (the flow stress is constant below a minimum rate) that give the
#: Newton iterations of the rate-dependent return mapping the right-hand rate derivative there
#: (``vp_newton_dstress_dep``, D-018.19b); the others use ``vp_dstress_dep``.
VISCOPLASTIC_FLOOR = ('jc',)

#: Registered overstress laws \( \dot\lambda = g(\sigma_{eq}, \epsilon_p) \) (D-017): one template
#: ``src/plasticity/viscoplastic_laws/muscle_os_<name>.fypp`` each (the directory of the flow stress laws),
#: exported by ``muscle_os_laws``.
OVERSTRESS_MODELS = ('chaboche', 'perzyna')
OVERSTRESS_COMPOSED = ('chaboche', 'perzyna')

#: Registered yield criteria \( \sigma_{eq}(\boldsymbol\sigma) \) (D-017): one template
#: ``src/plasticity/yield_criteria/muscle_yield_<name>.fypp`` each, exported by
#: ``muscle_yield_laws``. ``YIELD_JOINT``: those that write their own joint evaluation of the
#: equivalent stress and its derivatives (the others get the default, three calls in turn).
YIELD_MODELS = ('vonmises', 'hill48', 'druckerprager', 'cpb06', 'cpb06ex2', 'yld2004')
YIELD_JOINT = ('yld2004',)
#: Isotropic yield criteria (the others have anisotropy axes): with a spatial strain measure, a material
#: with another criterion runs its return mapping in the frame of the polar rotation (D-018.15f).
YIELD_ISOTROPIC = ('vonmises', 'druckerprager')

#: Registered elastic models (D-017.12): one template
#: ``src/elasticity/muscle_elas_<name>.fypp`` each (the frozen elasticity:
#: ``src/viscoelasticity/muscle_elas_frozen.fypp``), exported by ``muscle_elas_laws``.
ELASTICITY_MODELS = ('linear', 'frozen')
#: Isotropic elastic models: a return mapping in a rotated frame (spatial strain measure with an
#: anisotropic yield criterion, D-018.15f) needs one (the elasticity is then frame independent).
ELASTICITY_ISOTROPIC = ('linear', 'frozen')

#: Registered linear viscoelastic models (D-017.13): one template
#: ``src/viscoelasticity/muscle_visco_<name>.fypp`` each, exported by ``muscle_visco_laws``. Their
#: step map is affine in the elastic strain: the algorithmic ("frozen") elasticity of the step is
#: the elastic model ``frozen`` (``elas_frozen_<k>``), which a return mapping uses as any other.
VISCOELASTIC_MODELS = ('maxwell', 'kelvin')

#: Maximum number of branches (history tensors) of a viscoelastic state (``MAX_VISCO_BRANCHES`` of
#: muscle_visco_types, D-017.13c): the size of the fixed state and the largest
#: ``Material.visco_branches``.
MAX_VISCO_BRANCHES = 10

#: Registered kinematics (D-017.14): one template ``src/kinematics/muscle_kin_<name>.fypp`` each,
#: exported by ``muscle_kin_laws``. They have no material parameters: the parameter type of a
#: kinematics is empty and only selects it at compile time (generic resolution, D-008 renames).
KINEMATICS_MODELS = ('small_strain', 'total_lagrangian', 'updated_lagrangian', 'material_logarithmic',
                     'spatial_logarithmic', 'corotational')
#: Kinematics updated from F (``kin_update``); the others (the corotational one) from a strain and a
#: rotation (``kin_corot_update``, D-017.14d).
KINEMATICS_F_MODELS = tuple(m for m in KINEMATICS_MODELS if m != 'corotational')
#: Kinematics whose material kernel takes the strain already rotated to the corotational frame by the
#: FE code, as Abaqus and ANSYS do (D-018.15f): the small strain signature, Jaumann tangent; the
#: rotation of the stored plastic strain is the job of the FE adapter.
KINEMATICS_FE_ROTATED = ('corotational',)
#: Kinematics with a spatial strain measure (Euler-Almansi, spatial Hencky): a material rotates its
#: stored plastic strain with the rotation increment of the polar decomposition (D-018.15f).
KINEMATICS_SPATIAL = ('updated_lagrangian', 'spatial_logarithmic')

#: Kinematics whose own family functions read the polar rotation R of the state (``kin_to_cauchy`` and
#: the tangents of the material logarithmic kinematics): their ``kin_update_norot`` is ``kin_update`` (D-019.7c).
KINEMATICS_ROTATION = ('material_logarithmic',)


#: Readable names of the models in the generated FORD docs (``model_label``); a model not listed
#: (e.g. a test-only law) is named by its identifier.
MODEL_LABELS = {('hard', 'bilinear'): 'bilinear', ('hard', 'ludwik'): 'Ludwik', ('hard', 'swift'): 'Swift',
                ('hard', 'vocemod'): 'modified Voce', ('vp', 'jc'): 'Johnson-Cook',
                ('vp', 'rk'): 'Rusinek-Klepaczko', ('vp', 'jcfull'): 'full Johnson-Cook',
                ('vp', 'mrk'): 'modified Rusinek-Klepaczko', ('vp', 'nnl'): 'Nemat-Nasser-Li',
                ('vp', 'va'): 'VA', ('os', 'perzyna'): 'Perzyna',
                ('os', 'chaboche'): 'Chaboche', ('yield', 'vonmises'): 'von Mises',
                ('yield', 'hill48'): 'Hill (1948)', ('yield', 'cpb06'): 'CPB06 (Cazacu-Plunkett-Barlat 2006)',
                ('yield', 'druckerprager'): 'extended linear Drucker-Prager',
                ('yield', 'cpb06ex2'): 'CPB06ex2 (Plunkett-Cazacu-Barlat 2008)',
                ('yield', 'yld2004'): 'Yld2004-18p (Barlat et al. 2005)',
                ('elas', 'linear'): 'linear isotropic elastic',
                ('elas', 'frozen'): 'frozen (algorithmic) viscoelastic',
                ('visco', 'maxwell'): 'generalized Maxwell', ('visco', 'kelvin'): 'Kelvin-Voigt',
                ('kin', 'small_strain'): 'small strain', ('kin', 'total_lagrangian'): 'total Lagrangian',
                ('kin', 'material_logarithmic'): 'material logarithmic',
                ('kin', 'updated_lagrangian'): 'updated Lagrangian',
                ('kin', 'spatial_logarithmic'): 'spatial logarithmic', ('kin', 'corotational'): 'corotational'}

# ---------------------------------------------------------------------------------------------
# Tensor layouts and term helpers
# ---------------------------------------------------------------------------------------------

_AXES = 'xyz'


def join_terms(terms, op='+', indent=12, lead=None):
    """Joins Fortran terms with ``op``, one term per line (free-form continuation), so that long
    expressions stay below the 132-column limit: ``join_terms(['a', 'b'])`` ->
    ``'a &\n            + b'``. ``lead``: columns used on the line before the expression; if
    given, an expression that fits in the line is kept on it (``join_terms(['a', 'b'], lead=8)``
    -> ``'a + b'``)."""
    line = ' {} '.format(op).join(terms)
    if lead is not None and lead + len(line) <= 132:
        return line
    return ' &\n{}{} '.format(' ' * indent, op).join(terms)


def assign_terms(lhs, terms, col, op='+'):
    """Fortran assignment ``lhs = t1 + t2 + ...`` written at column ``col`` (the indentation of
    the line in the template): on one line if it fits in 132 columns, else one term per line with
    the operators aligned under the first term (``join_terms``):
    ``assign_terms('W11', ['a', 'b'], 8)`` -> ``'W11 = a + b'``."""
    start = col + len(lhs) + len(' = ')
    return '{} = {}'.format(lhs, join_terms(terms, op=op, indent=start, lead=start))


def family_module(family):
    """Fortran module of a tensor family (``FAMILIES``) or of the identities (``'iden'``):
    ``family_module('2Osym')`` -> ``'muscle_tensor_2osym'``."""
    if family not in FAMILIES and family != 'iden':
        raise ValueError(f"unknown tensor family '{family}', expected one of {FAMILIES} or 'iden'")
    return 'muscle_tensor_' + family.lower()


def norm_terms(lay, var):
    """Terms of the family ``norm`` (modified L1 norm) of a 2nd-order tensor whose ``vals`` is
    ``var``, in storage order, a shear component of ``2Osym`` weighted by its multiplicity 2:
    ``norm_terms(layout('2Osym', '2D'), 'a%vals')`` -> ``['abs(a%vals(1))', ..., '2*abs(a%vals(4))']``.
    Shared by the ``norm`` functions of the families and the finite-difference macros, so that
    both give the same bits."""
    if lay.base is not None:
        raise ValueError(f"norm_terms: '{lay.dim}{lay.family}' is not a 2nd-order layout")
    return [('{}*'.format(c.mult) if c.mult > 1 else '') + 'abs({}({}))'.format(var, c.pos[0])
            for c in lay.comps]


def tname(base, kind):
    """Name of the concrete tensor type: ``tname('3D2Osym', 'dp')`` -> ``'ten_3D2Osym_dp'``."""
    if kind not in KINDS:
        raise ValueError(f"unknown kind '{kind}', expected one of {KINDS}")
    return f"ten_{base}_{kind}"


def iname(base, kind=None):
    """Name of an identity type: ``iname('2O')`` -> ``'iden_2O'`` (unscaled, no kind);
    ``iname('2OS', 'dp')`` -> ``'iden_2OS_dp'`` (scaled, one type per kind)."""
    if kind is None:
        return f"iden_{base}"
    if kind not in KINDS:
        raise ValueError(f"unknown kind '{kind}', expected one of {KINDS}")
    return f"iden_{base}_{kind}"


def iden_coef(name, idx):
    """Tensor component ``idx`` (``(i, j)`` or ``(i, j, k, l)``) of the unscaled identity
    ``iden_<name>``: 0, 1 or 0.5. ``iden_coef('4O4T', (1, 2, 1, 2))`` -> ``0.5``."""
    def delta(a, b):
        return 1 if a == b else 0

    if name not in IDENTITIES:
        raise ValueError(f"unknown identity '{name}', expected one of {IDENTITIES}")
    if len(idx) != (2 if name == '2O' else 4):
        raise ValueError(f"identity '{name}' has no component {idx}")
    if name == '2O':
        return delta(*idx)
    i, j, k, l = idx
    if name == '4O3T':
        return delta(i, j) * delta(k, l)
    return (delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k)) / 2


def iden_combination(lay, terms, kind):
    """Stored components of the layout ``lay`` where a sum of identities is not zero, in storage
    order, as ``[(subscript, expr)]``: ``subscript`` is the Fortran subscript of ``vals``
    (``'4'`` or ``'1, 1'``) and ``expr`` the value of the component.

    ``terms`` is a list of ``(name, var, sign)``: identity ``iden_<name>`` (``IDENTITIES``) scaled
    by the Fortran expression ``var`` (``None`` for the unscaled one, factor 1), with ``sign``
    ``'+'`` or ``'-'``. Each term contributes ``coef * var`` (``coef`` from ``iden_coef``), written
    as ``var`` or ``0.5_<kind> * var``; the terms keep their order and are joined with their
    signs, and unscaled terms are literals with ``kind`` (folded into one literal if every term is
    unscaled): ``iden_combination(layout('4O3sym', '3D'), [('4O3T', None, '+'),
    ('4O4T', 'b%val', '-')], 'dp')[0]`` -> ``('1', '1.0_dp - b%val')``."""
    def lit(x):
        return f'{float(x)!r}_{kind}'

    if kind not in KINDS:
        raise ValueError(f"unknown kind '{kind}', expected one of {KINDS}")
    out = []
    for c in lay.comps:
        parts = []   # (sign, number) for unscaled terms, (sign, text) for scaled ones
        for name, var, sign in terms:
            if sign not in '+-' or len(sign) != 1:
                raise ValueError(f"sign '{sign}' must be '+' or '-'")
            coef = iden_coef(name, c.idx)
            if coef == 0:
                continue
            if var is None:
                parts.append((sign, coef))
            else:
                parts.append((sign, var if coef == 1 else f'{lit(coef)} * {var}'))
        if all(isinstance(t, (int, float)) for _, t in parts):
            total = sum(t if s == '+' else -t for s, t in parts)
            if total:
                out.append((', '.join(map(str, c.pos)), lit(total)))
            continue
        expr = ''
        for n, (sign, t) in enumerate(parts):
            text = lit(t) if isinstance(t, (int, float)) else t
            if n == 0:
                expr = text if sign == '+' else f'-{text}'
            else:
                expr += f' {sign} {text}'
        out.append((', '.join(map(str, c.pos)), expr))
    return out


def iden_matrix(name, one, half, zero):
    """Rows of the matrix of the identity ``iden_<name>`` (``IDENTITIES``) shown by its formatted
    output, as Fortran expressions: ``one`` for the entries equal to the factor, ``half`` for half
    of it and ``zero`` for zeros. ``2O``: the 3x3 matrix \\delta_ij. ``4O3T`` / ``4O4T``: the 6x6
    Voigt matrix C_IJ (3D ``2Osym`` order 11, 22, 33, 12, 23, 13; tensorial components, no shear
    factors) of \\delta_ij \\delta_kl (ones in the normal 3x3 block) and of
    (\\delta_ik \\delta_jl + \\delta_il \\delta_jk)/2 (``diag(1, 1, 1, 1/2, 1/2, 1/2)``)."""
    if name == '2O':
        coef = [[iden_coef(name, (i, j)) for j in range(1, 4)] for i in range(1, 4)]
    else:
        pairs = [c.idx for c in layout('2Osym', '3D').comps]
        coef = [[iden_coef(name, r + c) for c in pairs] for r in pairs]
    expr = {0: zero, 1: one, 0.5: half}
    return [[expr[v] for v in row] for row in coef]


Comp = namedtuple('Comp', 'name pos idx mult voigt')
Comp.__doc__ = """One stored component.

name:  component name, e.g. ``'xy'`` or ``'xxyz'`` (the keyword of ``init2``).
pos:   storage position in ``vals``, ``(I,)`` for vector storage or ``(I, J)`` for matrix storage.
idx:   tensor indices ``(i, j)`` or ``(i, j, k, l)`` of the representative component.
mult:  number of tensor components equal to this one by symmetry (1 if unsymmetric).
voigt: 4th order: ``(I, J)`` in the matrix of the base layout (equal to ``pos`` for matrix
       storage); 2nd order: ``pos``.
"""


class Layout(namedtuple('Layout', 'family dim shape comps base')):
    """Storage of one tensor type: ``vals`` has shape ``shape``; ``comps`` is in storage order
    (column-major for matrix storage). ``base`` is the 2nd-order layout that indexes the rows
    and columns of a 4th-order family (``None`` for 2nd order)."""

    __slots__ = ()

    def index(self, *idx):
        """Storage position of the tensor component ``idx`` or ``None`` if it is not stored
        (structural zero of a 2D type)."""
        return _index_of(self)[_canonical(self.family, idx)]

    def ref(self, var, *idx, zero='0'):
        """Fortran reference to the tensor component ``idx`` stored in the array ``var``, or
        ``zero`` (a symbolic zero) if the component is not stored:
        ``layout('2Osym', '2D').ref('a%vals', 1, 3, zero='0.0_dp')`` -> ``'0.0_dp'``."""
        pos = self.index(*idx)
        if pos is None:
            return zero
        return '{}({})'.format(var, ', '.join(str(p) for p in pos))

    def matrix(self, var, zero='0'):
        """Rows of the matrix shown to the user, as Fortran references (``ref``): the full 3x3
        matrix for 2nd order; for 4th order the matrix of the base layout (``C(I, J)`` with I, J
        the rows/columns of the 2nd-order storage, e.g. 6x6 Voigt for ``4O3sym`` in 3D)."""
        if self.base is None:
            return [[self.ref(var, i, j, zero=zero) for j in range(1, 4)] for i in range(1, 4)]
        rows = [c.idx for c in self.base.comps]
        return [[self.ref(var, *r, *c, zero=zero) for c in rows] for r in rows]


def _name(*idx):
    return ''.join(_AXES[i - 1] for i in idx)


def _canonical(family, idx):
    """Representative of the symmetry class of ``idx`` within ``family``."""
    pairs = [tuple(idx[n:n + 2]) for n in range(0, len(idx), 2)]
    if family in ('2Osym', '4O2sym', '4O3sym'):
        pairs = [tuple(sorted(p)) for p in pairs]
    if family == '4O3sym':
        pairs = sorted(pairs)
    return tuple(i for p in pairs for i in p)


def _orbit(family, idx):
    """All tensor indices equal to ``idx`` by the symmetries of ``family``."""
    order = len(idx)
    return {t for t in _all_indices(order) if _canonical(family, t) == _canonical(family, idx)}


def _all_indices(order):
    return itertools.product(range(1, 4), repeat=order)


_INDEX_CACHE = {}


def _index_of(lay):
    """Map canonical tensor indices -> storage position (``None`` for unstored ones)."""
    key = (lay.family, lay.dim)
    if key not in _INDEX_CACHE:
        table = {_canonical(lay.family, t): None for t in _all_indices(len(lay.comps[0].idx))}
        for c in lay.comps:
            table[_canonical(lay.family, c.idx)] = c.pos
        _INDEX_CACHE[key] = table
    return _INDEX_CACHE[key]


# Representative (i, j) of each stored 2nd-order component, copied from the .f90 sources.
_SECOND_ORDER = {
    ('2O', '3D'): [(1, 1), (2, 1), (3, 1), (1, 2), (2, 2), (3, 2), (1, 3), (2, 3), (3, 3)],
    ('2O', '2D'): [(1, 1), (2, 1), (1, 2), (2, 2), (3, 3)],
    ('2Osym', '3D'): [(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (1, 3)],
    ('2Osym', '2D'): [(1, 1), (2, 2), (3, 3), (1, 2)],
}

# 2nd-order layout whose storage index addresses the rows/columns of each 4th-order family
_BASE_OF = {'4O': '2O', '4O2sym': '2Osym', '4O3sym': '2Osym'}


def _build(family, dim):
    def comp(idx, pos, voigt):
        return Comp(_name(*idx), pos, idx, len(_orbit(family, idx)), voigt)

    if (family, dim) in _SECOND_ORDER:
        ij = _SECOND_ORDER[(family, dim)]
        comps = [comp(p, (s,), (s,)) for s, p in enumerate(ij, start=1)]
        return Layout(family, dim, (len(ij),), tuple(comps), None)

    base_layout = _build(_BASE_OF[family], dim)
    base = _SECOND_ORDER[(_BASE_OF[family], dim)]
    n = len(base)
    if family == '4O3sym':
        pairs = [(row, row + off) for off in range(n) for row in range(1, n - off + 1)]
        comps = [comp(base[r - 1] + base[c - 1], (s,), (r, c))
                 for s, (r, c) in enumerate(pairs, start=1)]
        return Layout(family, dim, (len(pairs),), tuple(comps), base_layout)

    comps = [comp(base[r - 1] + base[c - 1], (r, c), (r, c))
             for c in range(1, n + 1) for r in range(1, n + 1)]
    return Layout(family, dim, (n, n), tuple(comps), base_layout)


_LAYOUTS = {(f, d): _build(f, d) for f in FAMILIES for d in DIMS}


def layout(family, dim):
    """Storage layout of the tensor type ``ten_<dim><family>``, e.g. ``layout('2Osym', '2D')``."""
    try:
        return _LAYOUTS[(family, dim)]
    except KeyError:
        raise ValueError(f"unknown layout '{dim}{family}', families {FAMILIES}, dims {DIMS}")


# ---------------------------------------------------------------------------------------------
# Operators between tensor types (Phase 4, D-015)
# ---------------------------------------------------------------------------------------------

OpPair = namedtuple('OpPair', 'a b res mirror prefix', defaults=(False, None))
OpPair.__doc__ = """One row of the operator table ``OP_TABLE``: ``a <op> b -> res``.

a, b:   operands. A tensor family (``'2Osym'``: the dimension is that of the expansion), a
        tensor with an explicit dimension (``'3D4O3sym'``) or an identity tag (``'I2O'``,
        ``'I2OS'``, ``'I4O3T'``, ...: ``'I'`` + unscaled name, ``'S'`` appended if scaled).
        ``b`` is ``None`` for a unary operator. For ``=``, ``a`` is the left-hand side.
res:    result, written like an operand, or ``'real'`` (scalar of the kind of the expansion);
        ``None`` for ``=``.
mirror: the swapped row ``b <op> a -> res`` exists too (structural only: it says nothing about
        whether both give the same value).
prefix: prefix of the specific procedure name if it differs from that of the generic
        (``OP_PREFIX``). The specific names are private to the operator modules (only the
        generics are public): the prefixes only keep the historical names.
"""

OpEntry = namedtuple('OpEntry', 'generic name a b res swapped tag_a tag_b tag_res dim kind')
OpEntry.__doc__ = """One specific procedure of a generic, for one dimension and kind
(``op_entries``).

generic: key of ``OP_GENERICS`` (``'+'``, ``'.ddot.'``, ``'='``, ...).
name:    specific procedure name, historical name + ``_<kind>`` (``sum_I2O_3D2Osym_dp``).
a, b:    concrete type names of the arguments (``ten_3D2Osym_dp``, ``iden_2O``,
         ``iden_2OS_sp``); ``b`` is ``None`` for a unary operator.
res:     concrete result type (``ten_3D2Osym_dp``, ``real(dp)``), ``None`` for ``=``.
swapped: the entry is the mirror (``b <op> a``) of a row of the table.
tag_a, tag_b, tag_res: the operands and result as written in the table, with the dimension
         resolved (``'3D2Osym'``, ``'I2OS'``, ``'real'``; ``None`` where absent), in the order of
         ``a``, ``b``. The template picks the body from them (``op_iden`` tells an identity).
dim, kind: dimension and kind of the expansion (``'3D'``, ``'dp'``) that produced the entry.
"""

#: Module of the operators of each generic (``src/tensors/operators/muscle_tensor_ops_<m>``).
OP_MODULES = {
    'addition_subtraction': ('+', '-'),
    'assignment': ('=',),
    'contraction_double': ('.ddot.',),
    'contraction_single': ('*',),
    'dyadic': ('.tdot.', '.tdotsym.'),
    'transform': ('.transform.',),
}

#: Generics in module order.
OP_GENERICS = tuple(g for gens in OP_MODULES.values() for g in gens)

#: Prefix of the specific procedure names of each generic (historical names).
OP_PREFIX = {'+': 'sum', '-': 'sub', '=': 'assign', '.ddot.': 'ddot', '*': 'dot',
             '.tdot.': 'tdot', '.tdotsym.': 'tdotsym', '.transform.': 'transform'}

# Sum and difference: identity +- tensor of the matching order, both orders; 2O +- 2Osym. The
# sums of two 4th-order identities (one of each kind, 4O3T(S) with 4O4T(S); never two of the same
# identity) have no tensor operand and give a 3D tensor in the original: they exist only in 3D.
_ADD_SUB = (
    [OpPair(i, f, f, True) for i in ('I2O', 'I2OS') for f in ('2Osym', '2O')]
    + [OpPair(i, f, f, True) for i in ('I4O3TS', 'I4O4TS', 'I4O3T', 'I4O4T')
       for f in ('4O2sym', '4O3sym')]
    + [OpPair(a, b, '3D4O3sym', True) for a, b in (('I4O3TS', 'I4O4T'), ('I4O3TS', 'I4O4TS'),
                                                   ('I4O3T', 'I4O4T'), ('I4O3T', 'I4O4TS'))]
    + [OpPair('2O', '2Osym', '2O', True)]
)

#: Operator table: generic -> rows (``OpPair``). Reproduces exactly the historical double
#: precision operators (``tools/fypp/ref/ops_signatures_dp.txt``, checked by
#: ``tests/tools/test_muscle_meta.py``), gaps included: e.g. ``iden_4O4T .ddot. 2Osym`` does not
#: exist (only the scaled ``iden_4O4TS``), ``.tdot.``/``.tdotsym.`` take only ``iden_2O``, the
#: 4th-order assignments of identities only some of the unscaled/scaled ones.
OP_TABLE = {
    '+': _ADD_SUB,
    '-': _ADD_SUB,
    '=': [
        OpPair('2O', '2Osym', None),
        OpPair('2Osym', 'I2O', None),
        OpPair('2Osym', 'I2OS', None),
        OpPair('4O2sym', '4O3sym', None),
        OpPair('4O2sym', 'I4O4T', None),
        OpPair('4O2sym', 'I4O4TS', None),
        OpPair('4O3sym', 'I4O3TS', None),
        OpPair('4O3sym', 'I4O4TS', None),
    ],
    '.ddot.': [
        OpPair('I2O', '2Osym', 'real', True),
        OpPair('I2O', '2O', 'real', True),
        OpPair('I2OS', '2Osym', 'real', True),
        OpPair('I2OS', '2O', 'real', True),
        OpPair('2O', '2Osym', 'real', True),
        OpPair('I4O4T', '2O', '2Osym', True),
        OpPair('I4O4TS', '2O', '2Osym', True),
        OpPair('I4O4TS', '2Osym', '2Osym', True),
        OpPair('4O2sym', '2Osym', '2Osym', True),
        OpPair('4O3sym', '2Osym', '2Osym', True),
        OpPair('4O3sym', '4O3sym', '4O2sym'),
    ],
    # Single contraction: 'dot_' between tensors, 'mul_' with an identity (historical names).
    '*': [
        OpPair('2Osym', '2Osym', '2O'),
        OpPair('2O', '2Osym', '2O', True),
        OpPair('2O', '2O', '2O'),
        OpPair('I2O', '2Osym', '2Osym', True, 'mul'),
        OpPair('I2O', '2O', '2O', True, 'mul'),
        OpPair('I2OS', '2Osym', '2Osym', True, 'mul'),
        OpPair('I2OS', '2O', '2O', True, 'mul'),
    ],
    '.tdot.': [
        OpPair('2Osym', '2Osym', '4O2sym'),
        OpPair('2Osym', 'I2O', '4O2sym', True),
    ],
    '.tdotsym.': [
        OpPair('2Osym', '2Osym', '4O3sym'),
        OpPair('2Osym', None, '4O3sym'),
        OpPair('2Osym', 'I2O', '4O3sym', True),
    ],
    '.transform.': [
        OpPair('2O', '2Osym', '2Osym'),
        OpPair('2O', '4O2sym', '4O2sym'),
        OpPair('2O', '4O3sym', '4O3sym'),
    ],
}

#: Kind of the result of an operator whose operands carry no kind (only unscaled identities,
#: e.g. ``iden_4O3T + iden_4O4T``): double precision, as for ``I + I`` (D-014, point 3). Those
#: operators exist only in this kind (one per argument types, or the generic is ambiguous).
#: The same rule is written as ``_dp`` literals in ``src/tensors/types/muscle_tensor_iden.fypp``.
OP_KINDLESS_KIND = 'dp'


def op_interface(generic):
    """Fortran generic spec of ``generic``: ``'operator(+)'``, ``'assignment(=)'``."""
    if generic not in OP_TABLE:
        raise ValueError(f"unknown generic '{generic}', expected one of {OP_GENERICS}")
    return 'assignment(=)' if generic == '=' else f'operator({generic})'


def _op_dim(operand):
    """Explicit dimension of a table operand (``'3D4O3sym'`` -> ``'3D'``), else ``None``."""
    if operand and operand[:2] in DIMS:
        return operand[:2]
    return None


def op_iden(operand):
    """``(name, scaled)`` of an identity tag of the table (``'I4O4TS'`` -> ``('4O4T', True)``),
    ``None`` for anything else (tensor tag, ``'real'``, ``None``)."""
    if not operand or not operand.startswith('I'):
        return None
    if operand[1:] in IDENTITIES:
        return operand[1:], False
    if operand.endswith('S') and operand[1:-1] in IDENTITIES:
        return operand[1:-1], True
    raise ValueError(f"unknown identity '{operand}'")


def _op_has_kind(operand):
    """Whether the concrete type of a table operand depends on the kind."""
    iden = op_iden(operand)
    return operand is not None and (iden is None or iden[1])


def _op_concrete(operand, dim, kind):
    """``(tag, type)`` of a table operand: the tag used in the procedure name and the type name."""
    if operand is None:
        return None, None
    if operand == 'real':
        return None, f'real({kind})'
    iden = op_iden(operand)
    if iden is not None:
        name, scaled = iden
        return operand, iname(name + 'S', kind) if scaled else iname(name)
    base = operand if _op_dim(operand) else dim + operand
    layout(base[2:], base[:2])  # validates the family
    return base, tname(base, kind)


def _op_row_dims(row):
    """Dimensions in which a row exists: the explicit one of its operands/result, else all."""
    dims = {_op_dim(x) for x in (row.a, row.b, row.res)} - {None}
    if len(dims) > 1:
        raise ValueError(f'row {row} mixes dimensions {sorted(dims)}')
    return tuple(dims) or DIMS


def op_entries(generic, dim, kind):
    """Specific procedures (``OpEntry``) of ``generic`` for one dimension and kind, in table order
    (each row followed by its mirror). Rules (D-014, D-015):

    - tensors and scaled identities take ``kind`` (no mixed kinds); unscaled identities have no
      kind and combine with tensors of every kind;
    - a row whose operands carry no kind exists only in ``OP_KINDLESS_KIND``;
    - a row with an explicit dimension exists only in that dimension (sums of identities).

    ``op_entries('+', '3D', 'dp')[0]`` -> ``OpEntry('+', 'sum_I2O_3D2Osym_dp', 'iden_2O',
    'ten_3D2Osym_dp', 'ten_3D2Osym_dp', False, 'I2O', '3D2Osym', '3D2Osym', '3D', 'dp')``."""
    if generic not in OP_TABLE:
        raise ValueError(f"unknown generic '{generic}', expected one of {OP_GENERICS}")
    if dim not in DIMS:
        raise ValueError(f"unknown dimension '{dim}', expected one of {DIMS}")
    if kind not in KINDS:
        raise ValueError(f"unknown kind '{kind}', expected one of {KINDS}")
    out = []
    for row in OP_TABLE[generic]:
        if dim not in _op_row_dims(row):
            continue
        if not (_op_has_kind(row.a) or _op_has_kind(row.b)) and kind != OP_KINDLESS_KIND:
            continue
        tag_res, res = _op_concrete(row.res, dim, kind)
        tag_res = tag_res or row.res
        prefix = row.prefix or OP_PREFIX[generic]
        orders = [(row.a, row.b, False)] + ([(row.b, row.a, True)] if row.mirror else [])
        for a, b, swapped in orders:
            tag_a, type_a = _op_concrete(a, dim, kind)
            tag_b, type_b = _op_concrete(b, dim, kind)
            name = '_'.join(t for t in (prefix, tag_a, tag_b, kind) if t)
            out.append(OpEntry(generic, name, type_a, type_b, res, swapped, tag_a, tag_b, tag_res,
                               dim, kind))
    return out


def op_layout(tag):
    """Layout of a tensor tag of an ``OpEntry`` (``'3D2Osym'`` -> ``layout('2Osym', '3D')``)."""
    if _op_dim(tag) is None:
        raise ValueError(f"'{tag}' is not a tensor tag with a dimension")
    return layout(tag[2:], tag[:2])


def op_uses(module, kinds):
    """Types used by the operator module ``module`` (key of ``OP_MODULES``) for ``kinds``, as
    ``[(fortran_module, [type names])]``, one entry per module and kind (the kindless identities
    in their own entry) so that each ``use`` line stays short: ``op_uses('assignment', ['dp'])``
    -> ``[('muscle_tensor_2o', ['ten_3D2O_dp', 'ten_2D2O_dp']), ...]``."""
    if module not in OP_MODULES:
        raise ValueError(f"unknown operator module '{module}', expected one of {tuple(OP_MODULES)}")
    groups = {}
    for generic in OP_MODULES[module]:
        for kind in kinds:
            for dim in DIMS:
                for e in op_entries(generic, dim, kind):
                    for tag, typ in ((e.tag_a, e.a), (e.tag_b, e.b), (e.tag_res, e.res)):
                        if typ is None or typ.startswith('real('):
                            continue
                        iden = op_iden(tag)
                        if iden is None:
                            key = (family_module(tag[2:]), kind)
                        else:
                            key = (family_module('iden'), kind if iden[1] else None)
                        names = groups.setdefault(key, [])
                        if typ not in names:
                            names.append(typ)
    return [(mod, names) for (mod, _), names in groups.items()]


#: Argument name of each identity in the operator procedures (tensors are ``a``, ``b``).
_OP_IDEN_ARG = {'2O': 'I2', '4O3T': 'I3', '4O4T': 'I4'}

#: LaTeX symbol of each identity (unscaled).
_OP_IDEN_SYMBOL = {'2O': r'\mathbf{I}', '4O3T': r'\mathbb{I}_{3T}', '4O4T': r'\mathbb{I}^S'}


def _op_operands(e):
    """Tags of the operands of the entry ``e`` (one for a unary operator)."""
    return [t for t in (e.tag_a, e.tag_b) if t is not None]


def _op_order(tag):
    """Tensor order (2 or 4) of a tensor or identity tag."""
    iden = op_iden(tag)
    family = iden[0] if iden else tag[2:]
    return int(family[0])


def op_args(e):
    """Dummy argument names of the operator entry ``e`` (``OpEntry``), one per operand: the
    tensors are ``a`` and ``b`` (in order), the identities ``I2`` (``2O``), ``I3`` (``4O3T``) and
    ``I4`` (``4O4T``): ``sum_I2O_3D2Osym_dp`` -> ``['I2', 'a']``, ``sum_I4O3TS_I4O4T_dp`` ->
    ``['I3', 'I4']``, ``dot_3D2O_3D2Osym_dp`` -> ``['a', 'b']``."""
    letters = iter('ab')
    return [_OP_IDEN_ARG[op_iden(t)[0]] if op_iden(t) else next(letters) for t in _op_operands(e)]


def op_symbols(e, letters='AB'):
    r"""LaTeX symbols ``[sym_a, (sym_b,) sym_res]`` of the operands and the result of the entry
    ``e`` for the documentation. Tensors take ``letters`` in order, bold (``\mathbf``) for 2nd
    order and blackboard (``\mathbb``) for 4th order; identities their usual symbols
    (``\mathbf{I}``, ``\mathbb{I}_{3T}``, ``\mathbb{I}^S``), preceded by the factor ``c`` if
    scaled (``c_1``, ``c_2`` if both operands are scaled identities); the result is ``res``
    (``None`` for ``=``): ``sum_I2OS_3D2O_dp`` ->
    ``['c\mathbf{I}', '\mathbf{A}', '\mathbf{res}']``."""
    tags = _op_operands(e)
    scaled = [t for t in tags if op_iden(t) and op_iden(t)[1]]
    factors = iter(['c_1', 'c_2'] if len(scaled) == 2 else ['c'])
    names = iter(letters)
    out = []
    for t in tags:
        iden = op_iden(t)
        if iden:
            out.append((next(factors) if iden[1] else '') + _OP_IDEN_SYMBOL[iden[0]])
        else:
            out.append((r'\mathbf{{{}}}' if _op_order(t) == 2 else r'\mathbb{{{}}}').format(next(names)))
    if e.tag_res == 'real':
        out.append('res')
    elif e.tag_res is not None:
        out.append(r'\mathbf{res}' if _op_order(e.tag_res) == 2 else r'\mathbb{res}')
    else:
        out.append(None)
    return out


def op_doc(e, text):
    """One-line FORD summary of the entry ``e``: ``text`` (no final period; runs of spaces, e.g.
    from fypp line continuations, count as one) followed by the dimension and kind of the entry:
    ``op_doc(e, 'Sum')`` -> ``'Sum (3D, `real(dp)`).'``. The line is a comment, so it may run past
    the free-form limit."""
    return '{} ({}, `real({})`).'.format(' '.join(text.split()), e.dim, e.kind)


def op_iden_term(tag, var, sign='+'):
    """Term ``(name, var, sign)`` of ``iden_combination`` for the identity operand ``tag`` whose
    factor is the Fortran expression ``var`` (dropped, ``None``, if the identity is unscaled):
    ``op_iden_term('I2OS', 'I2%val')`` -> ``('2O', 'I2%val', '+')``,
    ``op_iden_term('I4O4T', 'I4%val', '-')`` -> ``('4O4T', None, '-')``."""
    iden = op_iden(tag)
    if iden is None:
        raise ValueError(f"'{tag}' is not an identity")
    return iden[0], var if iden[1] else None, sign


def op_mirror(e):
    """Entry of the row whose mirror is ``e`` (same generic, dimension and kind, operands
    swapped): for ``tdotsym_I2O_3D2Osym_dp`` -> the entry ``tdotsym_3D2Osym_I2O_dp``."""
    if not e.swapped:
        raise ValueError(f"'{e.name}' is not a mirrored entry")
    return [x for x in op_entries(e.generic, e.dim, e.kind)
            if not x.swapped and (x.tag_a, x.tag_b) == (e.tag_b, e.tag_a)][0]


def _op_dot(la, a, lb, b, i, j):
    """``[(k, 'a_ik*b_kj')]`` for k = 1, 2, 3 in order, without the structural zeros."""
    return [(k, '{}*{}'.format(la.ref(a, i, k), lb.ref(b, k, j)))
            for k in (1, 2, 3) if la.index(i, k) and lb.index(k, j)]


def op_dot_terms(la, a, lb, b, i, j):
    """Terms of the single contraction \\( A_{ik} B_{kj} \\) (k = 1, 2, 3 in order) of the 2nd-order
    tensors stored in ``a`` (layout ``la``) and ``b`` (layout ``lb``), without the structural
    zeros of 2D: ``op_dot_terms(layout('2O', '2D'), 'a%vals', layout('2Osym', '2D'), 'b%vals',
    3, 3)`` -> ``['a%vals(5)*b%vals(3)']``."""
    return [term for _, term in _op_dot(la, a, lb, b, i, j)]


def op_sym_sum(lay, var, i, j):
    """\\( A_{ij} + A_{ji} \\) of the 2nd-order tensor stored in ``var`` (layout ``lay``), in
    parentheses, or the single reference \\( A_{ii} \\) if ``i == j`` (the shear pairs of a double
    contraction with a symmetric tensor): ``op_sym_sum(layout('2O', '3D'), 'a%vals', 1, 2)`` ->
    ``'(a%vals(4) + a%vals(2))'``."""
    if i == j:
        return lay.ref(var, i, i)
    return '({} + {})'.format(lay.ref(var, i, j), lay.ref(var, j, i))


def op_weight(x, w, kind):
    """Factor ``x`` of a double contraction in Voigt storage times its shear weight ``w``
    (``Comp.mult``: 1 or 2): ``op_weight('b%vals(4)', 2, 'dp')`` -> ``'2.0_dp*b%vals(4)'``,
    ``op_weight('b%vals(1)', 1, 'dp')`` -> ``'b%vals(1)'``."""
    return x if w == 1 else f'{float(w)!r}_{kind}*{x}'


def op_weighted(x, y, w, kind):
    """Product ``x*y`` of a double contraction in Voigt storage with the shear weight ``w``
    (``Comp.mult``: 1 or 2) on the right factor, ``x*(w*y)`` (one form for every kernel and
    dimension, D-015.8); ``w == 1`` gives ``x*y``:
    ``op_weighted('a%vals(4)', 'b%vals(4)', 2, 'dp')`` -> ``'a%vals(4)*(2.0_dp*b%vals(4))'``."""
    return f'{x}*{y}' if w == 1 else f'{x}*({op_weight(y, w, kind)})'


def op_ddot_voigt_terms(l4, a, b, I, kind):
    """Terms of \\( (\\mathbb{A} : \\mathbf{b})_I = \\sum_J A_{IJ} w_J b_J \\) (J ascending) of the
    4th-order tensor stored in ``a`` (layout ``l4``, rows and columns in the storage of its base
    2nd-order layout) and the ``2Osym`` tensor stored in ``b``; ``w_J`` is the shear weight
    (``op_weighted``): ``op_ddot_voigt_terms(layout('4O3sym', '2D'), 'a%vals', 'b%vals', 1,
    'dp')[3]`` -> ``'a%vals(10)*(2.0_dp*b%vals(4))'``."""
    base = l4.base.comps
    return [op_weighted(l4.ref(a, *base[I - 1].idx, *c.idx), l4.base.ref(b, *c.idx), c.mult, kind)
            for c in base]


def op_ddot_matrix_terms(la, a, lb, b, I, J, kind):
    """Terms of \\( (\\mathbb{A} : \\mathbb{B})_{IJ} = \\sum_K A_{IK} w_K B_{KJ} \\) (K ascending) of the
    4th-order tensors stored in ``a`` (layout ``la``) and ``b`` (layout ``lb``), with Voigt rows and
    columns in the storage of their base 2nd-order layout and the shear weight ``w_K`` on ``B``
    (``op_weighted``): ``op_ddot_matrix_terms(layout('4O3sym', '2D'), 'a%vals', layout('4O3sym',
    '2D'), 'b%vals', 1, 1, 'dp')[3]`` -> ``'a%vals(10)*(2.0_dp*b%vals(10))'``."""
    base = la.base.comps
    return [op_weighted(la.ref(a, *base[I - 1].idx, *c.idx), lb.ref(b, *c.idx, *base[J - 1].idx),
                        c.mult, kind) for c in base]


def op_congruence(dim, f, s, t='T'):
    """Terms of the congruence \\( R = F S F^T \\) of the ``2Osym`` tensor stored in ``s`` by the
    ``2O`` tensor stored in ``f``, in the order of the historical unrolled code: an intermediate
    \\( T_{iL} = F_{iK} S_{KL} \\) (K = 1, 2, 3) stored in the scalars ``<t><i><L>``, then
    \\( R_{ij} = T_{iL} F_{jL} \\) (L = 1, 2, 3); structural zeros (2D) are dropped. A component
    made of a single product (the out-of-plane 33 in 2D) is written ``F_iK*F_jL*S_KL`` without the
    intermediate, as the original. Returns ``(tmp, res)``: ``tmp`` = ``[((i, L), expr), ...]`` (the
    entries of ``t`` that are used) and ``res`` = ``{(i, j): expr}`` for the stored components of
    the ``2Osym`` result."""
    lf, ls = layout('2O', dim), layout('2Osym', dim)

    def fs(i, L):
        return _op_dot(lf, f, ls, s, i, L)

    tmp, res = {}, {}
    for c in ls.comps:
        i, j = c.idx
        outer = [L for L in (1, 2, 3) if lf.index(j, L) and fs(i, L)]
        if len(outer) == 1 and len(fs(i, outer[0])) == 1:
            L = outer[0]
            K = fs(i, L)[0][0]
            res[c.idx] = '{}*{}*{}'.format(lf.ref(f, i, K), lf.ref(f, j, L), ls.ref(s, K, L))
            continue
        for L in outer:
            tmp[(i, L)] = ' + '.join(term for _, term in fs(i, L))
        res[c.idx] = ' + '.join('{}{}{}*{}'.format(t, i, L, lf.ref(f, j, L)) for L in outer)
    return sorted(tmp.items()), res


def op_voigt_transform(dim, var, kind):
    """Voigt matrix of the transformation of a 4th-order tensor with minor symmetries by the
    ``2O`` tensor stored in ``var``: \\( R = T C T^T \\) with
    \\( T_{(ij),(ab)} = A_{ia} A_{jb} + [a \\neq b] A_{ib} A_{ja} \\), rows (ij) and columns (ab)
    in the ``2Osym`` storage order of ``dim``. Returns ``{(I, J): expression}`` for the entries
    that are not structurally zero; for ``i == j`` the two equal products are written
    ``2.0_<kind>*A_ia*A_ib`` (the historical expression):
    ``op_voigt_transform('3D', 'A%vals', 'dp')[(1, 4)]`` -> ``'2.0_dp*A%vals(1)*A%vals(4)'``."""
    lay = layout('2O', dim)
    base = layout('2Osym', dim).comps
    out = {}
    for r in base:
        i, j = r.idx
        for c in base:
            a, b = c.idx
            pairs = [(i, a, j, b)] + ([(i, b, j, a)] if a != b else [])
            terms = ['{}*{}'.format(lay.ref(var, w, x), lay.ref(var, y, z))
                     for w, x, y, z in pairs if lay.index(w, x) and lay.index(y, z)]
            if not terms:
                continue
            if len(terms) == 2 and i == j:
                expr = '2.0_{}*{}'.format(kind, terms[0])
            else:
                expr = ' + '.join(terms)
            out[(r.pos[0], c.pos[0])] = expr
    return out


VoigtCongruence = namedtuple('VoigtCongruence', 'n T rows w res')
VoigtCongruence.__doc__ = """Scalar form of \\( R = T C T^T \\) on the Voigt matrices
(``op_voigt_congruence``).

n:    size of the Voigt matrices.
T:    ``op_voigt_transform``: ``{(I, J): expr}`` for the entries of T that are not zero.
rows: ``{I: [K, ...]}`` columns of the non-zero entries of row I of T, ascending.
w:    ``[((I, K), terms)]`` entries \\( W_{IK} = T_{IL} C_{LK} \\) (L in ``rows[I]``) needed by
      the stored components of R, sorted, with C read from its own storage.
res:  ``[(pos, terms)]`` stored components \\( R_{IJ} = W_{IK} T_{JK} \\) (K in ``rows[J]``), in
      storage order.
"""


def op_voigt_congruence(lay, a, c, kind, t='T', w='W'):
    """Terms of \\( R = T C T^T \\) (``op_voigt_transform``) for the 4th-order layout ``lay``
    (``4O2sym`` or ``4O3sym``) of C and R, with A stored in ``a`` and C in ``c``, as scalars:
    the entries of T are named ``<t><I><J>``, those of W = T C ``<w><I><K>``; the sums keep the
    order of the original kernels (L, K ascending). For ``4O3sym`` only the entries of W used by
    the upper triangle of R appear: ``op_voigt_congruence(layout('4O3sym', '2D'), 'A%vals',
    'C_in%vals', 'dp').res[0]`` -> ``('1', ['W11*T11', 'W12*T12', 'W14*T14'])``."""
    T = op_voigt_transform(lay.dim, a, kind)
    base = lay.base.comps
    n = len(base)
    rows = {I: [K for K in range(1, n + 1) if (I, K) in T] for I in range(1, n + 1)}
    need = sorted({(x.voigt[0], K) for x in lay.comps for K in rows[x.voigt[1]]})
    wterms = [((I, K), ['{}{}{}*{}'.format(t, I, L, lay.ref(c, *base[L - 1].idx, *base[K - 1].idx))
                        for L in rows[I]]) for I, K in need]
    res = [(', '.join(map(str, x.pos)),
            ['{}{}{}*{}{}{}'.format(w, x.voigt[0], K, t, x.voigt[1], K) for K in rows[x.voigt[1]]])
           for x in lay.comps]
    return VoigtCongruence(n, T, rows, wterms, res)


def _check_op_table():
    """Validate every row of ``OP_TABLE`` once, at import (families, identities, dimensions,
    mirrors), so a wrong row fails even if no template expands it."""
    for generic, rows in OP_TABLE.items():
        for row in rows:
            if row.mirror and (row.b is None or row.a == row.b):
                raise ValueError(f'{generic} row {row}: a mirror needs two different operands')
            for dim in _op_row_dims(row):
                for x in (row.a, row.b, row.res):
                    _op_concrete(x, dim, KINDS[0])


_check_op_table()


# ---------------------------------------------------------------------------------------------
# Constitutive model families (D-007, D-017.6)
# ---------------------------------------------------------------------------------------------
# A family fixes the names every model of it follows (D-007): the generics ``<prefix>_<func>``
# of its family module, the model modules and types, which functions each model has, the values a
# composed law receives and the joint evaluation. The signatures themselves are written by hand in
# the Fortran templates (the generated target tests keep their own copy, tests/gpu/target_cases.py).

#: Function of a family: ``name`` (generic ``<prefix>_<name>``), ``comp`` (names of the composition
#: values it receives when the model composes another law, D-017.6b, after the other arguments) and
#: ``models`` (the registered models that have this function; ``None``: every model; the corotational
#: kinematics is updated from a strain and a rotation, the others from F).
ModelFunc = namedtuple('ModelFunc', 'name comp models', defaults=((), None))

#: Parameter of a model: ``name``; ``type`` = ``'real'``, ``'integer'``, ``'logical'`` or a tensor
#: family (``'2Osym'``, or with its dimension, ``'3D2Osym'``); ``shape`` = ``''`` (scalar) or a fixed
#: shape such as ``'(3)'``; ``doc`` (FORD text) and ``default`` (``None`` or a number written without
#: kind, e.g. ``'0.0'``; the kind is added for real parameters).
#: ``private`` (optional, default ``False``): the component is private to the model module (derived
#: parameters such as the normalization factor of CPB06, D-017.10b): the type can then only be built
#: with the initializer of the model (``yield_<m>_init``).
ModelParam = namedtuple('ModelParam', 'name type shape doc default private', defaults=(False,))

#: Value of a composed law (D-017.6b): the ``name`` of the argument the free functions receive, and
#: the function of the composed family (``func``, e.g. ``'stress'`` for ``hard_stress``) that the
#: caller (the material, D-017.13d) evaluates at the plastic strain.
CompValue = namedtuple('CompValue', 'name func')

#: Joint evaluation of several functions of a family at the same arguments (yield criteria: the
#: equivalent stress and its two derivatives, which a return mapping needs together): a subroutine
#: ``<prefix>_<name>`` that returns the results of ``funcs``, and a ``flag`` function
#: ``<prefix>_<flag>(prm)`` (constant ``logical``: whether the model writes its own, cheaper joint
#: evaluation; the others get the default, the functions called in turn).
ModelJoint = namedtuple('ModelJoint', 'name flag funcs')

#: A model family: ``prefix`` of every generated name (``hard`` -> ``hard_stress``,
#: ``hard_bilinear_dp``), ``module`` (format of the model modules), ``family_module``, ``funcs``
#: (``ModelFunc``), ``models`` (the registry), ``composition`` (``CompValue`` of a family whose
#: models may compose another law: the viscoplastic laws with a hardening law), ``composed``
#: (registered models that compose the other law), ``dims`` (dimensions of a family with tensor
#: arguments: the specific function names then carry the dimension,
#: ``<prefix>_<func>_<model>_<dim>_<kind>``; ``None`` for a family of scalar arguments) and
#: ``joint`` (a ``ModelJoint`` or ``None``).
ModelFamily = namedtuple('ModelFamily', 'prefix module family_module funcs models composition composed dims joint',
                         defaults=((), (), None, None))

_F = ModelFunc

#: Composition with an isotropic hardening law (viscoplastic and overstress families): the
#: composer evaluates ``hard_stress`` and ``hard_dstress_dep`` at ``ep`` and passes the values.
_HARD_COMPOSITION = (CompValue('sy', 'stress'), CompValue('dsy_dep', 'dstress_dep'))

#: Model families by prefix.
MODEL_FAMILIES = {
    # Isotropic hardening laws sigma_y(ep).
    'hard': ModelFamily(
        prefix='hard', module='muscle_hard_{}', family_module='muscle_hard_laws', models=HARDENING_MODELS,
        funcs=(_F('stress'), _F('dstress_dep'), _F('ddstress_ddep'))),
    # Flow stress laws sigma(ep, epd): the partial derivatives, and dstress_dep, the derivative used by
    # the return mapping, d sigma(ep_n + x, x/dt)/dx (with the time step); newton_dstress_dep only in
    # the laws with a rate floor (D-018.19b).
    'vp': ModelFamily(
        prefix='vp', module='muscle_vp_{}', family_module='muscle_vp_laws', models=VISCOPLASTIC_MODELS,
        funcs=(_F('flow_stress', ('sy',)), _F('dflow_dep', ('sy', 'dsy_dep')), _F('dflow_depd', ('sy',)),
               _F('dstress_dep', ('sy', 'dsy_dep')),
               _F('newton_dstress_dep', ('sy', 'dsy_dep'), models=VISCOPLASTIC_FLOOR)),
        composition=_HARD_COMPOSITION, composed=VISCOPLASTIC_COMPOSED),
    # Overstress laws g(sigma_eq, ep) (plastic multiplier rate).
    'os': ModelFamily(
        prefix='os', module='muscle_os_{}', family_module='muscle_os_laws', models=OVERSTRESS_MODELS,
        funcs=(_F('flow_rate', ('sy',)), _F('dflow_dsigma', ('sy',)), _F('dflow_dep', ('sy', 'dsy_dep'))),
        composition=_HARD_COMPOSITION, composed=OVERSTRESS_COMPOSED),
    # Yield criteria sigma_eq(sigma) and their derivatives with respect to the stress tensor. Only 3D
    # for now (dims-2d-3d.md: 2D needs a 2D vs 3D test first).
    'yield': ModelFamily(
        prefix='yield', module='muscle_yield_{}', family_module='muscle_yield_laws', models=YIELD_MODELS,
        funcs=(_F('stress_eq'), _F('dstress_eq'), _F('ddstress_eq')), dims=('3D',),
        joint=ModelJoint('stress_eq_and_derivatives', 'joint_evaluation', ('stress_eq', 'dstress_eq', 'ddstress_eq'))),
    # Elastic models sigma(eps) (D-017.12), in 3D and 2D (plane strain/stress and axisymmetry: the 33
    # component is kept). A model with internal variables belongs to the viscoelastic family.
    'elas': ModelFamily(
        prefix='elas', module='muscle_elas_{}', family_module='muscle_elas_laws',
        models=ELASTICITY_MODELS, dims=('3D', '2D'),
        funcs=(_F('stress'), _F('stress_dev'), _F('dstress_dstrain'), _F('dstress_dstrain_dev'))),
    # Linear viscoelastic models (D-017.13), integrated over a time step, with an explicit state
    # (muscle_visco_types). `visco_step_coefficients` once per step and material (it validates the time
    # step), `visco_freeze` once per step and point (the frozen elasticity `elas_frozen_<k>`), `visco_update`
    # once per converged step and point; the `_branches` variants work on a state of any size (D-019.7b).
    # Only 3D.
    'visco': ModelFamily(
        prefix='visco', module='muscle_visco_{}', family_module='muscle_visco_laws', models=VISCOELASTIC_MODELS,
        funcs=(_F('step_coefficients'), _F('freeze'), _F('update'), _F('n_branches'), _F('freeze_branches'),
               _F('update_branches'))),
    # Kinematics (D-017.14): the update of the kinematic state (from F, or from a strain and a rotation
    # in the corotational kinematics, D-017.14d), the conversion of the constitutive stress to Cauchy and
    # of the constitutive tangent to the spatial one. No material parameters. Only 3D.
    'kin': ModelFamily(
        prefix='kin', module='muscle_kin_{}', family_module='muscle_kin_laws', models=KINEMATICS_MODELS,
        funcs=(_F('update', models=KINEMATICS_F_MODELS), _F('update_norot', models=KINEMATICS_F_MODELS),
               _F('corot_update', models=('corotational',)), _F('to_cauchy'), _F('to_spatial_tangent'),
               _F('geometric_tangent'), _F('is_finite_strain'))),
}


# ---------------------------------------------------------------------------------------------
# Parameters and initializers of the registered models (D-018.17e: single source)
# ---------------------------------------------------------------------------------------------
# The model templates write their parameter type and initializer by hand (a test checks them against
# ``MODEL_PARAMS``, ``MODEL_INITS``);
# the generated material initializers (``mat_init_<name>``, flat parameter array of an FE code) use
# the same inputs (``model_inputs``).

_P = ModelParam

#: Parameters of every registered model, ``MODEL_PARAMS[(prefix, model)]`` (``ModelParam``, storage
#: order = order of the positional structure constructor).
MODEL_PARAMS = MappingProxyType({
    ('elas', 'frozen'): (
        _P('k_alg', 'real', '', r'Algorithmic bulk modulus \( K^{alg} \).', '0.0'),
        _P('g_alg', 'real', '', r'Algorithmic shear modulus \( G^{alg} \).', '0.0'),
        _P('stress_hist', '3D2Osym', '',
           r'History stress \( \boldsymbol\sigma^{hist}_n \) (constant over the step).', '0.0'),
    ),
    ('elas', 'linear'): (
        _P('lambda', 'real', '', r'Lamé first parameter \( \lambda = E\nu/[(1+\nu)(1-2\nu)] \).', '0.0', True),
        _P('mu', 'real', '', r'Shear modulus \( \mu = E/[2(1+\nu)] \).', '0.0', True),
        _P('e1', 'real', '', r'Diagonal constant \( e_1 = \lambda + 2\mu \) of Hooke law.', '0.0', True),
    ),
    ('yield', 'cpb06'): (
        _P('C1', 'real', '(3,3)',
           r'Normal block \( C_{ij} \), i, j = 1..3, of the linear transformation \( \mathbf{C} \).', None),
        _P('C2', 'real', '(3)',
           r'Shear factors \( (C_{44}, C_{55}, C_{66}) \) of \( \mathbf{C} \) (xy, yz, xz).', None),
        _P('k', 'real', '', r'Strength differential parameter \( k \in [-1, 1] \).', None),
        _P('a', 'real', '', r'Exponent \( a \ge 1 \).', None),
        _P('B', 'real', '',
           r'Normalization factor \( B \) (Cazacu et al., 2006, Eq. 12), from C1, k and a.', None, True),
        _P('L', 'real', '(3,3)',
           r'Normal block of \( \mathbf{L} = \mathbf{C} \mathbb{P}_{dev} \) (Eq. 8), from C1.', None, True),
    ),
    ('yield', 'cpb06ex2'): (
        _P('L', 'real', '(3,3)',
           r'Normal block of \( \mathbf{L} = \mathbf{C} \mathbb{P}_{dev} \), which maps the normal stresses '
           r'to those of \( \boldsymbol\Sigma \) (Plunkett et al., 2008, Eq. 9).', None, True),
        _P('Lp', 'real', '(3,3)',
           r'Normal block of \( \mathbf{L}^\prime = \mathbf{C}^\prime \mathbb{P}_{dev} \), for \( '
           r'\boldsymbol\Sigma^\prime \) (Eq. 9; suffix p = primed in the paper).', None, True),
        _P('c_shear', 'real', '(3)',
           r'Shear factors of \( \mathbf{C} \) on the tensorial shears (xy, yz, xz).', None),
        _P('cp_shear', 'real', '(3)',
           r'Shear factors of \( \mathbf{C}^\prime \) on the tensorial shears (xy, yz, xz).', None),
        _P('k', 'real', '', r'Strength differential parameter \( k \) of \( \boldsymbol\Sigma \) (Eq. 8).', None),
        _P('kp', 'real', '',
           r'Strength differential parameter \( k^\prime \) of \( \boldsymbol\Sigma^\prime \) (Eq. 8).', None),
        _P('a', 'real', '', r'Degree of homogeneity \( a \) (Eq. 8).', None),
        _P('B', 'real', '',
           r'Normalization: \( \sigma_{eq} = s \) in uniaxial tension \( s \) along RD (Eq. 11 with \( '
           r'\theta = 0 \)).', None, True),
    ),
    ('yield', 'druckerprager'): (
        _P('f_factor', 'real', '',
           r'Calibration factor \( f \) of the deviatoric term (hardening mode).', '0.0', True),
        _P('I1_factor', 'real', '', r'Coefficient \( \tan\beta \, f \) of the pressure \( p \).', '0.0', True),
        _P('K', 'real', '',
           r'Ratio \( K \) of the yield stresses in triaxial tension and compression (\( 0.778 \le K \le 1 '
           r'\)).', '0.0'),
        _P('hardening_mode', 'integer', '',
           r'Hardening calibration mode (`DP_HARDENING_TENSION`, `DP_HARDENING_COMPRESSION` or '
           r'`DP_HARDENING_SHEAR`).', 'DP_HARDENING_TENSION'),
    ),
    ('yield', 'hill48'): (
        _P('f', 'real', '', r'Anisotropy coefficient \( F \) (yy-zz).', None),
        _P('g', 'real', '', r'Anisotropy coefficient \( G \) (zz-xx).', None),
        _P('h', 'real', '', r'Anisotropy coefficient \( H \) (xx-yy).', None),
        _P('l', 'real', '', r'Anisotropy coefficient \( L \) (shear yz).', None),
        _P('m', 'real', '', r'Anisotropy coefficient \( M \) (shear xz).', None),
        _P('n', 'real', '', r'Anisotropy coefficient \( N \) (shear xy).', None),
    ),
    ('yield', 'vonmises'): (),
    ('yield', 'yld2004'): (
        _P('Lp', 'real', '(3,3)',
           r'Normal block of \( \mathbf{L}^\prime = \mathbf{C}^\prime \mathbb{P}_{dev} \), which maps (sxx, '
           r'syy, szz) to (Spxx, Spyy, Spzz).', None, True),
        _P('cp_shear', 'real', '(3)',
           r'Shear factors of \( \mathbf{C}^\prime \) on (xy, yz, xz): (cp44, cp55, cp66).', None),
        _P('Lpp', 'real', '(3,3)',
           r'Normal block of \( \mathbf{L}^{\prime\prime} = \mathbf{C}^{\prime\prime} \mathbb{P}_{dev} \).', None, True),
        _P('cpp_shear', 'real', '(3)',
           r'Shear factors of \( \mathbf{C}^{\prime\prime} \) on (xy, yz, xz): (cpp44, cpp55, cpp66).', None),
        _P('a', 'real', '', r'Exponent \( a \ge 2 \) (the analytical derivatives need it).', None),
    ),
    ('hard', 'bilinear'): (
        _P('y0', 'real', '', r'Initial yield stress \( \sigma_{y0} \).', None),
        _P('K', 'real', '', r'Hardening modulus \( K \).', None),
    ),
    ('hard', 'ludwik'): (
        _P('sigma0', 'real', '', r'Initial yield stress \( \sigma_0 \).', None),
        _P('k', 'real', '', r'Strength coefficient \( K \).', None),
        _P('n', 'real', '', r'Hardening exponent \( n \).', None),
    ),
    ('hard', 'swift'): (
        _P('k', 'real', '', r'Strength coefficient \( K \).', None),
        _P('e0', 'real', '', r'Initial strain offset \( \epsilon_0 \) (> 0, see the module documentation).', None),
        _P('n', 'real', '', r'Hardening exponent \( n \).', None),
    ),
    ('hard', 'vocemod'): (
        _P('Sy', 'real', '', r'Initial yield stress \( \sigma_{y0} \).', None),
        _P('k', 'real', '', r'Linear hardening modulus \( k \) (asymptotic slope).', None),
        _P('q', 'real', '', r'Saturation stress range \( q \) of the exponential term.', None),
        _P('n', 'real', '', r'Saturation rate \( n \).', None),
    ),
    ('vp', 'jc'): (
        _P('C', 'real', '', r'Strain rate sensitivity \( C \).', None),
        _P('epdmax', 'real', '', r'Reference equivalent plastic strain rate \( \dot\epsilon_{p0} \).', None),
        _P('epmin', 'real', '',
           r'Minimum strain rate \( \dot\epsilon_{min} \) (floor of the rate, avoids \( \ln 0 \)).', '1.0e-7'),
    ),
    ('vp', 'jcfull'): (
        _P('C', 'real', '', r'Strain rate sensitivity \( C \).', '0.0'),
        _P('epdmax', 'real', '', r'Reference equivalent plastic strain rate \( \dot\epsilon_{p0} \).', '1.0'),
        _P('m', 'real', '', r'Thermal softening exponent \( m \) (0: no softening).', '0.0'),
        _P('Troom', 'real', '', r'Reference (room) temperature \( T_{room} \) (K).', '293.15'),
        _P('Tmelt', 'real', '', r'Melting temperature \( T_{melt} \) (K).', '1.0e3'),
        _P('T_current', 'real', '', r'Current temperature \( T \) (K), set by the caller.', '0.0'),
    ),
    ('vp', 'mrk'): (
        _P('B01', 'real', '', r'Linear hardening modulus \( B_{01} \).', None),
        _P('B02', 'real', '', r'Saturation hardening modulus \( B_{02} \).', None),
        _P('epdmax', 'real', '', r'Upper limit rate \( \dot\epsilon_{max} \).', None),
        _P('nu1', 'real', '', r'Rate sensitivity of \( B_1 \), \( \nu_1 \).', None),
        _P('nu2', 'real', '', r'Rate sensitivity of \( B_2 \), \( \nu_2 \).', None),
        _P('n0', 'real', '', r'Saturation rate at the reference rate \( n_0 \).', None),
        _P('D2', 'real', '', r'Rate sensitivity of the saturation rate \( D_2 \).', None),
        _P('epdmin', 'real', '', r'Lower limit rate \( \dot\epsilon_{min} \).', None),
        _P('chi1', 'real', '', r'Rate sensitivity of the rate factor \( \chi_1 \).', None),
        _P('chi2', 'real', '', r'Exponent of the rate factor \( \chi_2 \).', None),
        _P('sig_u', 'real', '', r'Athermal stress \( \sigma_u \).', None),
    ),
    ('vp', 'nnl'): (
        _P('sig_a', 'real', '', r'Athermal stress coefficient \( \sigma_a \).', None),
        _P('n1', 'real', '', r'Athermal hardening exponent \( n_1 \).', None),
        _P('sig_0', 'real', '', r'Thermal stress coefficient \( \sigma_0 \).', None),
        _P('KG0', 'real', '', r'Normalized activation coefficient \( k T/G_0 \).', None),
        _P('epd0', 'real', '', r'Reference strain rate \( \dot\epsilon_0 \).', None),
        _P('at', 'real', '', r'Microstructure evolution coefficient \( a \).', None),
        _P('n0', 'real', '', r'Microstructure evolution exponent \( n_0 \).', None),
        _P('q', 'real', '', r'Barrier profile exponent \( q \).', None),
        _P('p', 'real', '', r'Barrier profile exponent \( p \).', None),
        _P('eps_log', 'real', '', r'Floor of the logarithm arguments and of the rate.', '1.0e-12'),
        _P('eps_base', 'real', '', r'Floor of \( 1 - X^{1/q} \).', '1.0e-12'),
    ),
    ('vp', 'rk'): (
        _P('B0', 'real', '', r'Modulus of plasticity \( B_0 \).', None),
        _P('ep0', 'real', '', r'Strain offset \( \epsilon_0 \).', None),
        _P('n0', 'real', '', r'Hardening exponent at the reference rate \( n_0 \).', None),
        _P('D2', 'real', '', r'Rate sensitivity of the hardening exponent \( D_2 \).', None),
        _P('D1', 'real', '', r'Rate sensitivity of the effective stress \( D_1 \).', None),
        _P('nu', 'real', '', r'Rate sensitivity of the modulus \( \nu \).', None),
        _P('epdmin', 'real', '', r'Lower limit rate \( \dot\epsilon_{min} \).', None),
        _P('epdmax', 'real', '', r'Upper limit rate \( \dot\epsilon_{max} \).', None),
        _P('sig0', 'real', '', r'Effective stress at 0 K \( \sigma_0^* \).', None),
        _P('m', 'real', '', r'Exponent of the effective stress \( m \).', None),
        _P('epmin', 'real', '', r'Floor of the rate in the derivatives.', '1.0e-6'),
        _P('epmax', 'real', '', r'Unused (historical component).', '10000.0'),
    ),
    ('vp', 'va'): (
        _P('B', 'real', '', r'Hardening modulus \( B \).', None),
        _P('n', 'real', '', r'Hardening exponent \( n \).', None),
        _P('B1', 'real', '', r'Rate coefficient \( B_1 \).', None),
        _P('B2', 'real', '', r'Rate offset \( B_2 \).', None),
        _P('m', 'real', '', r'Rate exponent \( m \) (0: no rate term).', None),
        _P('sig_u', 'real', '', r'Athermal stress \( \sigma_u \).', None),
    ),
    ('os', 'chaboche'): (
        _P('k', 'real', '', r'Viscous resistance \( K \) (stress time\(^{1/n}\)).', '1.0'),
        _P('n', 'real', '', r'Positive rate exponent \( n \).', '1.0'),
    ),
    ('os', 'perzyna'): (
        _P('eta', 'real', '', r'Viscosity (relaxation time) \( \eta \).', '1.0'),
        _P('n', 'real', '', r'Rate sensitivity exponent \( N \).', '1.0'),
    ),
    ('visco', 'kelvin'): (
        _P('k', 'real', '', r'Elastic bulk modulus \( K \).', '0.0', True),
        _P('g', 'real', '', r'Elastic shear modulus \( G \).', '0.0', True),
        _P('zeta', 'real', '', r'Volumetric viscosity of the dashpot \( \zeta \).', '0.0', True),
        _P('mu_v', 'real', '', r'Deviatoric viscosity of the dashpot \( \mu_v \).', '0.0', True),
    ),
    ('visco', 'maxwell'): (
        _P('k', 'real', '', r'Bulk modulus \( K \).', '0.0', True),
        _P('g_inf', 'real', '', r'Long-term shear modulus \( G_\infty \).', '0.0', True),
        _P('n', 'integer', '', r'Number of Maxwell branches.', '0', True),
        _P('g', 'real', '(MAX_VISCO_BRANCHES)',
           r'Branch shear moduli \( G_k \) (0 in the unused slots).', '0.0', True),
        _P('tau', 'real', '(MAX_VISCO_BRANCHES)',
           r'Branch relaxation times \( \tau_k \) (1 in the unused slots).', '1.0', True),
    ),

    ('kin', 'small_strain'): (),
    ('kin', 'total_lagrangian'): (),
    ('kin', 'updated_lagrangian'): (),
    ('kin', 'material_logarithmic'): (),
    ('kin', 'spatial_logarithmic'): (),
    ('kin', 'corotational'): (),
})

#: Input of the initializer of a model (``ModelInit``): ``name`` (dummy argument, also the keyword of a
#: generated call), ``type`` (``'real'`` or ``'integer'``), ``shape`` (``''`` or ``'(:)'``: one value per
#: branch of a viscoelastic model, assumed shape), ``doc`` (FORD text) and ``domain`` (optional: Fortran
#: condition on the value ``{x}``, literals of kind ``{k}``, that the material initializer from a parameter
#: array checks, ``STATUS_INVALID_PARAMETER`` otherwise; the model initializer itself does not, D-018.17e).
ModelInput = namedtuple('ModelInput', 'name type shape doc domain', defaults=('', None))

#: Initializer ``<prefix>_<model>_init`` of a model whose parameters are not its inputs (derived or
#: validated private components, D-017.10b): ``doc`` (lines of the FORD summary, ``{k}`` = kind),
#: ``inputs`` (``ModelInput``, in the order of the dummy arguments after the parameters ``prm``),
#: ``status`` (FORD text of the final ``integer, intent(out) :: status``: every initializer has one, so that
#: all are called the same way, ``_STATUS_OK_ONLY`` when it validates nothing; ``tests/tools`` checks it).
#: The template writes the signature and the body.
ModelInit = namedtuple('ModelInit', 'doc inputs status', defaults=(None,))

_R = lambda name, doc, shape='': ModelInput(name, 'real', shape, doc)


_SHEAR_AXES = {'44': 'xy', '55': 'yz', '66': 'xz'}


def _matrix_inputs(prefix, rows_cols, symbol):
    """Real inputs ``<prefix><i><j>``, coefficient ``ij`` of the matrix ``symbol`` (LaTeX), for the index
    pairs ``rows_cols`` (44, 55, 66: the shears xy, yz, xz, library order):
    ``[i.name for i in _matrix_inputs('c', ['11', '44'], 'C')]`` -> ``['c11', 'c44']``."""
    return tuple(_R(f'{prefix}{ij}', rf'Coefficient \( {symbol}_{{{ij}}} \)'
                    + (f' (shear {_SHEAR_AXES[ij]})' if ij in _SHEAR_AXES else '') + '.') for ij in rows_cols)


_NORMAL = ('11', '12', '13', '21', '22', '23', '31', '32', '33')
_SHEAR = ('44', '55', '66')
_OFF_DIAG = ('12', '13', '21', '23', '31', '32')
_STATUS = r'`STATUS_OK` or `STATUS_INVALID_PARAMETER`.'
_STATUS_OK_ONLY = (r'Always `STATUS_OK` (no input is validated): the argument exists so that every model '
                   r'initializer is called the same way.')

#: Initializers of the models built with one (``MODEL_INITS[(prefix, model)]``).
MODEL_INITS = MappingProxyType({
    ('elas', 'linear'): ModelInit(
        doc=(r'Parameters of the linear elastic model in `real({k})` from the Young modulus and the Poisson',
             r'ratio: the Lamé parameters \( \lambda \), \( \mu \) and \( e_1 = \lambda + 2\mu \).',
             r'Needs \( \nu < 1/2 \) (D-017.12b).'),
        inputs=(ModelInput('young', 'real', '', r'Young modulus \( E > 0 \).', '{x} > 0.0_{k}'),
                ModelInput('poisson', 'real', '', r'Poisson ratio \( -1 < \nu < 1/2 \).',
                           '{x} > -1.0_{k} .and. {x} < 0.5_{k}')),
        status=_STATUS_OK_ONLY),
    ('yield', 'druckerprager'): ModelInit(
        doc=(r'Parameters of the Drucker-Prager model in `real({k})` matching the Abaqus material',
             r'options: the calibration factor f of the hardening mode and the pressure coefficient',
             r'tan(beta) f. An unknown mode keeps the previous f (as the original `init`; 0 for a',
             r'default-initialized `prm`), so `prm` is `intent(inout)`.'),
        inputs=(_R('beta_deg', r'Friction angle beta in degrees.'),
                _R('K', r'Yield stress ratio K (0.778 <= K <= 1).'),
                ModelInput('hardening_mode', 'integer', '',
                           r'Hardening calibration mode (`DP_HARDENING_*`: 1 tension, 2 compression, 3 shear).',
                           '{x} >= 1 .and. {x} <= 3')),
        status=_STATUS_OK_ONLY),
    ('yield', 'cpb06'): ModelInit(
        doc=(r'Parameters of the CPB06 model in `real({k})` from the components of the linear',
             r'transformation C, k and a; computes the normalization factor B (Cazacu et al., 2006,',
             r'Eq. 12) and the normal block L of C P_dev, so that they are not computed in every call.'),
        inputs=_matrix_inputs('c', _NORMAL + _SHEAR, 'C')
        + (_R('k', r'Strength differential parameter.'), _R('a', r'Exponent.')),
        status=_STATUS_OK_ONLY),
    ('yield', 'cpb06ex2'): ModelInit(
        doc=(r'Parameters of the CPB06ex2 model in `real({k})`: the maps of Plunkett et al. (2008),',
             r'Eqs. (9)-(10), and the normalization B of Eq. (8) to unit uniaxial tension in RD, as in',
             r'Eq. (11), so that they are not computed in every call.'),
        inputs=_matrix_inputs('c', _NORMAL + _SHEAR, 'C')
        + (_R('k', r'Strength differential parameter of Sigma.'),)
        + _matrix_inputs('cp', _NORMAL + _SHEAR, r'C^\prime')
        + (_R('kp', r'Strength differential parameter of Sigmap.'), _R('a', r'Degree of homogeneity.')),
        status=_STATUS_OK_ONLY),
    ('yield', 'yld2004'): ModelInit(
        doc=(r'Parameters of the Yld2004-18p model in `real({k})`: Lp = Cp P_dev and Lpp = Cpp P_dev',
             r'(Barlat et al., 2005, Eqs. 11 and 15), so that each call applies them to the stress',
             r'without forming the deviator. cpIJ = cp_IJ and cppIJ = cpp_IJ; the shears follow the',
             r'library order: 44 -> xy, 55 -> yz, 66 -> xz.'),
        inputs=_matrix_inputs('cp', _OFF_DIAG + _SHEAR, r'C^\prime')
        + _matrix_inputs('cpp', _OFF_DIAG + _SHEAR, r'C^{\prime\prime}')
        + (_R('a', r'Exponent.'),),
        status=_STATUS_OK_ONLY),
    ('visco', 'maxwell'): ModelInit(
        doc=(r'Parameters of the generalized Maxwell model in `real({k})` from the Prony series with',
             r'absolute moduli. `STATUS_INVALID_PARAMETER` (and the default parameters) if `shear` and',
             r'`tau` differ in size, have more than `MAX_VISCO_BRANCHES` branches or a `tau <= 0`.'),
        inputs=(_R('bulk', r'Bulk modulus \( K \).'), _R('shear_inf', r'Long-term shear modulus \( G_\infty \).'),
                _R('shear', r'Branch shear moduli \( G_k \).', '(:)'),
                _R('tau', r'Branch relaxation times \( \tau_k \).', '(:)')),
        status=_STATUS),
    ('visco', 'kelvin'): ModelInit(
        doc=(r'Parameters of the Kelvin-Voigt solid in `real({k})`: elastic moduli and viscosities.',
             r'`STATUS_INVALID_PARAMETER` (and the default parameters) if any of them is negative.'),
        inputs=(_R('bulk', r'Bulk modulus \( K \).'), _R('shear', r'Shear modulus \( G \).'),
                _R('bulk_viscosity', r'Volumetric viscosity \( \zeta \).'),
                _R('shear_viscosity', r'Deviatoric viscosity \( \mu_v \).')),
        status=_STATUS),
})


def model_param_list(prefix, model):
    """Parameters of a registered model (``MODEL_PARAMS``, as ``ModelParam``): what its template passes
    to the type components: ``[p.name for p in model_param_list('hard', 'bilinear')]`` ->
    ``['y0', 'K']``."""
    try:
        return MODEL_PARAMS[(prefix, model)]
    except KeyError:
        raise ValueError(f'no parameters registered for model {model!r} of family {prefix!r}') from None


def model_init(prefix, model):
    """``ModelInit`` of a model built with an initializer, or ``None`` (built by setting its public
    parameter components)."""
    return MODEL_INITS.get((prefix, model))


def model_inputs(prefix, model):
    """Inputs that build a model (``ModelInput``): those of its initializer (``MODEL_INITS``), or else its
    parameters, every one public (set as components): ``[i.name for i in model_inputs('elas',
    'linear')]`` -> ``['young', 'poisson']``."""
    init = model_init(prefix, model)
    if init is not None:
        return init.inputs
    params = model_param_list(prefix, model)
    private = [p.name for p in params if p.private]
    if private:
        raise ValueError(f'model {model!r} of family {prefix!r} has private parameters {private} but no '
                         f'initializer in MODEL_INITS')
    return tuple(ModelInput(p.name, p.type, p.shape, p.doc) for p in params)


def model_family(prefix):
    """The ``ModelFamily`` of a prefix: ``model_family('hard').family_module`` -> ``'muscle_hard_laws'``."""
    return MODEL_FAMILIES[prefix]


def model_func(prefix, name):
    """The ``ModelFunc`` ``name`` of a family."""
    for f in MODEL_FAMILIES[prefix].funcs:
        if f.name == name:
            return f
    raise KeyError(f'family {prefix!r} has no function {name!r}')


def model_label(prefix, model):
    """Readable name of a model for the docs (``MODEL_LABELS``, else the identifier):
    ``model_label('hard', 'vocemod')`` -> ``'modified Voce'``."""
    return MODEL_LABELS.get((prefix, model), model)


def model_module(prefix, model):
    """Module of a model (``ModelFamily.module``): ``model_module('elas', 'linear')`` -> ``muscle_elas_linear``,
    ``model_module('os', 'perzyna')`` -> ``muscle_os_perzyna``."""
    return MODEL_FAMILIES[prefix].module.format(model)


def model_type(prefix, model, kind):
    """Parameter type of a model: ``model_type('hard', 'bilinear', 'dp')`` -> ``hard_bilinear_dp``."""
    return f'{prefix}_{model}_{kind}'


def model_specific(prefix, func, model, kind, dim=None):
    """Specific procedure name, generic + model (+ dimension) + kind (D-011.4):
    ``model_specific('hard', 'stress', 'bilinear', 'dp')`` -> ``hard_stress_bilinear_dp``. In a
    family with dimensions (``ModelFamily.dims``) the procedures with tensor arguments carry it,
    ``model_specific('yield', 'stress_eq', 'hill48', 'dp', '3D')`` ->
    ``yield_stress_eq_hill48_3D_dp``; the flag of the joint evaluation does not (its only argument
    is the parameter type, which does not depend on the dimension)."""
    fam = MODEL_FAMILIES.get(prefix)
    if fam is not None and fam.dims and not (fam.joint and func == fam.joint.flag):
        if dim not in fam.dims:
            raise ValueError(f'{prefix}_{func}: dimension {dim!r} not in {fam.dims}')
        return f'{prefix}_{func}_{model}_{dim}_{kind}'
    return f'{prefix}_{func}_{model}_{kind}'


def model_dims(prefix):
    """Dimensions to instantiate the functions of a family for: ``ModelFamily.dims``, or
    ``[None]`` for a family of scalar arguments (one instance, no dimension)."""
    return list(MODEL_FAMILIES[prefix].dims or [None])


def model_procs(prefix, model=None):
    """Kernel procedures of a family, as the suffixes of its generics ``<prefix>_<name>``: the
    functions, then the joint evaluation and its flag (if the family has one):
    ``model_procs('yield')`` -> ``['stress_eq', 'dstress_eq', 'ddstress_eq',
    'stress_eq_and_derivatives', 'joint_evaluation']``. With ``model``, only the functions that
    model has (``ModelFunc.models``): ``model_procs('kin', 'corotational')`` ->
    ``['corot_update', 'to_cauchy', 'to_spatial_tangent', 'is_finite_strain']``."""
    fam = MODEL_FAMILIES[prefix]
    res = [f.name for f in fam.funcs if model is None or f.models is None or model in f.models]
    if fam.joint:
        res += [fam.joint.name, fam.joint.flag]
    return res


def _dim_family(typ):
    """Whether ``typ`` is a tensor family with its dimension (``'3D4O3sym'``)."""
    return typ[:2] in DIMS and typ[2:] in FAMILIES


def model_params(params):
    """Parameters of a model as ``ModelParam`` from tuples ``(name, type, shape, doc, default[,
    private])``, checking them (types, unique names)."""
    res = tuple(ModelParam(*p) for p in params)
    names = [p.name.lower() for p in res]
    if len(set(names)) != len(names):
        raise ValueError(f'repeated parameter names: {names}')
    for p in res:
        if p.type not in ('real', 'integer', 'logical') and p.type not in FAMILIES and not _dim_family(p.type):
            raise ValueError(f'parameter {p.name!r}: unknown type {p.type!r}')
    return res


def model_composed(prefix, model):
    """Whether a model of the family composes another law (``ModelFamily.composed``):
    ``model_composed('vp', 'jc')`` -> ``True``, ``model_composed('vp', 'rk')`` -> ``False``."""
    return model in MODEL_FAMILIES[prefix].composed


def model_comp_values(prefix, func):
    """``CompValue`` of the composition values a function receives (``ModelFunc.comp``), in
    order: ``[v.name for v in model_comp_values('vp', 'dflow_dep')]`` -> ``['sy', 'dsy_dep']``."""
    values = {v.name: v for v in MODEL_FAMILIES[prefix].composition}
    return [values[n] for n in model_func(prefix, func).comp]


def literal_value(text):
    """Value of a real literal with or without kind (``'1.0'``, ``'1.0_dp'``, ``'2e-3_sp'``), or
    ``None`` if ``text`` is ``None`` or any other expression (a variable, ``max(a, b)``, ...)."""
    if text is None:
        return None
    try:
        return float(str(text).strip().split('_')[0])
    except ValueError:
        return None


# ---------------------------------------------------------------------------------------------------
# Declared materials (Phase 7, D-008, D-018)
# ---------------------------------------------------------------------------------------------------

#: A declared material (D-008.1, D-018.9): a combination of registered models solved by one return
#: mapping. ``name`` (lower case Fortran identifier: module ``muscle_mat_<name>``, type
#: ``mat_<name>_<k>``), ``solver`` (a key of ``MATERIAL_SOLVERS``) and the model of each family
#: the solver needs, by family prefix (``elas``, ``yield_`` -- ``yield`` is a Python keyword --,
#: ``hard``, ``vp``, ``os``, ``visco``; ``None`` = not used); ``doc``: one line for the FORD docs. The
#: parameters are not part of the declaration: they are components of the material type, set at
#: run time. An elasticity of ``MATERIAL_ELAS_PER_POINT`` (the frozen viscoelastic one, ``elas='frozen'``)
#: needs the model of its family (``visco='maxwell'``) and ``visco_branches``: the number of branches
#: the history array of a point keeps (fixed per material, at most ``MAX_VISCO_BRANCHES``; the kernel
#: returns ``STATUS_INVALID_PARAMETER`` for a model with more). ``kin``: a kinematics
#: (``KINEMATICS_MODELS``, D-018.15) makes a finite strain material, whose point kernel takes the
#: deformation gradients at t_n and t_n+1, or, with ``KINEMATICS_FE_ROTATED``, the strain rotated by the
#: FE code (any solver with a consistent tangent: not ``cutting_plane``; not with a viscoelastic
#: elasticity yet; a spatial measure with an anisotropic yield criterion needs an isotropic
#: elasticity, ``ELASTICITY_ISOTROPIC``), and
#: ``tangent_rate`` (``MATERIAL_TANGENT_RATES``, only with ``kin``) chooses its spatial tangent.
Material = namedtuple('Material', 'name solver elas yield_ hard vp os doc visco visco_branches kin tangent_rate',
                      defaults=(None, None, None, None, None, '', None, 0, None, 'none'))

#: Field of ``Material`` of each model family a solver may use.
MATERIAL_FIELDS = {'elas': 'elas', 'yield': 'yield_', 'hard': 'hard', 'vp': 'vp', 'os': 'os', 'visco': 'visco',
                   'kin': 'kin'}

#: Spatial tangents of a finite strain material (D-018.15f): ``none``, ``kin_to_spatial_tangent`` (the
#: part through the constitutive stress), or ``jaumann``, plus ``kin_geometric_tangent``: the tangent of
#: the Jaumann rate of the Kirchhoff stress over J (Abaqus UMAT, ANSYS USERMAT).
MATERIAL_TANGENT_RATES = ('none', 'jaumann')

#: Elastic models whose parameters are written at each point and step by a model of another family
#: (D-017.13d): ``{elas model: family prefix}``. The material has a component of that family in place
#: of ``elas`` and a state with the state of that family (``material_state``): the frozen elasticity of
#: a step comes from the viscoelastic model and its converged state (`visco_freeze`).
MATERIAL_ELAS_PER_POINT = MappingProxyType({'frozen': 'visco'})

#: A return mapping solver usable by a material: ``families`` (prefixes of the model families it
#: calls, in the order of the components of the material type; each one is required), ``nstate``
#: (number of state variables of a point in an FE history array, without the stress, D-018.7),
#: ``unsupported`` (``{prefix: models}`` registered in the family but not yet usable by the solver),
#: ``composes`` (``{prefix: composed prefix}``: a model of family ``prefix`` that composes another law,
#: ``ModelFamily.composed``, takes it from the component of the composed family, which the material
#: then has, right after the composer; with a model that composes nothing the material has no such
#: component, D-017.6b), ``dt`` (the solver uses the time step: rate-dependent; every point kernel
#: takes it, D-018.12c) and ``iter_nw`` (default iteration limit of the material, component
#: ``iter_nw``: the historical limit of the solver). The defaults are immutable (``MappingProxyType``).
#: Its body (``src/solvers/muscle_<solver>_body.fypp``) is included by
#: ``src/materials/muscle_material.fypp`` with ``k`` fixed.
MaterialSolver = namedtuple('MaterialSolver', 'families nstate unsupported composes dt iter_nw',
                            defaults=(MappingProxyType({}), MappingProxyType({}), False, 200))

#: Return mapping solvers of the declared materials (D-018.1). The plastic state of all of them is the
#: plastic strain (6) and the equivalent plastic strain; a material with the frozen viscoelastic
#: elasticity adds its viscoelastic state (``material_state_size``).
MATERIAL_SOLVERS = {
    'closest_point': MaterialSolver(families=('elas', 'yield', 'hard'), nstate=7),
    # Overstress viscoplasticity (Perzyna, Chaboche): every overstress law composes a hardening law
    'closest_point_perzyna': MaterialSolver(families=('elas', 'yield', 'os'), nstate=7, composes={'os': 'hard'},
                                            dt=True),
    # Consistency viscoplasticity (flow stress sigma(ep, epd)): Johnson-Cook composes a hardening law
    'closest_point_time': MaterialSolver(families=('elas', 'yield', 'vp'), nstate=7, composes={'vp': 'hard'},
                                         dt=True),
    # Monolithic Newton on (plastic strain, plastic multiplier): a 7x7 system solved by an own LU (D-018.1)
    'closest_point_2': MaterialSolver(families=('elas', 'yield', 'hard'), nstate=7),
    # Explicit cutting plane (D-018.1): 100 corrections at most, as the historical solver
    'cutting_plane': MaterialSolver(families=('elas', 'yield', 'hard'), nstate=7, iter_nw=100),
}

#: Procedures of a material module, ``mat_<proc>_<name>_<k>`` (``material_proc``): ``update`` is
#: the point kernel (public generic ``mat_update``, shared by every material module), ``pack`` /
#: ``unpack`` move the state to / from an FE history array (public generics ``mat_pack_<name>`` /
#: ``mat_unpack_<name>``: their arguments do not involve the material type), the others are private
#: pieces of the solver (``solve``: the return mapping with the elasticity of the step and the plastic
#: state, called by ``update``; ``merit``: merit function of a line search). ``step`` (public generic
#: ``mat_step``, only with a viscoelastic elasticity): coefficients of a time step, once per step for
#: every point, an argument of ``update``. ``batch``: the batch API over SoA arrays (public generic
#: ``material_batch``, ``muscle_update_<name>``), which runs ``update`` on every point of a batch.
#: ``point`` / ``pointdd``: the point API over flat arrays for FE codes without / with the tangent
#: ``ddsdde`` (specifics of the public generic ``material_point``, ``muscle_point_<name>``, D-018.17b; two
#: specifics instead of an optional argument, D-018.16b); ``rotate``: rotation of the tensorial parts of an FE history
#: array (public generic ``material_rotate_hsv``, ``mat_rotate_hsv_<name>``, D-018.17c; only for a
#: material whose input is a strain, ``material_strain_input``); ``init``: the material from a flat parameter
#: array (public generic ``mat_init_<name>``, D-018.17a; ``material_props``).
MATERIAL_PROCS = ('update', 'solve', 'step', 'iter', 'tangent', 'merit', 'pack', 'unpack', 'batch', 'batchdev',
                  'point', 'pointdd', 'rotate', 'init')

#: Maximum length of a Fortran identifier (F2008).
FORTRAN_NAME_MAX = 63

#: Maximum length of ``Material.doc`` (one line of the generated FORD docs).
MATERIAL_DOC_MAX = 80

_MATERIAL_NAME_RE = r'[a-z][a-z0-9_]*'


def material_solver_body(solver):
    """Body of a solver of ``MATERIAL_SOLVERS`` (D-018.5), relative to ``src/solvers/``, included by
    ``src/materials/muscle_material.fypp``: ``material_solver_body('closest_point')`` ->
    ``muscle_closest_point_body.fypp``. The CMake dependencies of the material modules come from it
    (``muscle_meta.py solver-bodies``)."""
    if solver not in MATERIAL_SOLVERS:
        raise ValueError(f'unknown solver {solver!r} (MATERIAL_SOLVERS: {", ".join(MATERIAL_SOLVERS)})')
    return f'muscle_{solver}_body.fypp'


def material_module(name):
    """Module of a material: ``material_module('cp_vm_swift')`` -> ``muscle_mat_cp_vm_swift``."""
    return f'muscle_mat_{name}'


def material_type(name, kind):
    """Type of a material: ``material_type('cp_vm_swift', 'dp')`` -> ``mat_cp_vm_swift_dp``."""
    return f'mat_{name}_{kind}'


def material_nstate(name):
    """Constant with the size of the state of a material in an FE history array:
    ``material_nstate('cp_vm_swift')`` -> ``MAT_CP_VM_SWIFT_NSTATE``."""
    return f'MAT_{name.upper()}_NSTATE'


def material_visco(mat):
    """Family prefix of the model that writes the elasticity of the material at each point
    (``MATERIAL_ELAS_PER_POINT``: ``'visco'`` for ``elas='frozen'``), or ``None`` (a material
    elasticity with constant parameters, a component of the material type)."""
    return MATERIAL_ELAS_PER_POINT.get(mat.elas)


def material_state_size(mat):
    """Number of state variables of a point of the material in an FE history array (value of
    ``material_nstate``): the plastic state of the solver (``MaterialSolver.nstate``) and, with a
    viscoelastic elasticity, the elastic strain (6) and ``visco_branches`` history tensors (6 each)."""
    n = MATERIAL_SOLVERS[mat.solver].nstate
    if material_visco(mat) is not None:
        n += 6 + 6*mat.visco_branches
    return n


def material_state(mat, kind):
    """``(module, type)`` of the state of a point of the material (argument of ``mat_update``):
    ``plastic_state_<k>`` of muscle_plastic_state, or with a viscoelastic elasticity the type of the
    material module (D-019.7b) ``mat_state_<name>_<k>``, components ``plastic`` and ``visco``
    (``material_visco_state``)."""
    if material_visco(mat) is not None:
        return material_module(mat.name), material_state_type(mat.name, kind)
    return 'muscle_plastic_state', f'plastic_state_{kind}'


def material_state_type(name, kind):
    """State type of a point of a material with a viscoelastic elasticity (D-019.7b), defined in the
    material module: ``mat_state_<name>_<k>``."""
    return f'mat_state_{name}_{kind}'


def material_visco_state(name, kind):
    """Viscoelastic state of a point of a material with a viscoelastic elasticity (D-019.7b): the
    elastic strain and ``visco_branches`` history tensors, ``mat_vstate_<name>_<k>``."""
    return f'mat_vstate_{name}_{kind}'


def material_spatial_kin(mat):
    """Whether the material has a kinematics with a spatial strain measure (``KINEMATICS_SPATIAL``): its
    kernel rotates the converged plastic strain with the rotation increment (D-018.15f)."""
    return mat.kin in KINEMATICS_SPATIAL


def material_kin_rotation(mat):
    """Whether the point kernel of a finite strain material reads the polar rotation R of its kinematic
    states (D-019.7c): a spatial strain measure (rotation increment, D-018.15f) or a kinematics of
    ``KINEMATICS_ROTATION``; otherwise it calls ``kin_update_norot`` (no polar decomposition)."""
    return material_spatial_kin(mat) or mat.kin in KINEMATICS_ROTATION


def material_rotated(mat):
    """Whether the return mapping of the material runs in the frame of the polar rotation: a spatial
    strain measure with an anisotropic yield criterion (D-018.15f)."""
    return material_spatial_kin(mat) and mat.yield_ is not None and mat.yield_ not in YIELD_ISOTROPIC


def material_proc(proc, name, kind):
    """Specific procedure of a material (generic ``mat_<proc>`` + material + kind, D-011.4):
    ``material_proc('update', 'cp_vm_swift', 'dp')`` -> ``mat_update_cp_vm_swift_dp``."""
    if proc not in MATERIAL_PROCS:
        raise ValueError(f'unknown material procedure {proc!r} (MATERIAL_PROCS: {MATERIAL_PROCS})')
    return f'mat_{proc}_{name}_{kind}'


def material_batch(name):
    """Public generic of the batch API of a material (every kind, specifics ``material_proc('batch', ...)``):
    ``material_batch('cp_vm_swift')`` -> ``muscle_update_cp_vm_swift``. Its length is that of
    ``material_proc('unpack', name, k)``, so ``material_identifiers`` already bounds it."""
    return f'muscle_update_{name}'


def material_batch_dev(name):
    """Public generic of the batch API of a material on arrays given by their device addresses (D-019.7a,
    every kind, specifics ``material_proc('batchdev', ...)``): ``material_batch_dev('cp_vm_swift')`` ->
    ``muscle_update_cp_vm_swift_dev``. The longest generated name: listed in ``material_identifiers``."""
    return material_batch(name) + '_dev'


def material_point(name):
    """Public generic of the point API of a material for FE codes (D-018.17b, every kind, specifics
    ``material_proc('point', ...)``): ``material_point('cp_vm_swift')`` -> ``muscle_point_cp_vm_swift``.
    Shorter than ``material_proc('tangent', name, k)``, so ``material_identifiers`` already bounds it."""
    return f'muscle_point_{name}'


def material_rotate_hsv(name):
    """Public generic of the rotation of the tensorial parts of an FE history array (D-018.17c, every
    kind, specifics ``material_proc('rotate', ...)``): ``material_rotate_hsv('cp_vm_swift')`` ->
    ``mat_rotate_hsv_cp_vm_swift``. As long as ``material_proc('tangent', name, k)``, so
    ``material_identifiers`` already bounds it."""
    return f'mat_rotate_hsv_{name}'


def material_strain_input(mat):
    """Whether the input of the point kernel of the material is a strain (small strain, or the
    corotational kinematics: the strain rotated by the FE code, D-018.15f) rather than the deformation
    gradients F_n, F_np1 of a finite strain kinematics. Only such a material has
    ``material_rotate_hsv``: with F the kernel itself keeps the state objective."""
    return mat.kin is None or mat.kin in KINEMATICS_FE_ROTATED


def material_identifiers(name):
    """Every Fortran identifier generated for a material name (module, types, procedures) in every
    kind of ``KINDS``: what the name length is validated against."""
    ids = [material_module(name), material_nstate(name), material_batch_dev(name)]
    for k in KINDS:
        ids += [material_type(name, k), material_state_type(name, k), material_visco_state(name, k)]
        ids += [material_proc(p, name, k) for p in MATERIAL_PROCS]
    return ids


def material_composed(mat, prefix):
    """Prefix of the family whose law the model of family ``prefix`` of the material composes
    (``MaterialSolver.composes``), or ``None``: ``material_composed(Material('a',
    'closest_point_time', 'linear', 'vonmises', 'swift', vp='jc'), 'vp')`` -> ``'hard'`` (``None``
    with ``vp='rk'``)."""
    target = MATERIAL_SOLVERS[mat.solver].composes.get(prefix)
    if target is None or not model_composed(prefix, getattr(mat, MATERIAL_FIELDS[prefix])):
        return None
    return target


def material_families(mat):
    """Model families of the material (prefixes, in the order of the components of its type): the
    solver families, each composer followed by the family it composes when its model composes one;
    ``elas`` is replaced by the family that writes it at each point (``MATERIAL_ELAS_PER_POINT``)."""
    res = []
    for p in MATERIAL_SOLVERS[mat.solver].families:
        if p == 'elas' and material_visco(mat) is not None:
            res.append(material_visco(mat))
            continue
        res.append(p)
        composed = material_composed(mat, p)
        if composed is not None:
            res.append(composed)
    if mat.kin is not None:
        res.append('kin')
    return res


def material_components(mat):
    """Components of the material type, ``[(prefix, model)]`` (``material_families``):
    ``material_components(Material('a', 'closest_point', 'linear', 'vonmises', 'swift'))`` ->
    ``[('elas', 'linear'), ('yield', 'vonmises'), ('hard', 'swift')]``. The component of family
    ``prefix`` is named ``prefix`` and has the parameter type ``model_type(prefix, model, k)``."""
    return [(p, getattr(mat, MATERIAL_FIELDS[p])) for p in material_families(mat)]


def material_comp_args(mat, prefix, func, ep, obj='mat'):
    """Trailing arguments of a call of the family function ``<prefix>_<func>`` of the material:
    the values of the composed law at ``ep`` (Fortran expression) that the function receives
    (``ModelFunc.comp``, D-017.6b), evaluated from the component of the composed family of ``obj``,
    or ``''`` if the model composes nothing: ``material_comp_args(m, 'vp', 'dstress_dep', 'ep')`` ->
    ``', hard_stress(mat%hard, ep), hard_dstress_dep(mat%hard, ep)'`` for Johnson-Cook."""
    composed = material_composed(mat, prefix)
    if composed is None:
        return ''
    return ''.join(f', {composed}_{v.func}({obj}%{composed}, {ep})' for v in model_comp_values(prefix, func))


def _article(word):
    """Indefinite article of a word in a message: ``_article('os')`` -> ``'an'``."""
    return 'an' if word[:1] in 'aeiou' else 'a'


def _validate_material_model(mat, prefix, used):
    """Checks the model of family ``prefix`` of a material (``validate_material``): required, registered
    and supported by the solver if the material uses the family, ``None`` otherwise."""
    solver = MATERIAL_SOLVERS[mat.solver]
    field = MATERIAL_FIELDS[prefix]
    model = getattr(mat, field)
    if not used:
        if model is not None:
            composer = [p for p, c in solver.composes.items() if c == prefix]
            why = (f' (its {composer[0]!r} model {getattr(mat, MATERIAL_FIELDS[composer[0]])!r} composes none)'
                   if composer else '')
            raise ValueError(f'material {mat.name!r}: solver {mat.solver!r} does not use {_article(prefix)} '
                             f'{prefix!r} model ({field}={model!r}){why}')
        return
    if model is None:
        raise ValueError(f'material {mat.name!r}: solver {mat.solver!r} needs {_article(prefix)} {prefix!r} model '
                         f'({field}=...)')
    if model not in MODEL_FAMILIES[prefix].models:
        raise ValueError(f'material {mat.name!r}: {model!r} is not a registered {prefix!r} model '
                         f'(registered: {", ".join(MODEL_FAMILIES[prefix].models)})')
    if model in solver.unsupported.get(prefix, ()):
        raise ValueError(f'material {mat.name!r}: the {prefix!r} model {model!r} is not supported yet '
                         f'by solver {mat.solver!r}')


def validate_material(mat):
    """Checks one ``Material`` against the solver and model family registries; raises
    ``ValueError`` with the reason. Returns ``mat``."""
    if not isinstance(mat, Material):
        raise ValueError(f'{mat!r} is not a muscle_meta.Material')
    if not isinstance(mat.name, str) or not re.fullmatch(_MATERIAL_NAME_RE, mat.name):
        raise ValueError(f'material name {mat.name!r}: a lower case Fortran identifier is required')
    longest = max(material_identifiers(mat.name), key=len)
    if len(longest) > FORTRAN_NAME_MAX:
        raise ValueError(f'material {mat.name!r}: name too long, {longest!r} has {len(longest)} > '
                         f'{FORTRAN_NAME_MAX} characters')
    if not isinstance(mat.doc, str) or len(mat.doc) > MATERIAL_DOC_MAX:
        raise ValueError(f'material {mat.name!r}: doc must be a string of at most {MATERIAL_DOC_MAX} characters')
    if mat.solver not in MATERIAL_SOLVERS:
        raise ValueError(f'material {mat.name!r}: unknown solver {mat.solver!r} (known: '
                         f'{", ".join(MATERIAL_SOLVERS)})')
    solver = MATERIAL_SOLVERS[mat.solver]
    # The solver families first: whether a composed family is used depends on their models
    for prefix in solver.families:
        _validate_material_model(mat, prefix, True)
    families = material_families(mat)
    for prefix in MATERIAL_FIELDS:
        if prefix not in solver.families:
            _validate_material_model(mat, prefix, prefix in families)
    nb = mat.visco_branches
    if material_visco(mat) is None:
        if nb != 0:
            raise ValueError(f'material {mat.name!r}: visco_branches={nb!r} without a viscoelastic elasticity')
    elif not isinstance(nb, int) or isinstance(nb, bool) or not 0 <= nb <= MAX_VISCO_BRANCHES:
        raise ValueError(f'material {mat.name!r}: visco_branches must be an integer in [0, {MAX_VISCO_BRANCHES}] '
                         f'(got {nb!r})')
    if mat.tangent_rate not in MATERIAL_TANGENT_RATES:
        raise ValueError(f'material {mat.name!r}: tangent_rate must be one of {", ".join(MATERIAL_TANGENT_RATES)} '
                         f'(got {mat.tangent_rate!r})')
    if mat.kin is None and mat.tangent_rate != 'none':
        raise ValueError(f'material {mat.name!r}: tangent_rate={mat.tangent_rate!r} without a kinematics (kin=...)')
    if mat.kin is not None and material_visco(mat) is not None:
        raise ValueError(f'material {mat.name!r}: a kinematics with a viscoelastic elasticity is not supported yet')
    if mat.kin is not None and mat.solver == 'cutting_plane':
        raise ValueError(f'material {mat.name!r}: a kinematics needs the consistent tangent of the return mapping '
                         f'(cutting_plane returns the continuum tangent)')
    if material_rotated(mat) and mat.elas not in ELASTICITY_ISOTROPIC:
        raise ValueError(f'material {mat.name!r}: a spatial strain measure with an anisotropic yield criterion '
                         f'needs an isotropic elasticity ({", ".join(ELASTICITY_ISOTROPIC)})')
    material_props(mat)
    return mat


def validate_materials(materials):
    """Checks a registry (sequence of ``Material``, unique names, at least one) and returns it as a
    tuple."""
    materials = tuple(materials)
    if not materials:
        raise ValueError('the material registry is empty')
    seen = set()
    for m in materials:
        validate_material(m)
        if m.name in seen:
            raise ValueError(f'material {m.name!r} declared twice')
        seen.add(m.name)
    return materials


def load_python_file(path, name=None, base=None):
    """Executes the Python file ``path`` as a module and returns it (data files of the templates: the
    material registry, the parameter sets of the material tests). A relative ``path`` is relative to
    the directory of the file ``base`` (e.g. fypp's ``_THIS_FILE_``) or, without it, to the current
    directory; ``name``: module name (default ``_muscle_data_<file stem>``)."""
    if base is not None and not os.path.isabs(path):
        path = os.path.join(os.path.dirname(os.path.abspath(base)), path)
    if name is None:
        name = '_muscle_data_' + os.path.splitext(os.path.basename(path))[0]
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None:
        raise ValueError(f'cannot load the Python file {path!r}')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_materials(path=None):
    """The validated material registry: the ``MATERIALS`` of the Python file ``path``
    (``-DMUSCLE_MATERIALS=<path>``), or of ``tools/fypp/muscle_materials.py`` (the test registry)
    if ``path`` is empty or ``None``. The file builds its entries with ``muscle_meta.Material``."""
    if not path:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'muscle_materials.py')
    module = load_python_file(path, '_muscle_material_registry')
    if not hasattr(module, 'MATERIALS'):
        raise ValueError(f'material registry {path!r} defines no MATERIALS')
    return validate_materials(module.MATERIALS)


def material(name, path=None):
    """The ``Material`` ``name`` of the registry ``path`` (see ``load_materials``)."""
    for m in load_materials(path):
        if m.name == name:
            return m
    raise ValueError(f'material {name!r} is not in the registry {path or "muscle_materials.py"!r}')


#: Parameter of the flat parameter array of a material (``material_props``, D-018.17a, e): ``index``
#: (1-based position in ``props``), ``component`` (prefix of the material component, ``None`` for a
#: setting of the solver), ``name`` (the input of the model or the component of the material, as in
#: Fortran), ``branch`` (branch of an input with one value per viscoelastic branch, else 0), ``type``
#: (``'real'``, or ``'integer'``: given as a real and rounded with ``nint``), ``doc`` and ``domain`` (the
#: ``ModelInput.domain`` checked by ``mat_init_<name>``, or ``None``).
MaterialProp = namedtuple('MaterialProp', 'index component name branch type doc domain')

#: Settings of the solver at the end of the parameter array of every material (components of the material
#: type), both > 0 (D-018.17e; the defaults are ``NEWTON_TOL_<k>`` and ``MaterialSolver.iter_nw``).
MATERIAL_SOLVER_PROPS = (
    ModelInput('tol', 'real', '', r'Relative tolerance of the Newton iterations (> 0; default of the type: '
               r'`NEWTON_TOL_<k>`).', '{x} > 0.0_{k}'),
    ModelInput('iter_nw', 'integer', '', r'Maximum number of iterations (> 0).', '{x} > 0'))


def material_nprops(name):
    """Constant with the size of the parameter array of a material: ``material_nprops('cp_vm_swift')`` ->
    ``MAT_CP_VM_SWIFT_NPROPS``."""
    return f'MAT_{name.upper()}_NPROPS'


def material_props(mat):
    """Parameter array of a material (``MaterialProp``, D-018.17e): the inputs of the model of each
    component (``model_inputs``, components in the order of ``material_components``) and then
    ``MATERIAL_SOLVER_PROPS``. An input with one value per branch (``'(:)'``) takes ``visco_branches``
    values, interleaved by branch with the other such inputs of the model (Maxwell: bulk, shear_inf,
    shear_1, tau_1, shear_2, tau_2, ...). Only real and integer inputs (an integer is rounded with
    ``nint``): ``[p.name for p in material_props(Material('a', 'closest_point', 'linear', 'vonmises',
    'bilinear'))]`` -> ``['young', 'poisson', 'y0', 'K', 'tol', 'iter_nw']``."""
    res = []

    def add(component, inp, branch=0):
        if inp.type not in ('real', 'integer'):
            raise ValueError(f'material {mat.name!r}: input {inp.name!r} of {component!r} has type {inp.type!r}: '
                             f'only real and integer inputs can be given in a parameter array')
        res.append(MaterialProp(len(res) + 1, component, inp.name, branch, inp.type, inp.doc, inp.domain))

    for p, m in material_components(mat):
        inputs = model_inputs(p, m)
        for inp in inputs:
            if inp.shape == '':
                add(p, inp)
            elif inp.shape != '(:)':
                raise ValueError(f'material {mat.name!r}: input {inp.name!r} of {p!r} has the fixed shape '
                                 f'{inp.shape!r}: not supported in a parameter array')
        for b in range(1, mat.visco_branches + 1):
            for inp in inputs:
                if inp.shape == '(:)':
                    add(p, inp, b)
    for inp in MATERIAL_SOLVER_PROPS:
        add(None, inp)
    return tuple(res)


def material_prop_value(prop, props='props'):
    """Fortran value of a parameter of the array: ``props(i)``, or ``nint(props(i))`` for an integer."""
    ref = f'{props}({prop.index})'
    return f'nint({ref})' if prop.type == 'integer' else ref


def material_prop_check(prop, kind, props='props'):
    """Fortran condition that the value of a parameter of the array is valid (``MaterialProp.domain``), or
    ``None``: ``material_prop_check(MaterialProp(5, None, 'iter_nw', 0, 'integer', '', '{x} > 0'), 'dp')``
    -> ``'nint(props(5)) > 0'``."""
    if prop.domain is None:
        return None
    return prop.domain.format(x=material_prop_value(prop, props), k=kind)


#: Conditions on several inputs of a model that ``mat_init_<name>`` checks (D-018.19c),
#: ``MODEL_DOMAINS[(prefix, model)]``: ``(condition, doc)`` pairs, the condition a Fortran logical expression
#: with each input as ``{name}`` and literals of kind ``{k}``, ``doc`` its FORD text;
#: ``STATUS_INVALID_PARAMETER`` otherwise (as ``ModelInput.domain``: the model initializer and the kernel
#: do not check them). Johnson-Cook: the rate factor \( 1 + C\ln x \) stays positive on the rate floor
#: \( x = \dot\epsilon_{min}/\dot\epsilon_{p0} \) (else a negative flow stress there).
MODEL_DOMAINS = MappingProxyType({
    ('vp', 'jc'): (('{C}*abs(log({epmin}/{epdmax})) < 1.0_{k}',
                    r'\( C\,|\ln(\dot\epsilon_{min}/\dot\epsilon_{p0})| < 1 \) (positive flow stress on the rate '
                    r'floor)'),),
})


def material_domain_checks(mat, kind, props='props'):
    """``(condition, doc)`` of every ``MODEL_DOMAINS`` entry of the components of a material, the condition
    with the values of the parameter array (``material_prop_value``):
    ``[c for c, _ in material_domain_checks(Material('a', 'closest_point_time', 'linear', 'vonmises',
    'swift', vp='jc'), 'dp')]`` -> ``['props(3)*abs(log(props(5)/props(4))) < 1.0_dp']``."""
    res = []
    for p, m in material_components(mat):
        mine = {q.name: q for q in material_props(mat) if q.component == p}
        for cond, doc in MODEL_DOMAINS.get((p, m), ()):
            names = re.findall(r'\{(\w+)\}', cond)
            missing = [n for n in names if n != 'k' and n not in mine]
            if missing:
                raise ValueError(f'MODEL_DOMAINS of {p} {m}: {missing} are not inputs of the model')
            values = {n: material_prop_value(q, props) for n, q in mine.items()}
            res.append((cond.format(k=kind, **values), doc))
    return res


def material_init_call(mat, prefix, kind, obj='mat', props='props', indent=8, width=120):
    """Call of the initializer of the model of component ``prefix`` (``MODEL_INITS``) with keyword
    arguments from the parameter array: an input per branch as a strided section (``shear=props(3:5:2)``),
    then ``status=status`` if the initializer has one; continued with ``&`` before ``width`` columns."""
    model = getattr(mat, MATERIAL_FIELDS[prefix])
    init = model_init(prefix, model)
    mine = [q for q in material_props(mat) if q.component == prefix]
    args = []
    for inp in init.inputs:
        qs = [q for q in mine if q.name == inp.name]
        if inp.shape == '(:)':
            step = len([i for i in init.inputs if i.shape == '(:)'])
            first = qs[0].index if qs else 1
            last = qs[-1].index if qs else 0
            args.append(f'{inp.name}={props}({first}:{last}:{step})')
        else:
            args.append(f'{inp.name}={material_prop_value(qs[0], props)}')
    if init.status:
        args.append('status=status')
    head = ' ' * indent + f'call {prefix}_{model}_init({obj}%{prefix}, '
    lines, cur = [], head
    for n, a in enumerate(args):
        piece = a + (', ' if n < len(args) - 1 else ')')
        if len(cur) + len(piece.rstrip()) + 2 > width:
            lines.append(cur.rstrip() + ' &')
            cur = ' ' * (indent + 4) + piece
        else:
            cur += piece
    lines.append(cur.rstrip())
    return '\n'.join(lines)


def material_props_table(mat):
    """Markdown table of the parameter array of a material (``material_props``): position, component
    (model), parameter and description; printed by ``muscle_meta.py props <name> [PATH]``."""
    models = dict(material_components(mat))
    rows = [f'Parameters of material `{mat.name}` (`{material_nprops(mat.name)}` = {len(material_props(mat))}, '
            f'`mat_init_{mat.name}`):', '',
            '| # | Component | Parameter | Description |', '|---|-----------|-----------|-------------|']
    for q in material_props(mat):
        comp = f'{q.component} ({model_label(q.component, models[q.component])})' if q.component else 'solver'
        name = q.name + (f'({q.branch})' if q.branch else '')
        doc = q.doc.replace('|', r'\|') + (' Integer, given as a real.' if q.type == 'integer' else '')
        rows.append(f'| {q.index} | {comp} | `{name}` | {doc} |')
    for _, doc in material_domain_checks(mat, 'dp'):
        rows += ['', f'`mat_init_{mat.name}` also requires {doc} (D-018.19c).']
    return '\n'.join(rows)


def _main(argv):
    """Command line used by CMake: ``muscle_meta.py materials [PATH]`` prints the names of the
    validated material registry separated by ``;`` (a CMake list); ``muscle_meta.py props NAME [PATH]``
    prints the Markdown table of the parameter array of a material (``material_props_table``, for the
    documentation of an FE adapter); ``muscle_meta.py solver-bodies`` prints the body of each solver
    (``material_solver_body``, a CMake list). Run as a script this file is
    ``__main__``: the registry is loaded through the module ``muscle_meta`` (the one whose
    ``Material`` the registry files import)."""
    import muscle_meta
    if len(argv) in (1, 2) and argv[0] == 'materials':
        print(';'.join(m.name for m in muscle_meta.load_materials(argv[1] if len(argv) == 2 else None)))
        return 0
    if argv == ['solver-bodies']:
        print(';'.join(muscle_meta.material_solver_body(s) for s in muscle_meta.MATERIAL_SOLVERS))
        return 0
    if len(argv) in (2, 3) and argv[0] == 'props':
        print(muscle_meta.material_props_table(muscle_meta.material(argv[1], argv[2] if len(argv) == 3 else None)))
        return 0
    print('usage: muscle_meta.py materials [PATH] | props NAME [PATH] | solver-bodies', file=sys.stderr)
    return 2


if __name__ == '__main__':
    try:
        sys.exit(_main(sys.argv[1:]))
    except ValueError as err:
        print(f'muscle_meta: {err}', file=sys.stderr)
        sys.exit(1)
