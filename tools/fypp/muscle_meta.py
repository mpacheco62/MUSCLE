# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Python helpers and tensor layout metadata for the MUSCLE FYPP templates.

fypp loads this module with ``-M tools/fypp -m muscle_meta``; ``include/muscle_common.fypp``
re-exports the helpers under short names. Rules: ``.claude/rules/dims-2d-3d.md``.

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

import itertools
from collections import namedtuple

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
#: ``src/plasticity/hardening_laws/muscle_hard_<name>.fypp`` each (model contract D-017.6,
#: include/muscle_models_common.fypp), exported by the family module ``muscle_hard_laws``.
HARDENING_MODELS = ('bilinear', 'ludwik', 'swift', 'vocemod')

#: Registered viscoplastic flow stress laws \( \sigma(\epsilon_p, \dot\epsilon_p) \) (D-017): one
#: template ``src/plasticity/viscoplastic_laws/muscle_vp_<name>.fypp`` each, exported by
#: ``muscle_vp_laws``. ``VISCOPLASTIC_COMPOSED``: those that compose a hardening law (D-017.6b).
VISCOPLASTIC_MODELS = ('jc', 'jcfull', 'mrk', 'nnl', 'rk', 'va')
VISCOPLASTIC_COMPOSED = ('jc', 'jcfull')

#: Registered overstress laws \( \dot\lambda = g(\sigma_{eq}, \epsilon_p) \) (D-017): same
#: directory and module names as the flow stress laws, exported by ``muscle_vp_overstress_laws``.
OVERSTRESS_MODELS = ('chaboche', 'perzyna')
OVERSTRESS_COMPOSED = ('chaboche', 'perzyna')

#: Registered yield criteria \( \sigma_{eq}(\boldsymbol\sigma) \) (D-017): one template
#: ``src/plasticity/yield_criteria/muscle_yield_<name>.fypp`` each, exported by
#: ``muscle_yield_criteria``. ``YIELD_JOINT``: those that write their own joint evaluation of the
#: equivalent stress and its derivatives (the others get the default, three calls in turn).
YIELD_MODELS = ('vonmises', 'hill48', 'druckerprager', 'cpb06', 'cpb06ex2', 'yld2004')
YIELD_JOINT = ('yld2004',)

#: Registered elastic models (D-017.12): one template
#: ``src/elasticity/muscle_elasticity_<name>.fypp`` each, exported by ``muscle_elasticity_laws``.
ELASTICITY_MODELS = ('linear', 'frozen')

#: Registered linear viscoelastic models (D-017.13): one template
#: ``src/viscoelasticity/muscle_visco_<name>.fypp`` each, exported by ``muscle_visco_laws``. Their
#: step map is affine in the elastic strain: the algorithmic ("frozen") elasticity of the step is
#: the elastic model ``frozen`` (``elas_frozen_<k>``), which a return mapping uses as any other.
VISCOELASTIC_MODELS = ('maxwell', 'kelvin')

#: Registered kinematics (D-017.14): one template ``src/kinematics/muscle_kin_<name>.fypp`` each,
#: exported by ``muscle_kin_laws``. They have no material parameters: the parameter type of a
#: kinematics is empty and only selects it at compile time (generic resolution, D-008 renames).
KINEMATICS_MODELS = ('small_strain', 'total_lagrangian', 'updated_lagrangian', 'material_logarithmic',
                     'spatial_logarithmic', 'corotational')
#: Kinematics updated from F (``kin_update``); the others (the corotational one) from a strain and a
#: rotation (``kin_corot_update``, D-017.14d).
KINEMATICS_F_MODELS = tuple(m for m in KINEMATICS_MODELS if m != 'corotational')

#: Module of a model whose name does not follow ``ModelFamily.module`` (a model that keeps its
#: historical module, ``model_module``): the frozen elasticity of the viscoelastic models.
MODEL_MODULES = {('elas', 'frozen'): 'muscle_visco_frozen'}

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
# A family fixes the contract every model of it follows: the name of the parameter argument, the
# arguments and result of each function (same names for every model, D-007), the step scale and
# domain of the default finite differences, and the registered models. The macros of
# include/muscle_models_common.fypp generate the repetitive code from it.

#: Argument (or result) of a family function: ``name``; ``type`` = ``'real'``, ``'integer'``,
#: ``'logical'`` or a tensor family (``'2Osym'``, ...; needs the dimension); ``shape`` = ``''``
#: (scalar) or a fixed shape such as ``'(3)'``; ``doc`` = FORD text.
ModelArg = namedtuple('ModelArg', 'name type shape doc')

#: Function of a family: ``name`` (generic ``<prefix>_<name>``), ``doc`` (FORD summary, no final
#: period), ``result`` (a ``ModelArg``), ``deriv`` = ``(function, argument)`` this one is the
#: derivative of, or ``None``; ``extra`` (``ModelArg`` after the family arguments, only in this
#: function, e.g. the time step), ``comp`` (names of the composition values it receives when the
#: model composes another law, D-017.6b, after ``extra``) and ``default`` (formula of the default
#: body in terms of the other functions, ``{name}`` = call of function ``name`` of the same model,
#: generated by MODEL_DEFAULT; or ``None``).
#: ``outs`` (optional): ``ModelArg`` of the ``intent(out)`` arguments of a procedure that is a
#: subroutine (e.g. results with a ``status``, or a new state): then ``result`` is ``None`` and the
#: outputs follow the other arguments. A subroutine has no finite difference nor default body.
#: ``models`` (optional): the registered models that have this function (``None``: every model),
#: for a family whose models do not all share it (the corotational kinematics is updated from a
#: strain and a rotation, the others from F); ``model_procs`` of a model lists only its own.
#: ``opts`` (optional): ``ModelArg`` of ``intent(in), optional`` arguments, last (after ``outs``).
ModelFunc = namedtuple('ModelFunc', 'name doc result deriv extra comp default outs models opts',
                       defaults=((), (), None, (), None, ()))

#: Parameter of a model: ``name``, ``type`` and ``shape`` as in ``ModelArg``, ``doc`` (FORD text)
#: and ``default`` (``None`` or a number written without kind, e.g. ``'0.0'``; the kind is added
#: for real parameters).
#: ``private`` (optional, default ``False``): the component is private to the model module (derived
#: parameters such as the normalization factor of CPB06, D-017.10b): the type can then only be built
#: with the initializer of the model (``yield_<m>_init``); the historical OO type keeps it public.
ModelParam = namedtuple('ModelParam', 'name type shape doc default private', defaults=(False,))

#: Value of a composed law (D-017.6b): the ``arg`` (``ModelArg``) the free functions receive, and
#: the type-bound function of the composed OO law (``binding``) that the OO adapter evaluates at
#: the family argument ``ModelComposition.arg`` to obtain it.
CompValue = namedtuple('CompValue', 'arg binding')

#: Composition of a family whose models may contain another law (viscoplastic laws with a hardening
#: law): ``component`` (name of the polymorphic component of the OO adapter), ``base_type`` (its
#: declared type), ``arg`` (family argument the composed values depend on: a finite difference
#: in that argument cannot be taken with the values frozen) and ``values`` (``CompValue``).
ModelComposition = namedtuple('ModelComposition', 'component base_type arg values')

#: Joint evaluation of several functions of a family at the same arguments (yield criteria: the
#: equivalent stress and its two derivatives, which a return mapping needs together): a
#: subroutine ``<prefix>_<name>`` with the family arguments and one ``intent(out)`` argument per
#: function of ``funcs`` (named ``outs``, of the type of the function result), and a ``flag``
#: function ``<prefix>_<flag>(par)`` (constant ``logical``: whether the model writes its own,
#: cheaper joint evaluation; the others get the default, the functions called in turn).
#: ``models``: the registered models that write their own. ``oo_name`` / ``oo_flag``: bindings of
#: the historical OO base (overridden by the adapter only for ``models``); ``oo_args``: OO dummy
#: names of the family arguments.
#: ``optional``: outputs that are ``optional`` arguments (not computed when absent, D-017.10a).
ModelJoint = namedtuple('ModelJoint', 'name doc funcs outs flag flag_doc models oo_name oo_flag oo_args optional')

#: A model family: ``prefix`` of every generated name (``hard`` -> ``hard_stress``,
#: ``hard_bilinear_dp``), ``module`` (format of the model modules), ``family_module``, ``par``
#: (name of the parameter argument: chosen by the family so that it cannot clash with a physical
#: argument such as a pressure ``p``), ``args`` (``ModelArg`` of every function, in order),
#: ``funcs`` (``ModelFunc``), ``fd_scale`` (scale of the argument in the default finite
#: differences, ``h = rel max(|x|, scale)``, D-017.5a, or ``None``), ``fd_lower`` (lower bound of
#: the domain of the differentiated argument, or ``None``: below ``lower + h`` the default finite
#: differences are one-sided), ``models`` (the registry), ``composition`` (a
#: ``ModelComposition`` or ``None``), ``composed`` (registered models that compose the other law)
#: ``oo_funcs`` (functions bound by default in the OO adapter, ``None`` = all), ``oo_bindings``
#: (default OO bindings ``(binding, function, OO dummy names)`` when the historical names differ
#: from the family ones, or ``None``), ``dims`` (dimensions of a family with tensor arguments: the
#: specific function names then carry the dimension, ``<prefix>_<func>_<model>_<dim>_<kind>``;
#: ``None`` for a family of scalar arguments), ``joint`` (a ``ModelJoint`` or ``None``) and
#: ``oo_class`` (the tensor arguments of the historical OO base are polymorphic, ``class(...)``,
#: so the OO adapters declare them so to override its bindings; ``False``: ``type(...)``).
#: A default OO binding of ``oo_bindings`` may carry a fourth element, the dimension of its tensor
#: arguments (default ``OO_DIM``), for a historical API with bindings per dimension
#: (``stress_3D``, ``stress_2D``).
#: ``support``: templates of the family that follow ``module`` but are not models (shared types,
#: macros included by the model templates).
ModelFamily = namedtuple('ModelFamily',
                         'prefix module family_module par args funcs fd_scale fd_lower models '
                         'composition composed oo_funcs oo_bindings dims joint oo_class support',
                         defaults=(None, (), None, None, None, None, False, ()))

#: Dimension of the historical OO API (the OO adapters of a family with ``dims``).
OO_DIM = '3D'

_EP = ModelArg('ep', 'real', '', r'Equivalent plastic strain \( \epsilon_p \ge 0 \).')
_EPD = ModelArg('epd', 'real', '', r'Equivalent plastic strain rate \( \dot\epsilon_p \).')
_DT = ModelArg('dt', 'real', '', r'Time step \( \Delta t \) (\( \dot\epsilon_p = \Delta\epsilon_p/\Delta t \)).')
_SEQ = ModelArg('sigma_eq', 'real', '', r'Equivalent stress \( \sigma_{eq} \).')
_REAL = lambda doc: ModelArg('res', 'real', '', doc)
_SIGMA = ModelArg('sigma', '2Osym', '', r'Stress tensor \( \boldsymbol\sigma \).')
_STRAIN = ModelArg('strain', '2Osym', '', r'Strain tensor \( \boldsymbol\varepsilon \).')
_COEFS = ModelArg('coefs', 'type:visco_coefs', '', r'Coefficients of the time step (`visco_step_coefficients`).')
_STATE_N = ModelArg('state_n', 'type:visco_state', '', r'Converged state at \( t_n \).')
_KIN_STATE = ModelArg('state', 'type:kin_state', '', r'Kinematic state (`kin_update` or `kin_corot_update`).')

#: Composition with an isotropic hardening law (viscoplastic and overstress families): the
#: composer evaluates ``hard_stress`` and ``hard_dstress_dep`` at ``ep`` and passes the values.
_HARD_COMPOSITION = ModelComposition(
    component='hard_law', base_type='Base_hardening_laws', arg='ep',
    values=(CompValue(ModelArg('sy', 'real', '', r'Static yield stress \( \sigma_y(\epsilon_p) \) of the '
                               r'composed hardening law (`hard_stress`).'), 'stress'),
            CompValue(ModelArg('dsy_dep', 'real', '', r'Hardening modulus \( d\sigma_y/d\epsilon_p \) of '
                               r'the composed hardening law (`hard_dstress_dep`).'), 'dstress_dep')))

#: Model families by prefix.
MODEL_FAMILIES = {
    'hard': ModelFamily(
        prefix='hard', module='muscle_hard_{}', family_module='muscle_hard_laws', par='prm',
        args=(_EP,),
        funcs=(ModelFunc('stress', r'Flow stress \( \sigma_y(\epsilon_p) \)',
                         _REAL(r'\( \sigma_y \).'), None),
               ModelFunc('dstress_dep', r'Hardening modulus \( d\sigma_y/d\epsilon_p \)',
                         _REAL(r'\( d\sigma_y/d\epsilon_p \).'), ('stress', 'ep')),
               ModelFunc('ddstress_ddep', r'Second derivative \( d^2\sigma_y/d\epsilon_p^2 \)',
                         _REAL(r'\( d^2\sigma_y/d\epsilon_p^2 \).'),
                         ('dstress_dep', 'ep'))),
        fd_scale='1.0', fd_lower='0.0', models=HARDENING_MODELS),
    # Flow stress laws sigma(ep, epd). The partial derivatives are the kernel API; dstress_dep is
    # the historical derivative used by the return mapping, d sigma(ep_n + x, x/dt)/dx (with the
    # time step), by default dflow_dep + dflow_depd/dt (the formula of the historical laws and of
    # the numeric default of Base_viscoplastic_law). Default finite differences as the historical
    # `derivative` (no scale, no lower bound: the laws are evaluated at epd < 0 too).
    'vp': ModelFamily(
        prefix='vp', module='muscle_vp_{}', family_module='muscle_vp_laws', par='prm',
        args=(_EP, _EPD),
        funcs=(ModelFunc('flow_stress', r'Viscoplastic flow stress \( \sigma(\epsilon_p, \dot\epsilon_p) \)',
                         _REAL(r'\( \sigma \).'), None, comp=('sy',)),
               ModelFunc('dflow_dep', r'Partial derivative \( \partial\sigma/\partial\epsilon_p \) '
                         r'(at constant \( \dot\epsilon_p \))',
                         _REAL(r'\( \partial\sigma/\partial\epsilon_p \).'), ('flow_stress', 'ep'),
                         comp=('sy', 'dsy_dep')),
               ModelFunc('dflow_depd', r'Partial derivative \( \partial\sigma/\partial\dot\epsilon_p \) '
                         r'(at constant \( \epsilon_p \))',
                         _REAL(r'\( \partial\sigma/\partial\dot\epsilon_p \).'), ('flow_stress', 'epd'),
                         comp=('sy',)),
               ModelFunc('dstress_dep', r'Derivative of the flow stress with respect to the plastic '
                         r'strain increment, \( \partial\sigma/\partial\epsilon_p + '
                         r'\partial\sigma/\partial\dot\epsilon_p/\Delta t \)',
                         _REAL(r'\( d\sigma/d\Delta\epsilon_p \).'), None, extra=(_DT,),
                         comp=('sy', 'dsy_dep'), default='{dflow_dep} + {dflow_depd}/dt')),
        fd_scale=None, fd_lower=None, models=VISCOPLASTIC_MODELS,
        composition=_HARD_COMPOSITION, composed=VISCOPLASTIC_COMPOSED,
        oo_funcs=('flow_stress', 'dstress_dep')),
    # Overstress laws g(sigma_eq, ep) (plastic multiplier rate). Default finite differences as the
    # historical `derivative` of Base_overstress_law.
    'os': ModelFamily(
        prefix='os', module='muscle_vp_{}', family_module='muscle_vp_overstress_laws', par='prm',
        args=(_SEQ, _EP),
        funcs=(ModelFunc('flow_rate', r'Plastic multiplier rate \( g(\sigma_{eq}, \epsilon_p) \)',
                         _REAL(r'\( g \).'), None, comp=('sy',)),
               ModelFunc('dflow_dsigma', r'Derivative \( \partial g/\partial\sigma_{eq} \)',
                         _REAL(r'\( \partial g/\partial\sigma_{eq} \).'), ('flow_rate', 'sigma_eq'),
                         comp=('sy',)),
               ModelFunc('dflow_dep', r'Derivative \( \partial g/\partial\epsilon_p \)',
                         _REAL(r'\( \partial g/\partial\epsilon_p \).'), ('flow_rate', 'ep'),
                         comp=('sy', 'dsy_dep'))),
        fd_scale=None, fd_lower=None, models=OVERSTRESS_MODELS,
        composition=_HARD_COMPOSITION, composed=OVERSTRESS_COMPOSED),
    # Yield criteria sigma_eq(sigma) and their derivatives with respect to the stress tensor. Default
    # finite differences as the numeric procedures of Base_yield_critera (FD_GRAD, and FD_HESS of
    # stress_eq for the second derivative; no scale). Only 3D for now (dims-2d-3d.md: 2D needs a
    # 2D vs 3D test first). The OO adapters keep the historical binding names.
    'yield': ModelFamily(
        prefix='yield', module='muscle_yield_{}', family_module='muscle_yield_criteria', par='prm',
        args=(_SIGMA,),
        funcs=(ModelFunc('stress_eq', r'Equivalent stress \( \sigma_{eq}(\boldsymbol\sigma) \)',
                         _REAL(r'\( \sigma_{eq} \).'), None),
               ModelFunc('dstress_eq', r'Gradient \( \partial\sigma_{eq}/\partial\boldsymbol\sigma \)',
                         ModelArg('res', '2Osym', '', r'\( \partial\sigma_{eq}/\partial\boldsymbol\sigma \).'),
                         ('stress_eq', 'sigma')),
               ModelFunc('ddstress_eq', r'Second derivative \( \partial^2\sigma_{eq}/\partial\boldsymbol\sigma^2 \)',
                         ModelArg('res', '4O3sym', '', r'\( \partial^2\sigma_{eq}/\partial\boldsymbol\sigma^2 \).'),
                         ('dstress_eq', 'sigma'))),
        fd_scale=None, fd_lower=None, models=YIELD_MODELS,
        oo_bindings=(('stress_eq', 'stress_eq', ('stress',)),
                     ('dstressEq_dstress', 'dstress_eq', ('stress',)),
                     ('ddstressEq_ddstress', 'ddstress_eq', ('stress',))),
        dims=('3D',),
        joint=ModelJoint(
            name='stress_eq_and_derivatives',
            doc=r'Equivalent stress and its first and second derivatives at the same stress',
            funcs=('stress_eq', 'dstress_eq', 'ddstress_eq'), outs=('f', 'df', 'ddf'),
            flag='joint_evaluation',
            flag_doc=r'Whether the criterion writes its own joint evaluation (cheaper than the three '
                     r'functions in turn)',
            models=YIELD_JOINT, oo_name='stressEq_and_derivatives', oo_flag='joint_evaluation',
            oo_args=('stress',), optional=('ddf',))),
    # Elastic models sigma(eps): the stress, its deviatoric part and their tangents (D-017.12), in
    # 3D and 2D (plane strain/stress and axisymmetry: the 33 component is kept). The tangents of
    # 4O3sym have no default finite differences (tensor Jacobian): every model writes them. The
    # transient variants of the historical OO API (strain increment, time step and internal
    # variables, ignored by every elastic model) are not part of the kernel API: a model with
    # internal variables belongs to the viscoelastic family. The OO base declares its tensor
    # arguments `class(...)` and has one binding per dimension.
    'elas': ModelFamily(
        prefix='elas', module='muscle_elasticity_{}', family_module='muscle_elasticity_laws', par='prm',
        args=(_STRAIN,),
        funcs=(ModelFunc('stress', r'Stress \( \boldsymbol\sigma(\boldsymbol\varepsilon) \)',
                         ModelArg('res', '2Osym', '', r'\( \boldsymbol\sigma \).'), None),
               ModelFunc('stress_dev', r'Deviatoric stress \( \operatorname{dev}\boldsymbol\sigma(\boldsymbol\varepsilon) \)',
                         ModelArg('res', '2Osym', '', r'\( \operatorname{dev}\boldsymbol\sigma \).'), None),
               ModelFunc('dstress_dstrain', r'Tangent \( \partial\boldsymbol\sigma/\partial\boldsymbol\varepsilon \)',
                         ModelArg('res', '4O3sym', '', r'\( \partial\boldsymbol\sigma/\partial\boldsymbol\varepsilon \).'),
                         ('stress', 'strain')),
               ModelFunc('dstress_dstrain_dev', r'Deviatoric tangent '
                         r'\( \partial(\operatorname{dev}\boldsymbol\sigma)/\partial\boldsymbol\varepsilon \)',
                         ModelArg('res', '4O3sym', '',
                                  r'\( \partial(\operatorname{dev}\boldsymbol\sigma)/\partial\boldsymbol\varepsilon \).'),
                         ('stress_dev', 'strain'))),
        fd_scale=None, fd_lower=None, models=ELASTICITY_MODELS,
        oo_bindings=(('stress_3D', 'stress', ('strain',), '3D'),
                     ('stress_2D', 'stress', ('strain',), '2D'),
                     ('stress_dev_3D', 'stress_dev', ('strain',), '3D'),
                     ('stress_dev_2D', 'stress_dev', ('strain',), '2D'),
                     ('dstress_dstrain_3D', 'dstress_dstrain', ('strain',), '3D'),
                     ('dstress_dstrain_2D', 'dstress_dstrain', ('strain',), '2D'),
                     ('dstress_dstrain_dev_3D', 'dstress_dstrain_dev', ('strain',), '3D')),
        dims=('3D', '2D'), oo_class=True),
    # Linear viscoelastic models (D-017.13), integrated over a time step: the first family with a
    # state. The state is explicit (`visco_state_<k>`, fixed size, muscle_visco_types) and the
    # converged state is never modified: `visco_update` reads `state_n` and writes `state_np1`.
    # The work is split by how often it is done (as the historical Base_viscoelasticity):
    # `visco_step_coefficients` once per step and material (it validates the time step: the only
    # procedure with a `status`), `visco_freeze` once per step and point (the frozen elasticity
    # `elas_frozen_<k>` of the step, an elastic model of the `elas` family that any return mapping
    # uses), `visco_update` once per converged step and point. Only 3D (the state and the models
    # are 3D, D-017.13c). The OO adapters keep the historical bindings, written by the author.
    'visco': ModelFamily(
        prefix='visco', module='muscle_visco_{}', family_module='muscle_visco_laws', par='prm',
        args=(),
        funcs=(ModelFunc('step_coefficients', r'Coefficients of a time step (algorithmic moduli and the factors '
                         r'of each branch; the time step is validated)',
                         None, None,
                         extra=(ModelArg('dt', 'real', '', r'Time step \( \Delta t \).'),),
                         outs=(ModelArg('coefs', 'type:visco_coefs', '', r'Coefficients of the step (default '
                                        r'values when `status` is not `STATUS_OK`).'),
                               ModelArg('status', 'integer', '', r'`STATUS_OK` or `STATUS_INVALID_TIME_STEP` '
                                        r'(muscle_status).'))),
               ModelFunc('freeze', r'Frozen (algorithmic) elasticity of a time step at a point, '
                         r'\( \boldsymbol\sigma_{n+1} = \mathbb{C}^{alg} : \boldsymbol\varepsilon^e_{n+1} + '
                         r'\boldsymbol\sigma^{hist}_n \)',
                         ModelArg('elas', 'type:elas_frozen', '', r'Frozen elasticity (`elas_*` generics of '
                                  r'muscle_elasticity_laws).'), None,
                         extra=(_COEFS, _STATE_N)),
               ModelFunc('update', r'State at \( t_{n+1} \) from the converged state and the elastic strain '
                         r'at \( t_{n+1} \)',
                         None, None,
                         extra=(_COEFS, _STATE_N,
                                ModelArg('strain_e', '3D2Osym', '', r'Elastic strain at \( t_{n+1} \).')),
                         outs=(ModelArg('state_np1', 'type:visco_state', '', r'State at \( t_{n+1} \): '
                                        r'`state_n` with the branches of the model updated; the inactive '
                                        r'branches are copied from `state_n` (the historical OO update left '
                                        r'those of `state_np1` untouched: the same values, zero, through the '
                                        r'API).'),)),
               ModelFunc('n_branches', r'Number of history tensors of the state used by the model',
                         ModelArg('res', 'integer', '', r'Number of branches.'), None)),
        fd_scale=None, fd_lower=None, models=VISCOELASTIC_MODELS, oo_funcs=(),
        support=('muscle_visco_types', 'muscle_visco_adapters')),
    # Kinematics (D-017.14): the configuration of a material point from its deformation gradient F.
    # The state is explicit (`kin_state_<k>`, fixed size, muscle_kin_finite: F, F^-1, the polar
    # rotation R, J, C, b, b^-1 and the strain measure of the kinematics) and is written by
    # `kin_update`, the only procedure that validates (`status`, D-017.14a), or, in the corotational
    # kinematics, by `kin_corot_update` from a strain and a rotation (D-017.14d). Only what depends on
    # the kinematics is a family function: the update (which strain measure), the conversion of the
    # constitutive stress (work conjugate of that strain) to Cauchy, and of the constitutive tangent
    # to the spatial one. The strain measures, the stress conversions (PK1, PK2, Cauchy, Kirchhoff),
    # the push-forward/pull-back of tangents and the corotational rotations are the same for every
    # kinematics: functions of the state in muscle_kin_finite. No material parameters (empty
    # parameter types, D-017.14). Only 3D (as the historical kinematics). The OO adapters keep the
    # historical types and bindings, written with the macros of muscle_kin_adapters.
    'kin': ModelFamily(
        prefix='kin', module='muscle_kin_{}', family_module='muscle_kin_laws', par='prm',
        args=(),
        funcs=(ModelFunc('update', r'Kinematic state of a deformation gradient \( \mathbf{F} \) (the '
                         r'deformation gradient is validated)',
                         None, None,
                         extra=(ModelArg('F', '3D2O', '', r'Deformation gradient \( \mathbf{F} \).'),),
                         outs=(ModelArg('state', 'type:kin_state', '', r'Kinematic state (not valid when '
                                        r'`status` is not `STATUS_OK`).'),
                               ModelArg('status', 'integer', '', r'`STATUS_OK` or `STATUS_INVALID_DEFORMATION` '
                                        r'(muscle_status).')),
                         models=KINEMATICS_F_MODELS),
               ModelFunc('corot_update', r'Kinematic state of a corotational kinematics from the strain in the '
                         r'unrotated frame and the rotation (identity if absent)',
                         None, None,
                         extra=(ModelArg('strain', '3D2Osym', '', r'Strain in the unrotated frame.'),),
                         outs=(ModelArg('state', 'type:kin_state', '', r'Kinematic state: strain, '
                                        r'\( J = \exp(\operatorname{tr}\boldsymbol\varepsilon) \) and R; the '
                                        r'F-based components keep the reference configuration.'),
                               ModelArg('status', 'integer', '', r'`STATUS_OK` (no validation).')),
                         models=('corotational',),
                         opts=(ModelArg('R', '3D2O', '', r'Rotation \( \mathbf{R} \) (identity if absent).'),)),
               ModelFunc('to_cauchy', r'Cauchy stress from the constitutive stress (work conjugate of the '
                         r'strain of the kinematics)',
                         ModelArg('res', '3D2Osym', '', r'Cauchy stress \( \boldsymbol\sigma \).'), None,
                         extra=(_KIN_STATE, ModelArg('stress', '3D2Osym', '', r'Constitutive stress.'))),
               ModelFunc('to_spatial_tangent', r'Spatial tangent from the constitutive tangent (no geometric '
                         r'term, D-017.14b)',
                         ModelArg('res', '3D4O2sym', '', r'Spatial tangent.'), None,
                         extra=(_KIN_STATE, ModelArg('tangent', '3D4O2sym', '', r'Constitutive tangent.'))),
               ModelFunc('is_finite_strain', r'Whether the kinematics is a finite strain one',
                         ModelArg('res', 'logical', '', r'`.true.` for a finite strain kinematics.'), None)),
        fd_scale=None, fd_lower=None, models=KINEMATICS_MODELS, oo_funcs=(),
        support=('muscle_kin_finite', 'muscle_kin_adapters')),
}


def model_family(prefix):
    """The ``ModelFamily`` of a prefix: ``model_family('hard').par`` -> ``'prm'``."""
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
    """Module of a model: ``MODEL_MODULES`` if listed there, else ``ModelFamily.module``:
    ``model_module('elas', 'linear')`` -> ``muscle_elasticity_linear``,
    ``model_module('elas', 'frozen')`` -> ``muscle_visco_frozen``."""
    return MODEL_MODULES.get((prefix, model), MODEL_FAMILIES[prefix].module.format(model))


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


def model_proc_doc(prefix, name):
    """FORD summary of a kernel procedure of ``model_procs``."""
    fam = MODEL_FAMILIES[prefix]
    if fam.joint and name == fam.joint.name:
        return fam.joint.doc
    if fam.joint and name == fam.joint.flag:
        return fam.joint.flag_doc
    return model_func(prefix, name).doc


def model_proc_specifics(prefix, name, model, kinds):
    """Specific names of the kernel procedure ``name`` (of ``model_procs``) of a model, for every
    kind of ``kinds`` and every dimension of the family (the joint flag has no dimension):
    the module procedures of the family generic ``<prefix>_<name>`` for that model."""
    fam = MODEL_FAMILIES[prefix]
    dims = [None] if fam.joint and name == fam.joint.flag else model_dims(prefix)
    return [model_specific(prefix, name, model, k, d) for d in dims for k in kinds]


def model_specifics(prefix, model, kinds):
    """Specific names of every kernel procedure of a model (``model_procs``) in ``kinds``: what the
    model module makes public and the family module imports."""
    return [s for n in model_procs(prefix, model) for s in model_proc_specifics(prefix, n, model, kinds)]


def model_joint(prefix, model, joint=None):
    """Whether a model writes its own joint evaluation: ``ModelJoint.models`` for a registered
    model; a model not registered in the family (a test-only model) may set it with ``joint``
    (``None``: not joint). For a registered model ``joint`` must agree with the registry."""
    fam = MODEL_FAMILIES[prefix]
    reg = fam.joint is not None and model in fam.joint.models
    if joint is None:
        return reg
    if model in fam.models and bool(joint) != reg:
        raise ValueError(f'{prefix} {model}: joint={joint} disagrees with ModelJoint.models')
    if joint and fam.joint is None:
        raise ValueError(f'family {prefix!r} has no joint evaluation')
    return bool(joint)


def model_joint_outs(prefix):
    """``ModelArg`` of the ``intent(out)`` arguments of the joint evaluation of a family: the
    result of each function of ``ModelJoint.funcs`` named as in ``ModelJoint.outs``."""
    j = MODEL_FAMILIES[prefix].joint
    return [model_func(prefix, f).result._replace(name=o) for f, o in zip(j.funcs, j.outs)]


def model_decl(typ, shape, kind, dim=None):
    """Fortran type spec of an argument or parameter (``ModelArg``/``ModelParam`` ``type``):
    ``'real'`` -> ``real(dp)``, ``'integer'``, ``'logical'``, a tensor family -> ``type(ten_3D2Osym_dp)``
    (``dim`` required), or a tensor family with its dimension, ``'2D4O3sym'`` ->
    ``type(ten_2D4O3sym_dp)`` (a parameter of a fixed dimension, e.g. a tangent cached per dimension
    in a parameter type shared by both), or ``'type:<name>'``, a derived type with one instance per
    kind, ``'type:visco_state'`` -> ``type(visco_state_dp)``. ``shape`` is not part of the type spec
    (it follows the name)."""
    if typ.startswith('type:') and typ[5:].isidentifier():
        return f'type({typ[5:]}_{kind})'
    if typ == 'real':
        return f'real({kind})'
    if typ in ('integer', 'logical'):
        return typ
    if _dim_family(typ):
        return f'type({tname(typ, kind)})'
    if typ in FAMILIES:
        if dim is None:
            raise ValueError(f'tensor type {typ!r} needs a dimension')
        return f'type({tname(dim + typ, kind)})'
    raise ValueError(f'unknown model type {typ!r}')


def _dim_family(typ):
    """Whether ``typ`` is a tensor family with its dimension (``'3D4O3sym'``)."""
    return typ[:2] in DIMS and typ[2:] in FAMILIES


def model_default(param, kind):
    """Default initialization of a parameter component (``' = 0.0_dp'``) or ``''``: a
    real/integer/logical parameter, scalar or array (every element gets the scalar default), or a
    tensor of a fixed dimension (``'3D2Osym'``) of rank-1 storage, whose default is a value for
    every stored component (``' = ten_3D2Osym_dp([0.0_dp, ...])'``)."""
    if param.default is None:
        return ''
    if _dim_family(param.type) and not param.shape:
        lay = layout(param.type[2:], param.type[:2])
        if len(lay.shape) != 1:
            raise ValueError(f'parameter {param.name!r}: default of a tensor stored as a matrix')
        lit = f'{param.default}_{kind}'
        return ' = {}([{}])'.format(tname(param.type, kind), ', '.join([lit] * lay.shape[0]))
    if param.type not in ('real', 'integer', 'logical'):
        raise ValueError(f'parameter {param.name!r}: defaults only for intrinsic types and fixed-dimension tensors')
    return ' = ' + (f'{param.default}_{kind}' if param.type == 'real' else str(param.default))


def model_params_init_doc(prefix, model, params, validated=False):
    """FORD clause of the parameter type of a model with private parameters (``''`` if none is
    private): built only with the initializer ``<prefix>_<model>_init``, which computes the private
    derived components (D-017.10b), or, with ``validated``, which validates the parameters (every
    component private, none derived: e.g. the Prony series of Maxwell, D-017.13):
    ``model_params_init_doc('visco', 'kelvin', [('k', 'real', '', 'K.', '0.0', True)], True)`` ->
    `` (build it with `visco_kelvin_init`, which validates the private parameters)``."""
    ps = model_params(params)
    if not any(p.private for p in ps):
        return ''
    init = f'{prefix}_{model}_init'
    if validated:
        return f' (build it with `{init}`, which validates the private parameters)'
    return (f' (build it with `{init}`, which computes the private derived components; do not change the '
            f'inputs after `{init}`)')


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


def model_oo_extra(params, extra):
    """OO-only components of a historical OO adapter (``ModelParam`` from tuples as in
    ``model_params``): state that the historical type keeps but the kernel parameter type does not,
    e.g. the tangents cached by ``Elasticity_linear`` (D-017.12e). They follow the parameters in the
    OO type and are not copied to the kernel type; their names may not repeat a parameter name:
    ``[p.name for p in model_oo_extra([('mu', 'real', '', 'M.', None)],
    [('flag', 'logical', '', 'F.', '.false.')])]`` -> ``['flag']``."""
    res = model_params(extra or ())
    clash = {p.name.lower() for p in model_params(params)} & {p.name.lower() for p in res}
    if clash:
        raise ValueError(f'OO-only components repeat parameter names: {sorted(clash)}')
    return res


def model_composed(prefix, model):
    """Whether a model of the family composes another law (``ModelFamily.composed``):
    ``model_composed('vp', 'jc')`` -> ``True``, ``model_composed('vp', 'rk')`` -> ``False``."""
    return model in MODEL_FAMILIES[prefix].composed


def model_comp_values(prefix, func):
    """``CompValue`` of the composition values a function receives (``ModelFunc.comp``), in
    order: ``[v.arg.name for v in model_comp_values('vp', 'dflow_dep')]`` -> ``['sy', 'dsy_dep']``."""
    fam = MODEL_FAMILIES[prefix]
    names = model_func(prefix, func).comp
    if not names:
        return []
    values = {v.arg.name: v for v in fam.composition.values}
    return [values[n] for n in names]


def model_oo_args(prefix, func):
    """Arguments of a function in the OO adapter (and of every model that composes nothing):
    the family arguments and the ``extra`` ones of the function."""
    return list(MODEL_FAMILIES[prefix].args) + list(model_func(prefix, func).extra)


def model_args(prefix, func, model):
    """Arguments of a function of a model, after the parameter argument: the family arguments,
    the ``extra`` ones of the function and, if the model composes another law, the composition
    values the function receives (D-017.6b): ``[a.name for a in model_args('vp', 'dstress_dep',
    'jc')]`` -> ``['ep', 'epd', 'dt', 'sy', 'dsy_dep']``."""
    res = model_oo_args(prefix, func)
    if model_composed(prefix, model):
        res += [v.arg for v in model_comp_values(prefix, func)]
    return res


def model_call(prefix, func, model, kind, xarg=None, dim=None):
    """Call of a model function with its arguments (``model_args``), the argument ``xarg``
    replaced by ``{x}`` (the placeholder of the finite-difference macros): ``model_call('hard',
    'stress', 'bilinear', 'dp', 'ep')`` -> ``hard_stress_bilinear_dp(prm, {x})``."""
    fam = MODEL_FAMILIES[prefix]
    args = [fam.par] + ['{x}' if a.name == xarg else a.name for a in model_args(prefix, func, model)]
    return '{}({})'.format(model_specific(prefix, func, model, kind, dim), ', '.join(args))


def model_oo_bindings(prefix, bindings=None, model=None):
    """Bindings of an OO adapter as ``(binding, family function, OO dummy names or None,
    dimension)``: ``bindings`` if given, else ``ModelFamily.oo_bindings`` if the family has
    historical names, else every function of ``ModelFamily.oo_funcs`` (``None``: all) under its
    own name with the family argument names. The dimension (of the tensor arguments and result) is
    the optional fourth element of a given binding, by default ``OO_DIM`` in a family with
    dimensions and ``None`` otherwise. With ``model``, the default bindings only cover the
    functions that model has (``ModelFunc.models``), and a given binding of a function the model
    does not have is an error."""
    fam = MODEL_FAMILIES[prefix]
    own = model_procs(prefix, model)
    if bindings is None:
        if fam.oo_bindings is not None:
            bindings = [b for b in fam.oo_bindings if b[1] in own]
        else:
            bindings = [(f.name, f.name, None) for f in fam.funcs
                        if (fam.oo_funcs is None or f.name in fam.oo_funcs) and f.name in own]
    for b in bindings:
        if b[1] not in own:
            raise ValueError(f'OO binding {b[0]!r}: {prefix}_{b[1]} is not a function of {model!r}')
    dflt = OO_DIM if fam.dims else None
    res = []
    for b in bindings:
        if len(b) not in (3, 4):
            raise ValueError(f'OO binding {b!r}: expected (binding, function, names[, dimension])')
        dim = b[3] if len(b) == 4 else dflt
        if fam.dims and dim not in fam.dims or not fam.dims and dim is not None:
            raise ValueError(f'OO binding {b[0]!r}: dimension {dim!r} not in {fam.dims}')
        res.append((b[0], b[1], b[2], dim))
    return res


def model_oo_decl(prefix, arg, kind, dim):
    """Declaration of an argument of an OO adapter procedure (without attributes): as
    ``model_decl``, but a tensor argument is ``class(...)`` in a family whose historical base
    declares it so (``ModelFamily.oo_class``)."""
    decl = model_decl(arg.type, arg.shape, kind, dim)
    if MODEL_FAMILIES[prefix].oo_class and decl.startswith('type('):
        decl = 'class(' + decl[len('type('):]
    return decl


def model_default_expr(prefix, func, model, kind, dim=None):
    """Right-hand side of the default body of a function (``ModelFunc.default``) with every
    ``{name}`` replaced by the call of that function of the model: ``model_default_expr('vp',
    'dstress_dep', 'rk', 'dp')`` -> ``vp_dflow_dep_rk_dp(prm, ep, epd) +
    vp_dflow_depd_rk_dp(prm, ep, epd)/dt``."""
    f = model_func(prefix, func)
    if f.default is None:
        raise ValueError(f'{prefix}_{func} has no default formula')
    calls = {g.name: model_call(prefix, g.name, model, kind, dim=dim) for g in MODEL_FAMILIES[prefix].funcs
             if g.name != func}
    return f.default.format(**calls)


def isotropic_4o3sym(lay, diag, off, shear, zero):
    """Stored components of an isotropic ``4O3sym`` tensor
    \\( a\\,\\delta_{ij}\\delta_{kl} + b\\,(\\delta_{ik}\\delta_{jl} + \\delta_{il}\\delta_{jk}) \\)
    of the layout ``lay`` (3D or 2D), in storage order, as ``[(position, expr)]`` from its three
    distinct values, Fortran expressions: ``diag`` for \\( C_{iiii} \\) (a + 2b), ``off`` for
    \\( C_{iijj} \\), i != j (a), ``shear`` for \\( C_{ijij} \\), i != j (b), and ``zero`` for every
    other component: ``isotropic_4o3sym(layout('4O3sym', '2D'), 'e1', 'lam', 'mu', '0')[3]`` ->
    ``(4, 'mu')`` (xyxy)."""
    if lay.family != '4O3sym':
        raise ValueError(f'isotropic_4o3sym: layout of family {lay.family!r}')
    out = []
    for c in lay.comps:
        i, j, k, l = c.idx
        if i == j and k == l:
            expr = diag if i == k else off
        elif i != j and (i, j) == (k, l):
            expr = shear
        else:
            expr = zero
        out.append((c.pos[0], expr))
    return out


def literal_value(text):
    """Value of a real literal with or without kind (``'1.0'``, ``'1.0_dp'``, ``'2e-3_sp'``), or
    ``None`` if ``text`` is ``None`` or any other expression (a variable, ``max(a, b)``, ...)."""
    if text is None:
        return None
    try:
        return float(str(text).strip().split('_')[0])
    except ValueError:
        return None
