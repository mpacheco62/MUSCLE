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
