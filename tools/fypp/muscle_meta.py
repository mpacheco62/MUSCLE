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

_AXES = 'xyz'


def join_terms(terms, op='+', indent=12):
    """Joins Fortran terms with ``op``, one term per line (free-form continuation), so that long
    expressions stay below the 132-column limit: ``join_terms(['a', 'b'])`` ->
    ``'a &\n            + b'``."""
    return ' &\n{}{} '.format(' ' * indent, op).join(terms)


def tname(base, kind):
    """Name of the concrete tensor type: ``tname('3D2Osym', 'dp')`` -> ``'ten_3D2Osym_dp'``."""
    if kind not in KINDS:
        raise ValueError(f"unknown kind '{kind}', expected one of {KINDS}")
    return f"ten_{base}_{kind}"


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
