# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Unit tests of tools/fypp/muscle_meta.py.

The expected storage orders are copied from the original Fortran sources (`init2` bodies and
type documentation in src/tensors/types/). Run: python3 -m unittest tests/tools/test_muscle_meta.py
"""

import itertools
import pathlib
import random
import re
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / 'tools' / 'fypp'))
import muscle_meta as mm  # noqa: E402

# Inventory of the historical dp operator signatures (tools/fypp/ops_signatures.py)
INVENTORY = pathlib.Path(__file__).resolve().parents[2] / 'tools' / 'fypp' / 'ref' / 'ops_signatures_dp.txt'


def names(family, dim):
    return [c.name for c in mm.layout(family, dim).comps]


class TestSecondOrder(unittest.TestCase):

    def test_storage_orders(self):
        # init2_ten_3D2O / init2_ten_2D2O: "Internal storage is column-major"
        self.assertEqual(names('2O', '3D'), ['xx', 'yx', 'zx', 'xy', 'yy', 'zy', 'xz', 'yz', 'zz'])
        self.assertEqual(names('2O', '2D'), ['xx', 'yx', 'xy', 'yy', 'zz'])
        # inverse_ten_3D2Osym / inverse_ten_2D2Osym: (xx, yy, zz, xy, yz, xz) / (xx, yy, zz, xy)
        self.assertEqual(names('2Osym', '3D'), ['xx', 'yy', 'zz', 'xy', 'yz', 'xz'])
        self.assertEqual(names('2Osym', '2D'), ['xx', 'yy', 'zz', 'xy'])

    def test_positions_are_sequential(self):
        for family, dim in itertools.product(('2O', '2Osym'), mm.DIMS):
            lay = mm.layout(family, dim)
            self.assertEqual([c.pos for c in lay.comps], [(s,) for s in range(1, lay.shape[0] + 1)])


class TestFourthOrder(unittest.TestCase):

    def test_3d4o3sym_order(self):
        # init2_ten_3D4O3sym
        self.assertEqual(names('4O3sym', '3D'), [
            'xxxx', 'yyyy', 'zzzz', 'xyxy', 'yzyz', 'xzxz',
            'xxyy', 'yyzz', 'zzxy', 'xyyz', 'yzxz',
            'xxzz', 'yyxy', 'zzyz', 'xyxz',
            'xxxy', 'yyyz', 'zzxz',
            'xxyz', 'yyxz',
            'xxxz'])

    def test_2d4o3sym_order(self):
        # init2_ten_2D4O3sym; Voigt pairs (11, 22, 33, 44, 12, 23, 34, 13, 24, 14)
        self.assertEqual(names('4O3sym', '2D'), [
            'xxxx', 'yyyy', 'zzzz', 'xyxy', 'xxyy', 'yyzz', 'zzxy', 'xxzz', 'yyxy', 'xxxy'])

    def test_4o2sym_matrix(self):
        # init2_ten_3D4O2sym: column J holds (xx, yy, zz, xy, yz, xz) // name(J)
        lay = mm.layout('4O2sym', '3D')
        self.assertEqual(lay.shape, (6, 6))
        by_pos = {c.pos: c.name for c in lay.comps}
        self.assertEqual(by_pos[(1, 4)], 'xxxy')
        self.assertEqual(by_pos[(4, 1)], 'xyxx')
        self.assertEqual(by_pos[(5, 6)], 'yzxz')
        self.assertEqual(mm.layout('4O2sym', '2D').shape, (4, 4))

    def test_4o_matrix_uses_2o_index(self):
        # ten_3D4O: I(i,j) of ten_3D2O, e.g. map6 = (1, 5, 9, 4, 8, 7) for xx, yy, zz, xy, yz, xz;
        # ten_2D4O: I(1,1)=1, I(2,1)=2, I(1,2)=3, I(2,2)=4, I(3,3)=5
        by_idx = {c.idx: c.pos for c in mm.layout('4O', '3D').comps}
        self.assertEqual(by_idx[(1, 2, 2, 3)], (4, 8))
        self.assertEqual(by_idx[(1, 3, 3, 3)], (7, 9))
        by_idx = {c.idx: c.pos for c in mm.layout('4O', '2D').comps}
        self.assertEqual(by_idx[(2, 1, 3, 3)], (2, 5))
        self.assertEqual(mm.layout('4O', '2D').shape, (5, 5))


class TestConsistency(unittest.TestCase):

    def test_names_match_indices_and_are_unique(self):
        for family, dim in itertools.product(mm.FAMILIES, mm.DIMS):
            lay = mm.layout(family, dim)
            self.assertEqual(len({c.name for c in lay.comps}), len(lay.comps), (family, dim))
            self.assertEqual(len({c.pos for c in lay.comps}), len(lay.comps), (family, dim))
            for c in lay.comps:
                self.assertEqual(c.name, ''.join('xyz'[i - 1] for i in c.idx))

    def test_component_counts(self):
        expected = {('2O', '3D'): 9, ('2O', '2D'): 5, ('2Osym', '3D'): 6, ('2Osym', '2D'): 4,
                    ('4O', '3D'): 81, ('4O', '2D'): 25, ('4O2sym', '3D'): 36,
                    ('4O2sym', '2D'): 16, ('4O3sym', '3D'): 21, ('4O3sym', '2D'): 10}
        for key, n in expected.items():
            self.assertEqual(len(mm.layout(*key).comps), n, key)

    def test_2d_has_no_out_of_plane_shear(self):
        for family in mm.FAMILIES:
            for c in mm.layout(family, '2D').comps:
                for a, b in zip(c.idx[::2], c.idx[1::2]):
                    self.assertFalse(3 in (a, b) and (a, b) != (3, 3), (family, c))

    def test_errors(self):
        with self.assertRaises(ValueError):
            mm.layout('6O', '3D')
        with self.assertRaises(ValueError):
            mm.tname('3D2Osym', 'qp')
        self.assertEqual(mm.tname('3D2Osym', 'sp'), 'ten_3D2Osym_sp')


class TestQueries(unittest.TestCase):

    def test_index_of_stored_components(self):
        for family, dim in itertools.product(mm.FAMILIES, mm.DIMS):
            lay = mm.layout(family, dim)
            for c in lay.comps:
                self.assertEqual(lay.index(*c.idx), c.pos, (family, dim, c))

    def test_index_applies_symmetry(self):
        self.assertEqual(mm.layout('2Osym', '3D').index(2, 1), (4,))
        self.assertEqual(mm.layout('2O', '3D').index(2, 1), (2,))      # unsymmetric: yx
        self.assertEqual(mm.layout('4O2sym', '3D').index(2, 1, 3, 1), (4, 6))
        self.assertEqual(mm.layout('4O2sym', '3D').index(3, 1, 2, 1), (6, 4))  # no major symmetry
        lay = mm.layout('4O3sym', '3D')
        self.assertEqual(lay.index(3, 1, 2, 1), lay.index(1, 2, 1, 3))  # major + minor
        self.assertEqual(lay.index(1, 2, 1, 3), (15,))                  # xyxz

    def test_index_of_2d_structural_zeros(self):
        self.assertIsNone(mm.layout('2Osym', '2D').index(1, 3))
        self.assertIsNone(mm.layout('2O', '2D').index(3, 2))
        self.assertIsNone(mm.layout('4O', '2D').index(1, 1, 2, 3))
        self.assertIsNone(mm.layout('4O3sym', '2D').index(2, 3, 1, 1))
        self.assertEqual(mm.layout('4O3sym', '2D').index(3, 3, 2, 1), (7,))  # zzxy

    def test_multiplicities_cover_the_full_tensor(self):
        # 3D: every component of the full tensor; 2D: the in-plane ones plus 33
        full = {('2', '3D'): 9, ('2', '2D'): 5, ('4', '3D'): 81, ('4', '2D'): 25}
        for family, dim in itertools.product(mm.FAMILIES, mm.DIMS):
            total = sum(c.mult for c in mm.layout(family, dim).comps)
            self.assertEqual(total, full[(family[0], dim)], (family, dim))
        self.assertEqual([c.mult for c in mm.layout('2Osym', '2D').comps], [1, 1, 1, 2])
        mults = {c.name: c.mult for c in mm.layout('4O3sym', '3D').comps}
        self.assertEqual((mults['xxxx'], mults['xxyy'], mults['xyxy'], mults['xxxy'], mults['xyxz']),
                         (1, 2, 4, 4, 8))

    def test_voigt_and_base(self):
        for dim in mm.DIMS:
            self.assertIsNone(mm.layout('2Osym', dim).base)
            for family, base in (('4O', '2O'), ('4O2sym', '2Osym'), ('4O3sym', '2Osym')):
                lay = mm.layout(family, dim)
                self.assertEqual(lay.base, mm.layout(base, dim))
                names = [c.name for c in lay.base.comps]
                for c in lay.comps:
                    self.assertEqual(c.name, names[c.voigt[0] - 1] + names[c.voigt[1] - 1])
                    if len(c.pos) == 2:
                        self.assertEqual(c.voigt, c.pos)


class TestHelpers(unittest.TestCase):

    def test_ref(self):
        lay = mm.layout('2Osym', '2D')
        self.assertEqual(lay.ref('a%vals', 2, 1), 'a%vals(4)')
        self.assertEqual(lay.ref('a%vals', 1, 3, zero='0.0_dp'), '0.0_dp')
        self.assertEqual(mm.layout('4O2sym', '3D').ref('c', 1, 2, 3, 3), 'c(4, 3)')

    def test_matrix(self):
        m = mm.layout('2Osym', '2D').matrix('v', zero='0')
        self.assertEqual(m, [['v(1)', 'v(4)', '0'], ['v(4)', 'v(2)', '0'], ['0', '0', 'v(3)']])
        m = mm.layout('4O3sym', '3D').matrix('v')
        self.assertEqual((len(m), len(m[0])), (6, 6))
        self.assertEqual([m[i][i] for i in range(6)], ['v({})'.format(i) for i in range(1, 7)])
        self.assertEqual((m[0][1], m[1][0], m[0][5]), ('v(7)', 'v(7)', 'v(21)'))
        self.assertEqual(len(mm.layout('4O3sym', '2D').matrix('v')), 4)
        # 2D 4O3sym: the whole 4x4 Voigt matrix is stored (symmetric, no structural zeros)
        m = mm.layout('4O3sym', '2D').matrix('v')
        self.assertEqual(m, [['v(1)', 'v(5)', 'v(8)', 'v(10)'], ['v(5)', 'v(2)', 'v(6)', 'v(9)'],
                             ['v(8)', 'v(6)', 'v(3)', 'v(7)'], ['v(10)', 'v(9)', 'v(7)', 'v(4)']])

    def test_matrix_2o(self):
        # ten_3D2O: column-major 3x3; ten_2D2O: (xx, yx, xy, yy, zz), no out-of-plane coupling
        self.assertEqual(mm.layout('2O', '3D').matrix('v'),
                         [['v(1)', 'v(4)', 'v(7)'], ['v(2)', 'v(5)', 'v(8)'], ['v(3)', 'v(6)', 'v(9)']])
        self.assertEqual(mm.layout('2O', '2D').matrix('v', zero='0'),
                         [['v(1)', 'v(3)', '0'], ['v(2)', 'v(4)', '0'], ['0', '0', 'v(5)']])

    def test_matrix_4o_4o2sym(self):
        # Matrix storage: the matrix shown to the user is vals itself, C(I, J) = vals(I, J), with
        # I, J the storage indices of the base layout (2O for 4O, 2Osym for 4O2sym).
        for family, dim, n in [('4O', '3D', 9), ('4O', '2D', 5), ('4O2sym', '3D', 6), ('4O2sym', '2D', 4)]:
            m = mm.layout(family, dim).matrix('v', zero='0')
            self.assertEqual(m, [['v({}, {})'.format(i, j) for j in range(1, n + 1)]
                                 for i in range(1, n + 1)], (family, dim))
        # Spot checks against the original sources (see test_4o_matrix_uses_2o_index and
        # test_4o2sym_matrix): C_1223 of ten_3D4O at (4, 8), C_1123 of ten_3D4O2sym at (1, 5)
        self.assertEqual(mm.layout('4O', '3D').matrix('c')[3][7], 'c(4, 8)')
        self.assertEqual(mm.layout('4O', '3D').ref('c', 1, 2, 2, 3), 'c(4, 8)')
        self.assertEqual(mm.layout('4O2sym', '3D').matrix('c')[0][4], 'c(1, 5)')
        self.assertEqual(mm.layout('4O2sym', '3D').ref('c', 1, 1, 3, 2), 'c(1, 5)')

    def test_iname(self):
        self.assertEqual(mm.iname('2O'), 'iden_2O')
        self.assertEqual(mm.iname('4O4TS', 'sp'), 'iden_4O4TS_sp')
        with self.assertRaises(ValueError):
            mm.iname('2OS', 'qp')

    def test_iden_matrix(self):
        # 2O: delta_ij
        self.assertEqual(mm.iden_matrix('2O', 'c', 'h', '0'),
                         [['c', '0', '0'], ['0', 'c', '0'], ['0', '0', 'c']])
        # 4O3T: delta_ij delta_kl, ones in the normal block (assign_3D4O3sym_I4O3TS)
        m3 = mm.iden_matrix('4O3T', 'c', 'h', '0')
        self.assertEqual(m3, [['c' if i < 3 and j < 3 else '0' for j in range(6)] for i in range(6)])
        # 4O4T: symmetric identity, diag(1, 1, 1, 1/2, 1/2, 1/2) (assign_3D4O2sym_I4O4T)
        m4 = mm.iden_matrix('4O4T', 'c', 'h', '0')
        diag = ['c', 'c', 'c', 'h', 'h', 'h']
        self.assertEqual(m4, [[diag[i] if i == j else '0' for j in range(6)] for i in range(6)])
        with self.assertRaises(ValueError):
            mm.iden_matrix('4O1T', 'c', 'h', '0')

    def test_join_terms(self):
        self.assertEqual(mm.join_terms(['a', 'b', 'c'], indent=4), 'a &\n    + b &\n    + c')
        self.assertEqual(mm.join_terms(['x']), 'x')

    def test_family_module(self):
        self.assertEqual(mm.family_module('2Osym'), 'muscle_tensor_2osym')
        self.assertEqual(mm.family_module('4O3sym'), 'muscle_tensor_4o3sym')
        self.assertEqual(mm.family_module('iden'), 'muscle_tensor_iden')
        with self.assertRaises(ValueError):
            mm.family_module('3D2Osym')

    def test_norm_terms(self):
        self.assertEqual(mm.norm_terms(mm.layout('2Osym', '2D'), 'a%vals'),
                         ['abs(a%vals(1))', 'abs(a%vals(2))', 'abs(a%vals(3))', '2*abs(a%vals(4))'])
        terms = mm.norm_terms(mm.layout('2Osym', '3D'), 'v')
        self.assertEqual(terms[3:], ['2*abs(v(4))', '2*abs(v(5))', '2*abs(v(6))'])
        self.assertEqual(mm.norm_terms(mm.layout('2O', '2D'), 'v'), ['abs(v({}))'.format(i) for i in range(1, 6)])
        with self.assertRaises(ValueError):
            mm.norm_terms(mm.layout('4O', '3D'), 'v')


def all_entries(kinds=mm.KINDS):
    return [e for g in mm.OP_GENERICS for d in mm.DIMS for k in kinds for e in mm.op_entries(g, d, k)]


def signature(e):
    """Line of tools/fypp/ops_signatures.py for a dp entry (historical names, lower case)."""
    def norm(t):
        t = t.lower()
        if t.startswith('real('):
            return t
        return 'type({})'.format(t[:-3] if t.endswith('_dp') else t)
    args = ', '.join(norm(t) for t in (e.a, e.b) if t is not None)
    sig = '{} ({})'.format(mm.op_interface(e.generic), args)
    return sig if e.res is None else sig + ' -> ' + norm(e.res)


def kindless(e):
    return not any(t.endswith(('_sp', '_dp')) for t in (e.a, e.b) if t)


class TestOperatorTable(unittest.TestCase):

    def test_dp_is_the_historical_inventory(self):
        ref = [l.strip() for l in INVENTORY.read_text(encoding='utf-8').splitlines() if l.strip()]
        sigs = [signature(e) for e in all_entries(('dp',))]
        self.assertEqual(len(sigs), len(set(sigs)))
        self.assertEqual(sorted(sigs), sorted(ref))

    def test_sp_mirrors_dp(self):
        dp = all_entries(('dp',))
        sp = all_entries(('sp',))
        # Only dp: the rows without a kinded operand (iden_4O3T +- iden_4O4T and mirrors)
        n_kindless = sum(1 for g in ('+', '-') for row in mm.OP_TABLE[g]
                         if mm.op_iden(row.a) == (row.a[1:], False) and mm.op_iden(row.b) == (row.b[1:], False)
                         for _ in range(1 + row.mirror))
        self.assertEqual(sum(map(kindless, dp)), n_kindless)
        self.assertEqual(len(sp), len(dp) - n_kindless)
        for e in sp:
            for t in (e.name, e.a, e.b, e.res):
                self.assertNotIn('dp', t or '', e)
            self.assertTrue(e.name.endswith('_sp'), e)
            self.assertFalse(kindless(e), e)
            if e.res is not None:
                self.assertTrue(e.res.endswith('_sp') or e.res == 'real(sp)', e)
        to_sp = {tuple((t or '').replace('_dp', '_sp').replace('(dp)', '(sp)')
                       for t in (e.generic, e.name, e.a, e.b, e.res)) for e in dp if not kindless(e)}
        self.assertEqual(to_sp, {tuple(t or '' for t in (e.generic, e.name, e.a, e.b, e.res)) for e in sp})

    def test_unscaled_identities_have_no_kind(self):
        sp = mm.op_entries('+', '3D', 'sp')
        self.assertIn(mm.OpEntry('+', 'sum_I2O_3D2Osym_sp', 'iden_2O', 'ten_3D2Osym_sp', 'ten_3D2Osym_sp',
                                 False, 'I2O', '3D2Osym', '3D2Osym', '3D', 'sp'), sp)
        self.assertIn(mm.OpEntry('+', 'sum_I4O3TS_I4O4T_sp', 'iden_4O3TS_sp', 'iden_4O4T', 'ten_3D4O3sym_sp',
                                 False, 'I4O3TS', 'I4O4T', '3D4O3sym', '3D', 'sp'), sp)
        self.assertNotIn('sum_I4O3T_I4O4T_sp', [e.name for e in sp])
        self.assertEqual([e for e in mm.op_entries('+', '2D', 'dp') if e.name.startswith('sum_I4O3T_I4O4T')], [])

    def test_names_unique_and_unambiguous(self):
        # Specific names are unique; a generic never has two procedures with the same arguments.
        # Scope: the table only (the generics of the type modules are checked by the compiler).
        for g in mm.OP_GENERICS:
            entries = [e for d in mm.DIMS for k in mm.KINDS for e in mm.op_entries(g, d, k)]
            names = [e.name for e in entries]
            self.assertEqual(len(names), len(set(names)), g)
            args = [(e.a, e.b) for e in entries]
            self.assertEqual(len(args), len(set(args)), g)
        # ... and unique across generics too (one module-level namespace per operator module)
        self.assertEqual(len({e.name for e in all_entries()}), len(all_entries()))

    def test_historical_names(self):
        names = {e.name for e in all_entries(('dp',))}
        for n in ('sum_I2O_3D2Osym', 'sub_2D4O3sym_I4O4T', 'assign_3D4O3sym_I4O3TS', 'ddot_I4O4TS_2D2Osym',
                  'dot_2D2O_2D2Osym', 'mul_I2OS_3D2O', 'tdotsym_2D2Osym', 'tdot_I2O_3D2Osym',
                  'transform_3D2O_3D4O2sym', 'sum_I4O3T_I4O4T'):
            self.assertIn(n + '_dp', names)
        self.assertEqual({g for gs in mm.OP_MODULES.values() for g in gs}, set(mm.OP_TABLE))
        self.assertEqual(mm.op_interface('='), 'assignment(=)')
        self.assertEqual(mm.op_interface('.ddot.'), 'operator(.ddot.)')
        with self.assertRaises(ValueError):
            mm.op_entries('.cross.', '3D', 'dp')
        with self.assertRaises(ValueError):
            mm.op_entries('+', '1D', 'dp')
        with self.assertRaises(ValueError):
            mm.op_entries('+', '3D', 'qp')

    def test_mirror_and_tags(self):
        e = mm.op_entries('.ddot.', '2D', 'dp')
        direct = next(x for x in e if x.name == 'ddot_2D4O2sym_2D2Osym_dp')
        mirror = next(x for x in e if x.name == 'ddot_2D2Osym_2D4O2sym_dp')
        self.assertFalse(direct.swapped)
        self.assertTrue(mirror.swapped)
        self.assertEqual((mirror.a, mirror.b, mirror.res), ('ten_2D2Osym_dp', 'ten_2D4O2sym_dp', 'ten_2D2Osym_dp'))
        self.assertEqual((mirror.tag_a, mirror.tag_b, mirror.tag_res), ('2D2Osym', '2D4O2sym', '2D2Osym'))
        ddot_i = next(x for x in e if x.name == 'ddot_I2OS_2D2O_dp')
        self.assertEqual((ddot_i.tag_res, ddot_i.res), ('real', 'real(dp)'))
        assign = mm.op_entries('=', '3D', 'sp')[0]
        self.assertEqual((assign.b, assign.res, assign.tag_res), ('ten_3D2Osym_sp', None, None))
        self.assertEqual(mm.op_iden('I4O4TS'), ('4O4T', True))
        self.assertEqual(mm.op_iden('I2O'), ('2O', False))
        self.assertIsNone(mm.op_iden('3D2Osym'))
        # A bad row fails at validation, even if never expanded
        saved = mm.OP_TABLE['*']
        try:
            mm.OP_TABLE['*'] = saved + [mm.OpPair('2O', '2O', '2O', True)]
            with self.assertRaises(ValueError):
                mm._check_op_table()
            mm.OP_TABLE['*'] = saved + [mm.OpPair('2O', '5O', '2O')]
            with self.assertRaises(ValueError):
                mm._check_op_table()
        finally:
            mm.OP_TABLE['*'] = saved



# ---------------------------------------------------------------------------------------------
# Phase 4, batch A: helpers of the operator templates
# ---------------------------------------------------------------------------------------------

def feval(expr, env):
    """Evaluates a Fortran expression of the helpers in Python: kinded literals lose their kind,
    array references ``x%vals(i)`` / ``x%vals(i, j)`` become ``env['x'][(i,)]`` / ``[(i, j)]``."""
    py = re.sub(r'(\d+\.\d+)_(sp|dp)\b', r'\1', expr)
    py = re.sub(r'\b(\w+)%vals\(([\d, ]+)\)', lambda m: "V['{}'][({},)]".format(m.group(1), m.group(2)), py)
    py = py.replace('%', '_')
    return eval(py, {'__builtins__': {}}, dict(env))


def rnd_matrix(dim, rng, sym=False):
    """Random 3x3 matrix (0-based lists) with the 2D structural zeros (13, 23, 31, 32) if 2D."""
    m = [[rng.uniform(-2.0, 2.0) for _ in range(3)] for _ in range(3)]
    if sym:
        m = [[m[max(i, j)][min(i, j)] for j in range(3)] for i in range(3)]
    if dim == '2D':
        for a, b in ((0, 2), (1, 2), (2, 0), (2, 1)):
            m[a][b] = 0.0
    return m


def stored(lay, m):
    """Values of the stored components of a 2nd-order layout taken from a 3x3 matrix."""
    return {c.pos: m[c.idx[0] - 1][c.idx[1] - 1] for c in lay.comps}


class TestIdenCoef(unittest.TestCase):

    def test_values(self):
        self.assertEqual(mm.iden_coef('2O', (1, 1)), 1)
        self.assertEqual(mm.iden_coef('2O', (1, 2)), 0)
        self.assertEqual(mm.iden_coef('4O3T', (1, 1, 2, 2)), 1)
        self.assertEqual(mm.iden_coef('4O3T', (1, 2, 1, 2)), 0)
        self.assertEqual(mm.iden_coef('4O4T', (1, 1, 1, 1)), 1)
        self.assertEqual(mm.iden_coef('4O4T', (1, 1, 2, 2)), 0)
        self.assertEqual(mm.iden_coef('4O4T', (1, 2, 1, 2)), 0.5)
        self.assertEqual(mm.iden_coef('4O4T', (1, 2, 2, 1)), 0.5)
        self.assertEqual(mm.iden_coef('4O4T', (2, 3, 3, 2)), 0.5)

    def test_full_tensors(self):
        # delta_ij; delta_ij delta_kl; (delta_ik delta_jl + delta_il delta_jk) / 2
        for idx in itertools.product(range(1, 4), repeat=2):
            self.assertEqual(mm.iden_coef('2O', idx), 1 if idx[0] == idx[1] else 0, idx)
        for i, j, k, l in itertools.product(range(1, 4), repeat=4):
            self.assertEqual(mm.iden_coef('4O3T', (i, j, k, l)), 1 if i == j and k == l else 0)
            if i == j == k == l:
                expected = 1
            elif i != j and {(k, l), (l, k)} == {(i, j), (j, i)}:
                expected = 0.5
            else:
                expected = 0
            self.assertEqual(mm.iden_coef('4O4T', (i, j, k, l)), expected, (i, j, k, l))
        values = {mm.iden_coef(n, idx) for n in ('4O3T', '4O4T')
                  for idx in itertools.product(range(1, 4), repeat=4)}
        self.assertEqual(values, {0, 1, 0.5})

    def test_errors(self):
        for name, idx in (('4O1T', (1, 1, 1, 1)), ('2OS', (1, 1)), ('4O4TS', (1, 1, 1, 1)),
                          ('2O', (1, 1, 1, 1)), ('2O', (1,)), ('4O3T', (1, 1)), ('4O4T', (1, 1, 1))):
            with self.assertRaises(ValueError, msg=(name, idx)):
                mm.iden_coef(name, idx)

    def test_iden_matrix_matches_iden_coef(self):
        sym = {1: 'one', 0.5: 'half', 0: 'zero'}
        self.assertEqual(mm.iden_matrix('2O', 'one', 'half', 'zero'),
                         [[sym[mm.iden_coef('2O', (i, j))] for j in (1, 2, 3)] for i in (1, 2, 3)])
        voigt = [(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (1, 3)]
        for name in ('4O3T', '4O4T'):
            self.assertEqual(mm.iden_matrix(name, 'one', 'half', 'zero'),
                             [[sym[mm.iden_coef(name, r + c)] for c in voigt] for r in voigt], name)


def comb_env(a, b):
    return {'a': a, 'b': b}


class TestIdenCombination(unittest.TestCase):

    def test_second_order_single(self):
        I = [('2O', None, '+')]
        self.assertEqual(mm.iden_combination(mm.layout('2Osym', '3D'), I, 'dp'),
                         [('1', '1.0_dp'), ('2', '1.0_dp'), ('3', '1.0_dp')])
        self.assertEqual(mm.iden_combination(mm.layout('2Osym', '2D'), I, 'sp'),
                         [('1', '1.0_sp'), ('2', '1.0_sp'), ('3', '1.0_sp')])
        self.assertEqual(mm.iden_combination(mm.layout('2O', '3D'), I, 'dp'),
                         [('1', '1.0_dp'), ('5', '1.0_dp'), ('9', '1.0_dp')])
        self.assertEqual(mm.iden_combination(mm.layout('2O', '2D'), I, 'dp'),
                         [('1', '1.0_dp'), ('4', '1.0_dp'), ('5', '1.0_dp')])
        # scaled, both signs
        self.assertEqual(mm.iden_combination(mm.layout('2O', '2D'), [('2O', 'b%val', '+')], 'dp'),
                         [('1', 'b%val'), ('4', 'b%val'), ('5', 'b%val')])
        self.assertEqual(mm.iden_combination(mm.layout('2Osym', '3D'), [('2O', 'b%val', '-')], 'dp'),
                         [('1', '-b%val'), ('2', '-b%val'), ('3', '-b%val')])
        self.assertEqual(mm.iden_combination(mm.layout('2Osym', '2D'), [('2O', None, '-')], 'dp'),
                         [('1', '-1.0_dp'), ('2', '-1.0_dp'), ('3', '-1.0_dp')])

    def test_fourth_order_single(self):
        lay = mm.layout('4O3sym', '3D')
        self.assertEqual(mm.iden_combination(lay, [('4O4T', 'b%val', '+')], 'dp'),
                         [('1', 'b%val'), ('2', 'b%val'), ('3', 'b%val'),
                          ('4', '0.5_dp * b%val'), ('5', '0.5_dp * b%val'), ('6', '0.5_dp * b%val')])
        self.assertEqual(mm.iden_combination(lay, [('4O3T', None, '+')], 'dp'),
                         [(str(p), '1.0_dp') for p in (1, 2, 3, 7, 8, 12)])  # xxxx..zzzz, xxyy, yyzz, xxzz
        self.assertEqual(mm.iden_combination(mm.layout('4O3sym', '2D'), [('4O4T', None, '+')], 'sp'),
                         [('1', '1.0_sp'), ('2', '1.0_sp'), ('3', '1.0_sp'), ('4', '0.5_sp')])
        self.assertEqual(mm.iden_combination(mm.layout('4O3sym', '2D'), [('4O4T', 'b%val', '-')], 'sp'),
                         [('1', '-b%val'), ('2', '-b%val'), ('3', '-b%val'), ('4', '-0.5_sp * b%val')])
        # matrix storage: column-major subscripts 'I, J'
        self.assertEqual(mm.iden_combination(mm.layout('4O2sym', '3D'), [('4O3T', None, '+')], 'dp'),
                         [('{}, {}'.format(i, j), '1.0_dp') for j in (1, 2, 3) for i in (1, 2, 3)])
        self.assertEqual(mm.iden_combination(mm.layout('4O2sym', '2D'), [('4O4T', 'b%val', '+')], 'dp'),
                         [('1, 1', 'b%val'), ('2, 2', 'b%val'), ('3, 3', 'b%val'),
                          ('4, 4', '0.5_dp * b%val')])

    def test_two_identities(self):
        lay = mm.layout('4O3sym', '3D')
        # docstring example: unscaled + scaled keep their order and signs
        out = dict(mm.iden_combination(lay, [('4O3T', None, '+'), ('4O4T', 'b%val', '-')], 'dp'))
        self.assertEqual(out['1'], '1.0_dp - b%val')
        self.assertEqual(out['4'], '-0.5_dp * b%val')
        self.assertEqual(out['7'], '1.0_dp')
        self.assertEqual(sorted(out, key=int), ['1', '2', '3', '4', '5', '6', '7', '8', '12'])
        # both scaled
        out = dict(mm.iden_combination(mm.layout('4O2sym', '2D'),
                                       [('4O3T', 'a%val', '+'), ('4O4T', 'b%val', '+')], 'dp'))
        self.assertEqual(out['1, 1'], 'a%val + b%val')
        self.assertEqual(out['2, 1'], 'a%val')
        self.assertEqual(out['4, 4'], '0.5_dp * b%val')
        self.assertNotIn('4, 1', out)
        out = dict(mm.iden_combination(lay, [('4O4T', 'a%val', '-'), ('4O3T', 'b%val', '-')], 'sp'))
        self.assertEqual(out['1'], '-a%val - b%val')
        self.assertEqual(out['5'], '-0.5_sp * a%val')
        self.assertEqual(out['8'], '-b%val')

    def test_unscaled_terms_fold(self):
        lay = mm.layout('4O3sym', '3D')
        out = mm.iden_combination(lay, [('4O3T', None, '+'), ('4O4T', None, '+')], 'dp')
        self.assertEqual(out, [('1', '2.0_dp'), ('2', '2.0_dp'), ('3', '2.0_dp'), ('4', '0.5_dp'),
                               ('5', '0.5_dp'), ('6', '0.5_dp'), ('7', '1.0_dp'), ('8', '1.0_dp'),
                               ('12', '1.0_dp')])
        # 1 - 1 = 0 on the diagonal: the component disappears
        out = mm.iden_combination(lay, [('4O3T', None, '+'), ('4O4T', None, '-')], 'dp')
        self.assertEqual(out, [('4', '-0.5_dp'), ('5', '-0.5_dp'), ('6', '-0.5_dp'), ('7', '1.0_dp'),
                               ('8', '1.0_dp'), ('12', '1.0_dp')])
        # everything cancels
        self.assertEqual(mm.iden_combination(mm.layout('2Osym', '2D'),
                                             [('2O', None, '+'), ('2O', None, '-')], 'dp'), [])
        # a scaled term prevents the folding
        out = dict(mm.iden_combination(lay, [('4O3T', None, '+'), ('4O4T', None, '-'),
                                             ('4O4T', 'b%val', '+')], 'dp'))
        self.assertEqual(out['1'], '1.0_dp - 1.0_dp + b%val')

    def test_numeric_all_layouts(self):
        # The expressions evaluate to sum(sign * var * coef) for every stored component, and the
        # components left out are exactly the zero ones.
        a, b = 1.75, -0.375
        cases = [('2O', [('2O', 'a', '+')]), ('2O', [('2O', None, '-'), ('2O', 'b', '+')]),
                 ('4O', [('4O3T', 'a', '+'), ('4O4T', 'b', '-')]),
                 ('4O', [('4O4T', None, '+'), ('4O3T', 'b', '+')]),
                 ('4O', [('4O4T', None, '-'), ('4O3T', None, '+')])]
        for (order, terms), dim, kind in itertools.product(cases, mm.DIMS, mm.KINDS):
            for family in [f for f in mm.FAMILIES if f[0] == order[0]]:
                lay = mm.layout(family, dim)
                out = dict(mm.iden_combination(lay, terms, kind))
                for c in lay.comps:
                    sub = ', '.join(map(str, c.pos))
                    exp = sum((1 if s == '+' else -1) * (1.0 if v is None else {'a': a, 'b': b}[v])
                              * mm.iden_coef(n, c.idx) for n, v, s in terms)
                    if sub in out:
                        self.assertNotIn('_' + ('dp' if kind == 'sp' else 'sp'), out[sub])
                        self.assertAlmostEqual(feval(out[sub], comb_env(a, b)), exp, places=14,
                                               msg=(family, dim, terms, c))
                    else:
                        self.assertEqual(exp, 0.0, (family, dim, terms, c))

    def test_errors(self):
        lay = mm.layout('2Osym', '3D')
        for sign in ('*', '', '+-', '++', ' +'):
            with self.assertRaises(ValueError, msg=sign):
                mm.iden_combination(lay, [('2O', None, sign)], 'dp')
        with self.assertRaises(ValueError):
            mm.iden_combination(lay, [('2O', None, '+')], 'qp')
        with self.assertRaises(ValueError):
            mm.iden_combination(lay, [('4O1T', None, '+')], 'dp')
        with self.assertRaises(ValueError):  # 4th-order identity on a 2nd-order layout
            mm.iden_combination(lay, [('4O4T', None, '+')], 'dp')


class TestOpLayoutUses(unittest.TestCase):

    def test_op_layout(self):
        for family, dim in itertools.product(mm.FAMILIES, mm.DIMS):
            self.assertIs(mm.op_layout(dim + family), mm.layout(family, dim))
        for tag in ('2Osym', '4O3sym', 'I2O', 'I4O4TS', 'real', None, '3D6O'):
            with self.assertRaises(ValueError, msg=tag):
                mm.op_layout(tag)

    def test_op_uses_matches_entries(self):
        templated = {'muscle_tensor_' + f.lower() for f in mm.FAMILIES} | {'muscle_tensor_iden'}
        for module, kinds in itertools.product(mm.OP_MODULES, (['dp'], ['sp'], ['sp', 'dp'])):
            uses = mm.op_uses(module, kinds)
            expected = {t for g in mm.OP_MODULES[module] for d in mm.DIMS for k in kinds
                        for e in mm.op_entries(g, d, k) for t in (e.a, e.b, e.res)
                        if t is not None and not t.startswith('real(')}
            got = [t for _, ts in uses for t in ts]
            self.assertEqual(len(got), len(set(got)), (module, kinds))
            self.assertEqual(set(got), expected, (module, kinds))
            groups = []
            for mod, ts in uses:
                msg = (module, kinds, mod)
                self.assertIn(mod, templated, msg)   # never a compatibility module
                self.assertTrue(ts, msg)
                suffix = {t.rsplit('_', 1)[1] if t.endswith(('_sp', '_dp')) else None for t in ts}
                self.assertEqual(len(suffix), 1, msg)  # one kind (or kindless) per entry
                kind = suffix.pop()
                self.assertIn(kind, list(kinds) + [None], msg)
                if kind is None:
                    self.assertEqual(mod, 'muscle_tensor_iden', msg)
                    for t in ts:
                        self.assertIsNotNone(mm.op_iden('I' + t[len('iden_'):]), msg)
                        self.assertFalse(mm.op_iden('I' + t[len('iden_'):])[1], msg)
                groups.append((mod, kind))
            self.assertEqual(len(groups), len(set(groups)), (module, kinds))

    def test_op_uses_examples(self):
        uses = dict(((m, tuple(ts)) for m, ts in mm.op_uses('assignment', ['dp'])))
        self.assertEqual(uses['muscle_tensor_2o'], ('ten_3D2O_dp', 'ten_2D2O_dp'))
        self.assertIn(('muscle_tensor_iden', ['iden_2O', 'iden_4O3T', 'iden_4O4T']),
                      mm.op_uses('addition_subtraction', ['sp']))
        for module in mm.OP_MODULES:
            for mod, _ in mm.op_uses(module, ['sp', 'dp']):
                self.assertNotRegex(mod, r'muscle_tensor_(\dd|iden_)')
        with self.assertRaises(ValueError):
            mm.op_uses('contraction_triple', ['dp'])


class TestOpCongruence(unittest.TestCase):

    def check(self, dim, seed):
        rng = random.Random(seed)
        F = rnd_matrix(dim, rng)
        S = rnd_matrix(dim, rng, sym=True)
        R = [[sum(F[i][k] * S[k][l] * F[j][l] for k in range(3) for l in range(3)) for j in range(3)]
             for i in range(3)]
        lf, ls = mm.layout('2O', dim), mm.layout('2Osym', dim)
        env = {'V': {'f': stored(lf, F), 's': stored(ls, S)}}
        tmp, res = mm.op_congruence(dim, 'f%vals', 's%vals')
        defined = set()
        for (i, L), expr in tmp:
            self.assertIsNone(re.search(r'\bT\d\d\b', expr), expr)  # only F and S
            env['T{}{}'.format(i, L)] = feval(expr, env)
            defined.add('T{}{}'.format(i, L))
        used = set()
        self.assertEqual(set(res), {c.idx for c in ls.comps})
        for (i, j), expr in res.items():
            used |= set(re.findall(r'\bT\d\d\b', expr))
            self.assertAlmostEqual(feval(expr, env), R[i - 1][j - 1], delta=1e-12 * max(1.0, abs(R[i - 1][j - 1])),
                                   msg=(dim, (i, j), expr))
        self.assertEqual(used, defined)  # every temporary is used and every used one is defined
        return tmp, res

    def test_3d(self):
        for seed in range(5):
            tmp, res = self.check('3D', seed)
        self.assertEqual(len(tmp), 9)
        self.assertEqual(res[(1, 2)], 'T11*f%vals(2) + T12*f%vals(5) + T13*f%vals(8)')

    def test_2d(self):
        for seed in range(5):
            tmp, res = self.check('2D', 100 + seed)
        # out-of-plane 33: the single product F33 F33 S33, no temporary
        self.assertEqual(res[(3, 3)], 'f%vals(5)*f%vals(5)*s%vals(3)')
        self.assertTrue(all(i != 3 and L != 3 for (i, L), _ in tmp), tmp)
        for expr in [e for _, e in tmp] + list(res.values()):
            self.assertNotIn('0', re.sub(r'\(\d\)', '', expr).replace('T', ''), expr)

    def test_temp_name(self):
        tmp, res = mm.op_congruence('3D', 'F%vals', 'S%vals', t='tmp')
        self.assertTrue(all('tmp' in e for e in res.values()))


def full_4o_minor(dim, rng):
    """Random 4th-order tensor with minor symmetries (no major), {(i,j,k,l): value}, 1-based;
    in 2D the components with an out-of-plane shear pair are zero."""
    def ok(p):
        return dim == '3D' or 3 not in p or p == (3, 3)
    C = {}
    for i, j, k, l in itertools.product(range(1, 4), repeat=4):
        key = (min(i, j), max(i, j), min(k, l), max(k, l))
        if key not in C:
            C[key] = rng.uniform(-2.0, 2.0) if ok((i, j)) and ok((k, l)) else 0.0
    return {(i, j, k, l): C[(min(i, j), max(i, j), min(k, l), max(k, l))]
            for i, j, k, l in itertools.product(range(1, 4), repeat=4)}


class TestOpVoigtTransform(unittest.TestCase):

    def check(self, dim, kind, seed):
        rng = random.Random(seed)
        A = rnd_matrix(dim, rng)
        C = full_4o_minor(dim, rng)
        base = mm.layout('2Osym', dim).comps
        n = len(base)
        T = mm.op_voigt_transform(dim, 'A%vals', kind)
        env = {'V': {'A': stored(mm.layout('2O', dim), A)}}
        for expr in T.values():
            self.assertNotIn('_' + ('dp' if kind == 'sp' else 'sp'), expr)
        Tm = [[feval(T[(I, J)], env) if (I, J) in T else 0.0 for J in range(1, n + 1)]
              for I in range(1, n + 1)]
        Cv = [[C[r.idx + c.idx] for c in base] for r in base]
        TC = [[sum(Tm[I][K] * Cv[K][J] for K in range(n)) for J in range(n)] for I in range(n)]
        R = [[sum(TC[I][K] * Tm[J][K] for K in range(n)) for J in range(n)] for I in range(n)]
        for I, r in enumerate(base):
            for J, c in enumerate(base):
                i, j = r.idx
                k, l = c.idx
                ref = sum(A[i - 1][a - 1] * A[j - 1][b - 1] * A[k - 1][cc - 1] * A[l - 1][d - 1]
                          * C[(a, b, cc, d)] for a, b, cc, d in itertools.product(range(1, 4), repeat=4))
                self.assertAlmostEqual(R[I][J], ref, delta=1e-12 * max(1.0, abs(ref)),
                                       msg=(dim, kind, r.idx, c.idx))
        return T, Tm

    def test_3d(self):
        for kind, seed in itertools.product(mm.KINDS, range(3)):
            T, _ = self.check('3D', kind, seed)
            self.assertEqual(len(T), 36)
        T = mm.op_voigt_transform('3D', 'A%vals', 'dp')
        self.assertEqual(T[(1, 4)], '2.0_dp*A%vals(1)*A%vals(4)')
        self.assertEqual(mm.op_voigt_transform('3D', 'A%vals', 'sp')[(1, 4)], '2.0_sp*A%vals(1)*A%vals(4)')

    def test_2d_no_structural_zeros(self):
        for kind, seed in itertools.product(mm.KINDS, range(3)):
            T, Tm = self.check('2D', kind, 50 + seed)
            for (I, J), expr in T.items():
                self.assertNotEqual(Tm[I - 1][J - 1], 0.0, (I, J, expr))  # random A: not zero
                for idx in re.findall(r'A%vals\((\d+)\)', expr):
                    self.assertLessEqual(int(idx), 5, expr)
        # the entries left out are those that vanish for any 2D tensor A (33 row/column coupling)
        T = mm.op_voigt_transform('2D', 'A%vals', 'dp')
        self.assertEqual(sorted(T), sorted((I, J) for I in range(1, 5) for J in range(1, 5)
                                           if (I == 3) == (J == 3)))
        self.assertEqual(T[(3, 3)], 'A%vals(5)*A%vals(5)')


class TestJoinTerms(unittest.TestCase):

    def test_without_lead_is_unchanged(self):
        self.assertEqual(mm.join_terms(['a', 'b']), 'a &\n            + b')
        self.assertEqual(mm.join_terms(['a', 'b'], op='*', indent=2), 'a &\n  * b')
        self.assertEqual(mm.join_terms(['x']), 'x')

    def test_lead(self):
        self.assertEqual(mm.join_terms(['a', 'b'], lead=8), 'a + b')
        self.assertEqual(mm.join_terms(['a', 'b', 'c'], op='*', lead=0), 'a * b * c')
        self.assertEqual(mm.join_terms(['x'], lead=100), 'x')
        # exactly 132 columns fits; 133 does not
        terms = ['a' * 10] * 5   # 50 + 4 * 3 = 62 columns
        self.assertEqual(mm.join_terms(terms, lead=70), ' + '.join(terms))
        self.assertEqual(mm.join_terms(terms, indent=4, lead=71), ' &\n    + '.join(terms))
        long = ['W(:,{})*T(j,{})'.format(k, k) for k in range(1, 10)]
        out = mm.join_terms(long, indent=28, lead=28)
        self.assertEqual(out, ' &\n{}+ '.format(' ' * 28).join(long))
        self.assertTrue(all(len(line) <= 132 for line in out.splitlines()))


    def test_assign_terms(self):
        self.assertEqual(mm.assign_terms('W11', ['a', 'b'], 8), 'W11 = a + b')
        # the operators of the continuation lines go under the first term
        long = ['W(:,{0})*T(j,{0})'.format(k) for k in range(1, 10)]
        out = mm.assign_terms('res%vals(:,j)', long, 12)
        self.assertEqual(out, 'res%vals(:,j) = ' + ' &\n{}+ '.format(' ' * 28).join(long))
        self.assertTrue(all(12 * (n == 0) + len(line) <= 132 for n, line in enumerate(out.splitlines())))
        # exactly 132 columns fits: col + len('x = ') + len(term) = 132
        self.assertEqual(mm.assign_terms('x', ['a' * 120], 8), 'x = ' + 'a' * 120)
        self.assertIn('&', mm.assign_terms('x', ['a' * 60, 'b' * 60], 8))


def entry(generic, name):
    """The OpEntry ``name`` of ``generic`` (any dimension and kind)."""
    return [e for e in all_entries() if e.generic == generic and e.name == name][0]


class TestOpDocsAndArgs(unittest.TestCase):

    def test_dim_kind(self):
        for g in mm.OP_GENERICS:
            for d in mm.DIMS:
                for k in mm.KINDS:
                    for e in mm.op_entries(g, d, k):
                        self.assertEqual((e.dim, e.kind), (d, k))
                        self.assertTrue(e.name.endswith('_' + k), e)

    def test_op_args(self):
        self.assertEqual(mm.op_args(entry('+', 'sum_I2O_3D2Osym_dp')), ['I2', 'a'])
        self.assertEqual(mm.op_args(entry('+', 'sum_3D4O3sym_I4O4TS_sp')), ['a', 'I4'])
        self.assertEqual(mm.op_args(entry('+', 'sum_I4O3TS_I4O4T_dp')), ['I3', 'I4'])
        self.assertEqual(mm.op_args(entry('*', 'dot_2D2O_2D2Osym_sp')), ['a', 'b'])
        self.assertEqual(mm.op_args(entry('.tdotsym.', 'tdotsym_3D2Osym_dp')), ['a'])
        self.assertEqual(mm.op_args(entry('=', 'assign_2D4O2sym_I4O4TS_dp')), ['a', 'I4'])
        # one distinct name per operand everywhere
        for e in all_entries():
            args = mm.op_args(e)
            self.assertEqual(len(args), 1 if e.b is None else 2, e)
            self.assertEqual(len({a.lower() for a in args}), len(args), e)

    def test_op_symbols(self):
        self.assertEqual(mm.op_symbols(entry('+', 'sum_I2OS_3D2O_dp')),
                         [r'c\mathbf{I}', r'\mathbf{A}', r'\mathbf{res}'])
        self.assertEqual(mm.op_symbols(entry('-', 'sub_2D4O2sym_I4O3T_sp')),
                         [r'\mathbb{A}', r'\mathbb{I}_{3T}', r'\mathbb{res}'])
        self.assertEqual(mm.op_symbols(entry('+', 'sum_I4O3TS_I4O4TS_dp')),
                         [r'c_1\mathbb{I}_{3T}', r'c_2\mathbb{I}^S', r'\mathbb{res}'])
        self.assertEqual(mm.op_symbols(entry('+', 'sum_3D2O_3D2Osym_dp'), 'ab'),
                         [r'\mathbf{a}', r'\mathbf{b}', r'\mathbf{res}'])
        self.assertEqual(mm.op_symbols(entry('=', 'assign_3D2O_3D2Osym_dp')),
                         [r'\mathbf{A}', r'\mathbf{B}', None])
        self.assertEqual(mm.op_symbols(entry('.ddot.', 'ddot_I2O_3D2Osym_dp'))[-1], 'res')
        self.assertEqual(mm.op_symbols(entry('.tdotsym.', 'tdotsym_2D2Osym_sp'), 'ab'),
                         [r'\mathbf{a}', r'\mathbb{res}'])

    def test_op_doc(self):
        e = entry('+', 'sum_I2O_3D2Osym_dp')
        self.assertEqual(mm.op_doc(e, 'Sum'), 'Sum (3D, `real(dp)`).')
        # runs of spaces (fypp line continuations) count as one: always a single line
        self.assertEqual(mm.op_doc(entry('+', 'sum_I2O_2D2Osym_sp'), 'Sum \\(a  +\n   b\\)'),
                         r'Sum \(a + b\) (2D, `real(sp)`).')

    def test_op_iden_term(self):
        self.assertEqual(mm.op_iden_term('I2OS', 'I2%val'), ('2O', 'I2%val', '+'))
        self.assertEqual(mm.op_iden_term('I4O4T', 'I4%val', '-'), ('4O4T', None, '-'))
        with self.assertRaises(ValueError):
            mm.op_iden_term('3D2Osym', 'a%val')
        # the term feeds iden_combination
        lay = mm.layout('2Osym', '2D')
        self.assertEqual(mm.iden_combination(lay, [mm.op_iden_term('I2OS', 'I2%val')], 'dp'),
                         [('1', 'I2%val'), ('2', 'I2%val'), ('3', 'I2%val')])

    def test_op_mirror(self):
        for e in all_entries():
            if e.swapped:
                m = mm.op_mirror(e)
                self.assertFalse(m.swapped)
                self.assertEqual((m.generic, m.dim, m.kind, m.res), (e.generic, e.dim, e.kind, e.res))
                self.assertEqual((m.a, m.b, m.tag_a, m.tag_b), (e.b, e.a, e.tag_b, e.tag_a))
        self.assertEqual(mm.op_mirror(entry('.tdotsym.', 'tdotsym_I2O_2D2Osym_sp')).name,
                         'tdotsym_2D2Osym_I2O_sp')
        with self.assertRaises(ValueError):
            mm.op_mirror(entry('.tdotsym.', 'tdotsym_2D2Osym_I2O_sp'))


class TestOpDotTerms(unittest.TestCase):

    def test_examples(self):
        l2, s2 = mm.layout('2O', '2D'), mm.layout('2Osym', '2D')
        self.assertEqual(mm.op_dot_terms(l2, 'a%vals', s2, 'b%vals', 3, 3), ['a%vals(5)*b%vals(3)'])
        self.assertEqual(mm.op_dot_terms(l2, 'a%vals', s2, 'b%vals', 1, 1),
                         ['a%vals(1)*b%vals(1)', 'a%vals(3)*b%vals(4)'])
        self.assertEqual(mm.op_dot_terms(l2, 'a%vals', s2, 'b%vals', 1, 3), [])

    def test_numeric(self):
        rng = random.Random(7)
        for d in mm.DIMS:
            for fa, fb in itertools.product(('2O', '2Osym'), repeat=2):
                la, lb = mm.layout(fa, d), mm.layout(fb, d)
                ma, mb = rnd_matrix(d, rng, fa == '2Osym'), rnd_matrix(d, rng, fb == '2Osym')
                env = {'V': {'a': stored(la, ma), 'b': stored(lb, mb)}}
                for i, j in itertools.product((1, 2, 3), repeat=2):
                    terms = mm.op_dot_terms(la, 'a%vals', lb, 'b%vals', i, j)
                    got = feval(' + '.join(terms), env) if terms else 0.0
                    want = sum(ma[i - 1][k] * mb[k][j - 1] for k in range(3))
                    self.assertAlmostEqual(got, want, places=12, msg=(d, fa, fb, i, j))


def voigt_matrix(lay, rng):
    """Random symmetric Voigt matrix (0-based lists) of a 4O3sym layout and its stored values."""
    n = len(lay.base.comps)
    C = [[0.0] * n for _ in range(n)]
    for c in lay.comps:
        C[c.voigt[0] - 1][c.voigt[1] - 1] = C[c.voigt[1] - 1][c.voigt[0] - 1] = rng.uniform(-1, 1)
    return C, {c.pos: C[c.voigt[0] - 1][c.voigt[1] - 1] for c in lay.comps}


class TestOpDdotTerms(unittest.TestCase):

    def test_sym_sum(self):
        l3, l2 = mm.layout('2O', '3D'), mm.layout('2O', '2D')
        self.assertEqual(mm.op_sym_sum(l3, 'a%vals', 1, 2), '(a%vals(4) + a%vals(2))')
        self.assertEqual(mm.op_sym_sum(l2, 'a%vals', 1, 2), '(a%vals(3) + a%vals(2))')
        self.assertEqual(mm.op_sym_sum(l3, 'a%vals', 3, 3), 'a%vals(9)')

    def test_weight_form(self):
        # One form everywhere (D-015.8): the shear weight on the right factor
        self.assertEqual(mm.op_weight('b%vals(4)', 2, 'sp'), '2.0_sp*b%vals(4)')
        self.assertEqual(mm.op_weight('b%vals(1)', 1, 'dp'), 'b%vals(1)')
        self.assertEqual(mm.op_weighted('x', 'y', 1, 'dp'), 'x*y')
        self.assertEqual(mm.op_weighted('x', 'y', 2, 'dp'), 'x*(2.0_dp*y)')

    def test_voigt_terms_example(self):
        l4 = mm.layout('4O3sym', '2D')
        self.assertEqual(mm.op_ddot_voigt_terms(l4, 'a%vals', 'b%vals', 1, 'dp'),
                         ['a%vals(1)*b%vals(1)', 'a%vals(5)*b%vals(2)', 'a%vals(8)*b%vals(3)',
                          'a%vals(10)*(2.0_dp*b%vals(4))'])
        self.assertEqual(mm.op_ddot_voigt_terms(mm.layout('4O3sym', '3D'), 'a%vals', 'b%vals', 4,
                                                'dp')[3], 'a%vals(4)*(2.0_dp*b%vals(4))')

    def test_voigt_terms_numeric(self):
        # (A : b)_I = sum_J A_IJ w_J b_J against the full tensor contraction A_ijkl b_kl
        rng = random.Random(5)
        for d in mm.DIMS:
            lay = mm.layout('4O3sym', d)
            base = lay.base.comps
            C, vals = voigt_matrix(lay, rng)
            mb = rnd_matrix(d, rng, sym=True)
            env = {'V': {'a': vals, 'b': stored(lay.base, mb)}}
            for r in base:
                got = feval(' + '.join(mm.op_ddot_voigt_terms(lay, 'a%vals', 'b%vals', r.pos[0],
                                                               'dp')), env)
                want = sum(C[r.pos[0] - 1][lay.base.index(k, l)[0] - 1] * mb[k - 1][l - 1]
                           for k in (1, 2, 3) for l in (1, 2, 3) if lay.base.index(k, l))
                self.assertAlmostEqual(got, want, places=12, msg=(d, r))

    def test_matrix_terms(self):
        l2 = mm.layout('4O3sym', '2D')
        self.assertEqual(mm.op_ddot_matrix_terms(l2, 'a%vals', l2, 'b%vals', 1, 1, 'dp'),
                         ['a%vals(1)*b%vals(1)', 'a%vals(5)*b%vals(5)', 'a%vals(8)*b%vals(8)',
                          'a%vals(10)*(2.0_dp*b%vals(10))'])

    def test_matrix_terms_numeric(self):
        # Term structure of C_IJ = sum_K A_IK w_K B_KJ (K ascending, both dimensions). Same formula as
        # the helper: the independent check is test_muscle_tensor_ops_contraction_double_full (full tensors).
        rng = random.Random(9)
        for d in mm.DIMS:
            lay = mm.layout('4O3sym', d)
            n = len(lay.base.comps)
            A, va = voigt_matrix(lay, rng)
            B, vb = voigt_matrix(lay, rng)
            w = [c.mult for c in lay.base.comps]
            env = {'V': {'a': va, 'b': vb}}
            for I, J in itertools.product(range(1, n + 1), repeat=2):
                terms = mm.op_ddot_matrix_terms(lay, 'a%vals', lay, 'b%vals', I, J, 'dp')
                self.assertEqual(len(terms), n)
                want = sum(A[I - 1][K] * w[K] * B[K][J - 1] for K in range(n))
                self.assertAlmostEqual(feval(' + '.join(terms), env), want, places=12, msg=(d, I, J))


class TestOpVoigtCongruence(unittest.TestCase):

    def test_2d_example(self):
        V = mm.op_voigt_congruence(mm.layout('4O3sym', '2D'), 'A%vals', 'C_in%vals', 'dp')
        self.assertEqual(V.n, 4)
        self.assertEqual(V.T, mm.op_voigt_transform('2D', 'A%vals', 'dp'))
        self.assertEqual(V.rows, {1: [1, 2, 4], 2: [1, 2, 4], 3: [3], 4: [1, 2, 4]})
        self.assertEqual(V.res[0], ('1', ['W11*T11', 'W12*T12', 'W14*T14']))
        self.assertEqual(V.w[0], ((1, 1), ['T11*C_in%vals(1)', 'T12*C_in%vals(5)', 'T14*C_in%vals(10)']))

    def test_structure(self):
        for d in mm.DIMS:
            for f in ('4O2sym', '4O3sym'):
                lay = mm.layout(f, d)
                V = mm.op_voigt_congruence(lay, 'A%vals', 'C_in%vals', 'dp')
                self.assertEqual(V.n, len(lay.base.comps))
                for I, ks in V.rows.items():
                    self.assertEqual(ks, sorted(K for (r, K) in V.T if r == I))
                # one entry of R per stored component, in storage order
                self.assertEqual([p for p, _ in V.res], [', '.join(map(str, c.pos)) for c in lay.comps])
                # W holds exactly the entries the components of R read, sorted
                used = {tuple(map(int, t.split('*')[0][1:])) for _, terms in V.res for t in terms}
                self.assertEqual([IK for IK, _ in V.w], sorted(used))

    def test_numeric(self):
        # R_IJ = sum_K W_IK T_JK, W_IK = sum_L T_IL C_LK, with the entries of T as numbers
        rng = random.Random(11)
        for d in mm.DIMS:
            lay = mm.layout('4O3sym', d)
            n = len(lay.base.comps)
            V = mm.op_voigt_congruence(lay, 'A%vals', 'C_in%vals', 'dp')
            C = [[0.0] * n for _ in range(n)]
            for c in lay.comps:
                C[c.voigt[0] - 1][c.voigt[1] - 1] = C[c.voigt[1] - 1][c.voigt[0] - 1] = rng.uniform(-1, 1)
            Tn = {IJ: rng.uniform(-1, 1) for IJ in V.T}
            env = {'V': {'C_in': {c.pos: C[c.voigt[0] - 1][c.voigt[1] - 1] for c in lay.comps}}}
            env.update({'T{}{}'.format(*IJ): v for IJ, v in Tn.items()})
            for (I, K), terms in V.w:
                env['W{}{}'.format(I, K)] = feval(' + '.join(terms), env)
            for (pos, terms), c in zip(V.res, lay.comps):
                I, J = c.voigt
                want = sum(Tn.get((I, L), 0.0) * C[L - 1][K - 1] * Tn.get((J, K), 0.0)
                           for L in range(1, n + 1) for K in range(1, n + 1))
                self.assertAlmostEqual(feval(' + '.join(terms), env), want, places=12, msg=(d, pos))



class TestModelFamilies(unittest.TestCase):

    def test_model_label(self):
        self.assertEqual(mm.model_label('hard', 'vocemod'), 'modified Voce')
        self.assertEqual(mm.model_label('hard', 'ludwik'), 'Ludwik')
        self.assertEqual(mm.model_label('hard', 'fdpow'), 'fdpow')   # not listed: the identifier
        for (prefix, model) in mm.MODEL_LABELS:
            self.assertIn(model, mm.MODEL_FAMILIES[prefix].models)
    """Model contract metadata (D-017.6) and the registries of the model families."""

    ROOT = pathlib.Path(__file__).resolve().parents[2]

    #: Templates of a family that follow ``ModelFamily.module`` but are not models (shared types).
    SUPPORT = {'visco': {'muscle_visco_types'}, 'kin': {'muscle_kin_finite'}, 'yield': {'muscle_yield_anisotropic'}}

    def test_registry_matches_family_files(self):
        # Every model template of a family is registered, and every registered model has one
        # (the module of a model is muscle_<prefix>_<model>, D-027; support templates of a family are not
        # models)
        fams = mm.MODEL_FAMILIES.values()
        family_modules = {f.family_module for f in fams}
        for fam in fams:
            pattern = fam.module.format('*') + '.fypp'
            others = {mm.model_module(o.prefix, m) for o in fams if o is not fam for m in o.models}
            files = {p.stem for p in (self.ROOT / 'src').rglob(pattern)} - family_modules - others
            own = {fam.module.format(m) for m in fam.models}
            self.assertEqual(files, own | self.SUPPORT.get(fam.prefix, set()), fam.prefix)
            for m in fam.models:
                self.assertTrue(list((self.ROOT / 'src').rglob(mm.model_module(fam.prefix, m) + '.fypp')),
                                (fam.prefix, m))
            self.assertTrue(set(fam.composed) <= set(fam.models), fam.prefix)
            self.assertEqual(not fam.composition, not fam.composed, fam.prefix)
        self.assertEqual(mm.model_family('hard').models, mm.HARDENING_MODELS)
        self.assertEqual(mm.model_family('vp').models, mm.VISCOPLASTIC_MODELS)
        self.assertEqual(mm.model_family('os').models, mm.OVERSTRESS_MODELS)

    def test_names(self):
        self.assertEqual(mm.model_type('hard', 'bilinear', 'sp'), 'hard_bilinear_sp')
        self.assertEqual(mm.model_specific('hard', 'dstress_dep', 'bilinear', 'dp'),
                         'hard_dstress_dep_bilinear_dp')

    def test_composition(self):
        # Composed models receive the composition values each function needs (ModelFunc.comp)
        self.assertTrue(mm.model_composed('vp', 'jc'))
        self.assertFalse(mm.model_composed('vp', 'rk'))
        self.assertFalse(mm.model_composed('hard', 'bilinear'))
        self.assertEqual([v.func for v in mm.model_comp_values('vp', 'dflow_dep')], ['stress', 'dstress_dep'])
        self.assertEqual(mm.model_comp_values('hard', 'stress'), [])

    def test_families_are_consistent(self):
        for fam in mm.MODEL_FAMILIES.values():
            funcs = [f.name for f in fam.funcs]
            self.assertEqual(len(set(funcs)), len(funcs), fam.prefix)
            for f in fam.funcs:
                self.assertIs(mm.model_func(fam.prefix, f.name), f)
                self.assertTrue(set(f.comp) <= {v.name for v in fam.composition}, (fam.prefix, f.name))
            # the composed values are functions of the composed family (hardening laws)
            self.assertTrue({v.func for v in fam.composition} <= set(mm.model_procs('hard')), fam.prefix)
        with self.assertRaises(KeyError):
            mm.model_func('hard', 'nope')

    def test_model_params_validation(self):
        p = mm.model_params([('a', 'real', '', 'A.', '0.5'), ('n', 'integer', '', 'N.', 2),
                             ('v', 'real', '(3)', 'V.', None)])
        self.assertEqual([x.name for x in p], ['a', 'n', 'v'])
        with self.assertRaises(ValueError):
            mm.model_params([('a', 'real', '', '', None), ('A', 'real', '', '', None)])
        with self.assertRaises(ValueError):
            mm.model_params([('a', 'complex', '', '', None)])

    def test_yield_family(self):
        # Tensor argument: the specific names carry the dimension, except the joint flag
        self.assertEqual(mm.model_family('yield').models, mm.YIELD_MODELS)
        self.assertEqual(mm.model_specific('yield', 'stress_eq', 'hill48', 'dp', '3D'),
                         'yield_stress_eq_hill48_3D_dp')
        self.assertEqual(mm.model_specific('yield', 'joint_evaluation', 'hill48', 'sp'),
                         'yield_joint_evaluation_hill48_sp')
        with self.assertRaises(ValueError):
            mm.model_specific('yield', 'stress_eq', 'hill48', 'dp')
        self.assertEqual(mm.model_dims('yield'), ['3D'])
        self.assertEqual(mm.model_dims('hard'), [None])
        self.assertEqual(mm.model_procs('yield'), ['stress_eq', 'dstress_eq', 'ddstress_eq',
                                                   'stress_eq_and_derivatives', 'joint_evaluation'])
        self.assertEqual(mm.model_procs('hard'), ['stress', 'dstress_dep', 'ddstress_ddep'])
        # private parameters (6th element), public by default
        p = mm.model_params([('a', 'real', '', 'A.', None), ('B', 'real', '', 'B.', None, True)])
        self.assertEqual([x.private for x in p], [False, True])

    def test_families_dims_and_joint(self):
        for fam in mm.MODEL_FAMILIES.values():
            funcs = [f.name for f in fam.funcs]
            if fam.dims is not None:
                self.assertTrue(set(fam.dims) <= set(mm.DIMS), fam.prefix)
            if fam.joint is not None:
                j = fam.joint
                self.assertTrue(set(j.funcs) <= set(funcs), fam.prefix)
                self.assertNotIn(j.name, funcs)
                self.assertNotIn(j.flag, funcs)
        self.assertTrue(set(mm.YIELD_JOINT) <= set(mm.YIELD_MODELS))

    def test_kinematics_family(self):
        self.assertEqual(mm.model_family('kin').models, mm.KINEMATICS_MODELS)
        self.assertEqual(mm.model_procs('kin'), ['update', 'update_norot', 'corot_update', 'to_cauchy', 'to_spatial_tangent',
                                                 'geometric_tangent', 'is_finite_strain'])
        # Functions of only some models (ModelFunc.models): the corotational kinematics is updated
        # from a strain and a rotation, the others from F
        self.assertEqual(mm.model_procs('kin', 'corotational'),
                         ['corot_update', 'to_cauchy', 'to_spatial_tangent', 'geometric_tangent', 'is_finite_strain'])
        self.assertEqual(mm.model_procs('kin', 'total_lagrangian'),
                         ['update', 'update_norot', 'to_cauchy', 'to_spatial_tangent', 'geometric_tangent', 'is_finite_strain'])
        self.assertEqual(mm.model_procs('hard', 'bilinear'), mm.model_procs('hard'))
        # Newton derivative of the viscoplastic laws with a rate floor only (D-018.19b)
        self.assertEqual(mm.model_procs('vp', 'jc')[-1], 'newton_dstress_dep')
        self.assertNotIn('newton_dstress_dep', mm.model_procs('vp', 'rk'))
        for fam in mm.MODEL_FAMILIES.values():
            for f in fam.funcs:
                self.assertTrue(f.models is None or set(f.models) <= set(fam.models), (fam.prefix, f.name))
        self.assertEqual(mm.model_specific('kin', 'update', 'total_lagrangian', 'sp'),
                         'kin_update_total_lagrangian_sp')
        self.assertEqual(mm.model_params([]), ())

    def test_elastic_family(self):
        # 3D and 2D
        self.assertEqual(mm.model_dims('elas'), ['3D', '2D'])
        self.assertEqual(mm.model_specific('elas', 'stress', 'linear', 'sp', '2D'), 'elas_stress_linear_2D_sp')
        # tensor parameter of a fixed dimension
        self.assertEqual(mm.model_params([('t', '3D4O3sym', '', 'T.', None, True)])[0].type, '3D4O3sym')
        with self.assertRaises(ValueError):
            mm.model_params([('t', '4D4O3sym', '', 'T.', None)])

    def test_viscoelastic_family(self):
        # First family with a state (D-017.13): 3D only (no dimension in names)
        self.assertEqual(mm.model_procs('visco'), ['step_coefficients', 'freeze', 'update', 'n_branches',
                                                   'freeze_branches', 'update_branches'])
        self.assertEqual(mm.model_specific('visco', 'update', 'maxwell', 'sp'), 'visco_update_maxwell_sp')
        self.assertEqual(mm.model_dims('visco'), [None])
        # The frozen elasticity is an elastic model (its file is in src/viscoelasticity)
        self.assertIn('frozen', mm.model_family('elas').models)
        self.assertEqual(mm.model_module('elas', 'frozen'), 'muscle_elas_frozen')
        self.assertEqual(mm.model_module('visco', 'kelvin'), 'muscle_visco_kelvin')

    def test_target_case_signatures(self):
        # The signatures of the generated family target tests (tests/gpu/target_cases.py) cover the
        # functions of each family they call, and nothing else
        sys.path.insert(0, str(self.ROOT / 'tests' / 'gpu'))
        import target_cases as tc
        for prefix, (args, results) in tc.FAMILY_SIGNATURES.items():
            fam = mm.model_family(prefix)
            self.assertEqual(list(results), [f.name for f in fam.funcs], prefix)
            self.assertTrue(set(results.values()) <= {'real'} | set(mm.FAMILIES), prefix)
            self.assertEqual(fam.dims is None, all(tc.FAMILY_ARG_TYPES.get(a, 'real') == 'real' for a in args), prefix)
        self.assertEqual(set(tc.FAMILY_SIGNATURES), set(tc.FAMILY_PREFIXES) - {'visco'})
        self.assertTrue({p for p, _ in tc.FAMILY_EXTRA_ARGS} <= set(tc.FAMILY_SIGNATURES))

    def test_literal_value(self):
        self.assertEqual(mm.literal_value('1.0'), 1.0)
        self.assertEqual(mm.literal_value('2e-3_sp'), 2e-3)
        self.assertIsNone(mm.literal_value('sc'))
        self.assertIsNone(mm.literal_value('max(a, b)'))
        self.assertIsNone(mm.literal_value(None))

if __name__ == '__main__':
    unittest.main()
