# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Unit tests of tools/fypp/muscle_meta.py.

The expected storage orders are copied from the original Fortran sources (`init2` bodies and
type documentation in src/tensors/types/). Run: python3 -m unittest tests/tools/test_muscle_meta.py
"""

import itertools
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / 'tools' / 'fypp'))
import muscle_meta as mm  # noqa: E402


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


if __name__ == '__main__':
    unittest.main()
