# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Unit tests of the declared material registry (tools/fypp/muscle_materials.py) and of its validation
and helpers in tools/fypp/muscle_meta.py (D-008, D-018.9).

Run: python3 -m unittest tests/tools/test_muscle_materials.py
"""

import contextlib
import io
import pathlib
import subprocess
import sys
import tempfile
import textwrap
import unittest

FYPP_DIR = pathlib.Path(__file__).resolve().parents[2] / 'tools' / 'fypp'
sys.path.insert(0, str(FYPP_DIR))
import muscle_meta as mm  # noqa: E402

Material = mm.Material


def cp(name='cp_a', **kw):
    """A valid closest_point material, with fields overridden by ``kw``."""
    fields = dict(elas='linear', yield_='vonmises', hard='swift')
    fields.update(kw)
    return Material(name, 'closest_point', **fields)


class TestDefaultRegistry(unittest.TestCase):

    def test_loads_and_validates(self):
        mats = mm.load_materials()
        self.assertGreater(len(mats), 0)
        self.assertEqual(len({m.name for m in mats}), len(mats))

    def test_covers_every_model_of_each_solver(self):
        # D-008: the test registry covers every model a solver supports at least once
        mats = mm.load_materials()
        for solver, spec in mm.MATERIAL_SOLVERS.items():
            used = [m for m in mats if m.solver == solver]
            for prefix in spec.families:
                supported = set(mm.MODEL_FAMILIES[prefix].models) - set(spec.unsupported.get(prefix, ()))
                covered = {getattr(m, mm.MATERIAL_FIELDS[prefix]) for m in used}
                self.assertEqual(covered, supported, f'{solver}: {prefix} models not covered')

    def test_material_lookup(self):
        m = mm.material('cp_vm_swift')
        self.assertEqual((m.solver, m.elas, m.yield_, m.hard), ('closest_point', 'linear', 'vonmises', 'swift'))
        with self.assertRaises(ValueError):
            mm.material('no_such_material')


class TestValidation(unittest.TestCase):

    def assertInvalid(self, mats, text):
        with self.assertRaises(ValueError) as ctx:
            mm.validate_materials(mats)
        self.assertIn(text, str(ctx.exception))

    def test_valid(self):
        self.assertEqual(mm.validate_materials([cp('a'), cp('b', hard='bilinear')]),
                         (cp('a'), cp('b', hard='bilinear')))

    def test_empty_and_duplicates(self):
        self.assertInvalid([], 'empty')
        self.assertInvalid([cp('a'), cp('a', hard='ludwik')], 'declared twice')

    def test_not_a_material(self):
        self.assertInvalid([('a', 'closest_point')], 'not a muscle_meta.Material')

    def test_names(self):
        for bad in ('CP', '1cp', 'cp-a', '', 'cp a', '_cp'):
            self.assertInvalid([cp(bad)], 'lower case Fortran identifier')
        # Longest identifier: mat_tangent_<name>_sp / _dp (15 extra characters) <= 63
        longest_ok = 'a' * (mm.FORTRAN_NAME_MAX - len('mat_tangent__dp'))
        mm.validate_material(cp(longest_ok))
        self.assertInvalid([cp(longest_ok + 'b')], 'name too long')

    def test_solver(self):
        self.assertInvalid([Material('a', 'newton', 'linear', 'vonmises', 'swift')], 'unknown solver')

    def test_families(self):
        self.assertInvalid([cp(hard=None)], "needs a 'hard' model")
        self.assertInvalid([cp(yield_=None)], "needs a 'yield' model")
        self.assertInvalid([cp(vp='jc')], "does not use a 'vp' model")
        self.assertInvalid([cp(os='perzyna')], "does not use an 'os' model")

    def test_registered_models(self):
        self.assertInvalid([cp(yield_='tresca')], "'tresca' is not a registered 'yield' model")
        self.assertInvalid([cp(hard='vonmises')], "'vonmises' is not a registered 'hard' model")
        self.assertInvalid([cp(elas='maxwell')], "'maxwell' is not a registered 'elas' model")

    def test_composed_families(self):
        # A law that composes a hardening law takes it from the 'hard' component (MaterialSolver.composes)
        mm.validate_material(Material('a', 'closest_point_time', 'linear', 'vonmises', 'swift', vp='jc'))
        mm.validate_material(Material('a', 'closest_point_time', 'linear', 'vonmises', vp='rk'))
        self.assertInvalid([Material('a', 'closest_point_time', 'linear', 'vonmises', vp='jc')],
                           "needs a 'hard' model")
        self.assertInvalid([Material('a', 'closest_point_time', 'linear', 'vonmises', 'swift', vp='rk')],
                           "composes none")
        self.assertInvalid([Material('a', 'closest_point_perzyna', 'linear', 'vonmises', os='perzyna')],
                           "needs a 'hard' model")
        self.assertInvalid([Material('a', 'closest_point_perzyna', 'linear', 'vonmises', 'swift')],
                           "needs an 'os' model")
        self.assertInvalid([Material('a', 'closest_point_time', 'linear', 'vonmises', 'swift', vp='jc',
                                     os='perzyna')], "does not use an 'os' model")

    def test_viscoelastic(self):
        # The frozen viscoelastic elasticity comes from a viscoelastic model (Phase 7 step 6)
        self.assertInvalid([cp(elas='frozen')], "needs a 'visco' model")
        self.assertInvalid([cp(visco='maxwell')], "does not use a 'visco' model")
        self.assertInvalid([cp(elas='frozen', visco='bogus')], "not a registered 'visco' model")
        self.assertInvalid([cp(visco_branches=2)], 'without a viscoelastic elasticity')
        for bad in (-1, mm.MAX_VISCO_BRANCHES + 1, 2.0, True):
            self.assertInvalid([cp(elas='frozen', visco='maxwell', visco_branches=bad)], 'visco_branches must be')
        m = mm.validate_material(cp(elas='frozen', visco='maxwell', visco_branches=mm.MAX_VISCO_BRANCHES))
        self.assertEqual(mm.material_components(m), [('visco', 'maxwell'), ('yield', 'vonmises'), ('hard', 'swift')])
        self.assertEqual(mm.material_state_size(m), 7 + 6 + 6 * mm.MAX_VISCO_BRANCHES)
        self.assertEqual(mm.material_state(m, 'sp'), ('muscle_material_state', 'visco_plastic_state_sp'))
        self.assertEqual(mm.material_state_size(cp()), 7)
        self.assertEqual(mm.material_state(cp(), 'dp'), ('muscle_plastic_state', 'plastic_state_dp'))

    def test_kinematics(self):
        # Finite strain materials (D-018.15): a registered kinematics, the tangent rate only with one,
        # no viscoelastic elasticity yet
        m = mm.validate_material(cp(kin='updated_lagrangian', tangent_rate='jaumann'))
        self.assertEqual(mm.material_components(m)[-1], ('kin', 'updated_lagrangian'))
        self.assertTrue(mm.material_spatial_kin(m))
        self.assertFalse(mm.material_spatial_kin(cp(kin='total_lagrangian')))
        self.assertEqual(mm.material_state_size(m), 7)
        self.assertInvalid([cp(kin='bogus')], "not a registered 'kin' model")
        self.assertInvalid([cp(tangent_rate='jaumann')], 'without a kinematics')
        self.assertInvalid([cp(kin='total_lagrangian', tangent_rate='truesdell')], 'tangent_rate must be')
        self.assertInvalid([cp(kin='total_lagrangian', elas='frozen', visco='maxwell', visco_branches=2)],
                           'not supported yet')
        self.assertInvalid([Material('a', 'cutting_plane', 'linear', 'vonmises', 'swift', kin='total_lagrangian')],
                           'continuum tangent')
        self.assertTrue(mm.material_rotated(cp(kin='updated_lagrangian', yield_='hill48')))
        self.assertFalse(mm.material_rotated(cp(kin='updated_lagrangian')))
        self.assertFalse(mm.material_rotated(cp(kin='total_lagrangian', yield_='hill48')))
        try:
            iso = mm.ELASTICITY_ISOTROPIC
            mm.ELASTICITY_ISOTROPIC = ()
            self.assertInvalid([cp(kin='spatial_logarithmic', yield_='hill48')], 'isotropic elasticity')
        finally:
            mm.ELASTICITY_ISOTROPIC = iso
        # The test registry has a material of every kinematics
        kins = {m.kin for m in mm.load_materials()} - {None}
        self.assertEqual(kins, set(mm.KINEMATICS_MODELS))

    def test_solver_defaults_immutable(self):
        spec = mm.MaterialSolver(families=('elas',), nstate=1)
        with self.assertRaises(TypeError):
            spec.unsupported['elas'] = ('linear',)
        with self.assertRaises(TypeError):
            spec.composes['vp'] = 'hard'

    def test_solver_iteration_limit(self):
        # Default iteration limit of the material type: the historical limit of the solver
        self.assertEqual(mm.MaterialSolver(families=('elas',), nstate=1).iter_nw, 200)
        self.assertEqual(mm.MATERIAL_SOLVERS['closest_point_2'].iter_nw, 200)
        self.assertEqual(mm.MATERIAL_SOLVERS['cutting_plane'].iter_nw, 100)

    def test_doc(self):
        mm.validate_material(cp(doc='x' * mm.MATERIAL_DOC_MAX))
        self.assertInvalid([cp(doc='x' * (mm.MATERIAL_DOC_MAX + 1))], 'doc must be')


class TestNames(unittest.TestCase):

    def test_generated_names(self):
        self.assertEqual(mm.material_module('cp_a'), 'muscle_mat_cp_a')
        self.assertEqual(mm.material_type('cp_a', 'sp'), 'mat_cp_a_sp')
        self.assertEqual(mm.material_proc('update', 'cp_a', 'dp'), 'mat_update_cp_a_dp')
        self.assertEqual(mm.material_nstate('cp_a'), 'MAT_CP_A_NSTATE')
        with self.assertRaises(ValueError):
            mm.material_proc('bogus', 'cp_a', 'dp')

    def test_identifiers(self):
        ids = mm.material_identifiers('cp_a')
        self.assertIn('muscle_mat_cp_a', ids)
        self.assertIn('mat_tangent_cp_a_sp', ids)
        self.assertEqual(len(ids), 2 + len(mm.KINDS) * (1 + len(mm.MATERIAL_PROCS)))

    def test_components(self):
        self.assertEqual(mm.material_components(cp(yield_='hill48', hard='ludwik')),
                         [('elas', 'linear'), ('yield', 'hill48'), ('hard', 'ludwik')])
        jc = Material('a', 'closest_point_time', 'linear', 'vonmises', 'swift', vp='jc')
        self.assertEqual(mm.material_components(jc),
                         [('elas', 'linear'), ('yield', 'vonmises'), ('vp', 'jc'), ('hard', 'swift')])
        rk = Material('a', 'closest_point_time', 'linear', 'vonmises', vp='rk')
        self.assertEqual(mm.material_families(rk), ['elas', 'yield', 'vp'])

    def test_composed_arguments(self):
        jc = Material('a', 'closest_point_time', 'linear', 'vonmises', 'swift', vp='jc')
        self.assertEqual(mm.material_composed(jc, 'vp'), 'hard')
        self.assertIsNone(mm.material_composed(jc, 'yield'))
        self.assertEqual(mm.material_comp_args(jc, 'vp', 'flow_stress', 'ep'), ', hard_stress(mat%hard, ep)')
        self.assertEqual(mm.material_comp_args(jc, 'vp', 'dstress_dep', 'x', obj='m'),
                         ', hard_stress(m%hard, x), hard_dstress_dep(m%hard, x)')
        rk = Material('a', 'closest_point_time', 'linear', 'vonmises', vp='rk')
        self.assertEqual(mm.material_comp_args(rk, 'vp', 'dstress_dep', 'ep'), '')


class TestUserRegistry(unittest.TestCase):

    def write(self, text):
        tmp = tempfile.NamedTemporaryFile('w', suffix='.py', delete=False)
        self.addCleanup(pathlib.Path(tmp.name).unlink)
        tmp.write(textwrap.dedent(text))
        tmp.close()
        return tmp.name

    def test_load_path(self):
        path = self.write("""
            from muscle_meta import Material
            MATERIALS = [Material('steel', 'closest_point', elas='linear', yield_='hill48', hard='vocemod')]
            """)
        self.assertEqual([m.name for m in mm.load_materials(path)], ['steel'])
        self.assertEqual(mm.material('steel', path).yield_, 'hill48')

    def test_no_materials(self):
        path = self.write('X = 1\n')
        with self.assertRaises(ValueError):
            mm.load_materials(path)

    def test_invalid_entry(self):
        path = self.write("""
            from muscle_meta import Material
            MATERIALS = [Material('Steel', 'closest_point', elas='linear', yield_='hill48', hard='vocemod')]
            """)
        with self.assertRaises(ValueError):
            mm.load_materials(path)

    def test_command_line(self):
        # The command CMake runs at configure time (cmake/MuscleFypp.cmake, muscle_add_materials)
        out = subprocess.run([sys.executable, '-B', str(FYPP_DIR / 'muscle_meta.py'), 'materials'],
                             capture_output=True, text=True, check=True).stdout.strip()
        self.assertEqual(out.split(';'), [m.name for m in mm.load_materials()])
        path = self.write("""
            from muscle_meta import Material
            MATERIALS = [Material('a', 'closest_point', elas='frozen', yield_='hill48', hard='vocemod')]
            """)
        res = subprocess.run([sys.executable, '-B', str(FYPP_DIR / 'muscle_meta.py'), 'materials', path],
                             capture_output=True, text=True)
        self.assertEqual(res.returncode, 1)
        self.assertIn("needs a 'visco' model", res.stderr)
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(mm._main(['bogus']), 2)


class TestSolverBodies(unittest.TestCase):
    """The include chain of the material template (fypp takes only literal file names) and the bodies
    on disk match ``MATERIAL_SOLVERS``."""

    ROOT = FYPP_DIR.parents[1]

    def test_every_solver_has_a_body(self):
        for solver in mm.MATERIAL_SOLVERS:
            body = mm.material_solver_body(solver)
            self.assertTrue((self.ROOT / 'src' / 'solvers' / body).is_file(), body)

    def test_template_chain_matches_metadata(self):
        import re
        text = (self.ROOT / 'src' / 'materials' / 'muscle_material.fypp').read_text()
        chain = re.findall(r"^#:(?:el)?if MAT\.solver == '(\w+)'\n+#:include \"\.\./solvers/([\w.]+)\"", text,
                           re.MULTILINE)
        self.assertEqual(chain, [(s, mm.material_solver_body(s)) for s in mm.MATERIAL_SOLVERS])

    def test_unknown_solver(self):
        with self.assertRaises(ValueError):
            mm.material_solver_body('nope')

    def test_command_line(self):
        out = subprocess.run([sys.executable, '-B', str(FYPP_DIR / 'muscle_meta.py'), 'solver-bodies'],
                             capture_output=True, text=True, check=True).stdout.strip()
        self.assertEqual(out.split(';'), [mm.material_solver_body(s) for s in mm.MATERIAL_SOLVERS])


if __name__ == '__main__':
    unittest.main()
