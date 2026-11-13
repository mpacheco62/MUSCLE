# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Unit tests of the parameter array of the declared materials in tools/fypp/muscle_meta.py
(``material_props``, ``material_init_call``, ``material_props_table`` and the command
``muscle_meta.py props NAME``; D-018.17a, e). tests/materials/test_muscle_material_init.fypp checks that the Fortran
initializers ``mat_init_<name>`` follow the same order.

Run: python3 -m unittest tests/tools/test_muscle_material_props.py
"""

import pathlib
import re
import subprocess
import sys
import unittest

FYPP_DIR = pathlib.Path(__file__).resolve().parents[2] / 'tools' / 'fypp'
sys.path.insert(0, str(FYPP_DIR))
import muscle_meta as mm  # noqa: E402

Material = mm.Material


def table_rows(name):
    """Rows ``(index, component, parameter)`` of the table printed by ``muscle_meta.py props name``."""
    out = subprocess.run([sys.executable, '-B', str(FYPP_DIR / 'muscle_meta.py'), 'props', name],
                         check=True, capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        m = re.match(r'\| (\d+) \| (\w+)[^|]* \| `([^`]+)` \|', line)
        if m:
            rows.append((int(m.group(1)), m.group(2), m.group(3)))
    return out, rows


class TestMaterialProps(unittest.TestCase):

    def test_example(self):
        mat = Material('a', 'closest_point', 'linear', 'vonmises', 'bilinear')
        self.assertEqual([p.name for p in mm.material_props(mat)], ['young', 'poisson', 'y0', 'K', 'tol', 'iter_nw'])
        self.assertEqual([p.index for p in mm.material_props(mat)], [1, 2, 3, 4, 5, 6])
        self.assertEqual(mm.material_nprops('cp_vm_swift'), 'MAT_CP_VM_SWIFT_NPROPS')

    def test_maxwell_branches_interleaved(self):
        mat = mm.material('cp_maxwell_vm_bilinear')
        props = [(p.component, p.name, p.branch) for p in mm.material_props(mat)][:6]
        self.assertEqual(props, [('visco', 'bulk', 0), ('visco', 'shear_inf', 0), ('visco', 'shear', 1),
                                 ('visco', 'tau', 1), ('visco', 'shear', 2), ('visco', 'tau', 2)])
        call = mm.material_init_call(mat, 'visco', 'dp')
        self.assertIn('shear=props(3:5:2), tau=props(4:6:2)', call)
        self.assertIn('status=status', call)

    def test_integer_inputs_and_settings(self):
        mat = mm.material('cp_dp_ludwik')
        props = mm.material_props(mat)
        mode = [p for p in props if p.name == 'hardening_mode'][0]
        self.assertEqual(mode.type, 'integer')
        self.assertEqual(mm.material_prop_value(mode), f'nint(props({mode.index}))')
        self.assertEqual([(p.component, p.name, p.type) for p in props[-2:]],
                         [(None, 'tol', 'real'), (None, 'iter_nw', 'integer')])
        self.assertIn('hardening_mode=nint(props(5))', mm.material_init_call(mat, 'yield', 'dp'))

    def test_domains(self):
        props = {p.name: p for p in mm.material_props(mm.material('cp_dp_ludwik'))}
        self.assertEqual(mm.material_prop_check(props['young'], 'dp'), 'props(1) > 0.0_dp')
        self.assertEqual(mm.material_prop_check(props['poisson'], 'sp'), 'props(2) > -1.0_sp .and. props(2) < 0.5_sp')
        self.assertEqual(mm.material_prop_check(props['hardening_mode'], 'dp'),
                         'nint(props(5)) >= 1 .and. nint(props(5)) <= 3')
        self.assertEqual(mm.material_prop_check(props['tol'], 'dp'), 'props(9) > 0.0_dp')
        self.assertEqual(mm.material_prop_check(props['iter_nw'], 'dp'), 'nint(props(10)) > 0')
        self.assertIsNone(mm.material_prop_check(props['sigma0'], 'dp'))

    def test_model_domains(self):
        # Johnson-Cook (D-018.19c): a condition on several inputs, with their positions in the array
        mat = mm.material('tm_jc_vm_swift')
        self.assertEqual([c for c, _ in mm.material_domain_checks(mat, 'sp')],
                         ['props(3)*abs(log(props(5)/props(4))) < 1.0_sp'])
        self.assertEqual(mm.material_domain_checks(mm.material('tm_jcfull_vm_vocemod'), 'dp'), [])
        saved = mm.MODEL_DOMAINS
        try:
            mm.MODEL_DOMAINS = {('vp', 'jc'): (('{nope} > 0.0_{k}', 'X.'),)}
            with self.assertRaisesRegex(ValueError, 'not inputs of the model'):
                mm.material_domain_checks(mat, 'dp')
        finally:
            mm.MODEL_DOMAINS = saved

    def test_component_order_and_no_inputs(self):
        # Components in the order of the material type; kinematics and von Mises have no inputs
        for mat in mm.load_materials():
            with self.subTest(material=mat.name):
                comps = [p for p, _ in mm.material_components(mat)]
                seen = [p.component for p in mm.material_props(mat) if p.component]
                self.assertEqual(sorted(set(seen), key=comps.index), [c for c in comps if c in seen])
                self.assertNotIn('kin', seen)

    def test_long_init_call_is_wrapped(self):
        call = mm.material_init_call(mm.material('cp_cpb06ex2_swift'), 'yield', 'sp')
        lines = call.split('\n')
        self.assertGreater(len(lines), 2)
        self.assertTrue(all(len(line) <= 120 for line in lines))
        self.assertTrue(all(line.endswith(' &') for line in lines[:-1]))
        keywords = re.findall(r'(\w+)=props', call)
        self.assertEqual(keywords, [i.name for i in mm.model_inputs('yield', 'cpb06ex2')])

    def test_unsupported_input_is_rejected(self):
        saved = mm.MODEL_PARAMS
        try:
            params = dict(saved)
            params[('hard', 'bilinear')] = (mm.ModelParam('y', 'real', '(3)', 'Y.', None),)
            mm.MODEL_PARAMS = params
            with self.assertRaisesRegex(ValueError, 'fixed shape'):
                mm.material_props(Material('a', 'closest_point', 'linear', 'vonmises', 'bilinear'))
            params[('hard', 'bilinear')] = (mm.ModelParam('y', 'logical', '', 'Y.', None),)
            with self.assertRaisesRegex(ValueError, 'only real and integer'):
                mm.validate_material(Material('a', 'closest_point', 'linear', 'vonmises', 'bilinear'))
        finally:
            mm.MODEL_PARAMS = saved


class TestPropsTable(unittest.TestCase):

    def test_table_equals_metadata_for_every_material(self):
        for mat in mm.load_materials():
            with self.subTest(material=mat.name):
                out, rows = table_rows(mat.name)
                expected = [(p.index, p.component or 'solver', p.name + (f'({p.branch})' if p.branch else ''))
                            for p in mm.material_props(mat)]
                self.assertEqual(rows, expected)
                self.assertIn(f'`{mm.material_nprops(mat.name)}` = {len(expected)}', out)

    def test_usage(self):
        res = subprocess.run([sys.executable, '-B', str(FYPP_DIR / 'muscle_meta.py'), 'props'],
                             capture_output=True, text=True)
        self.assertEqual(res.returncode, 2)
        res = subprocess.run([sys.executable, '-B', str(FYPP_DIR / 'muscle_meta.py'), 'props', 'nope'],
                             capture_output=True, text=True)
        self.assertEqual(res.returncode, 1)


if __name__ == '__main__':
    unittest.main()
