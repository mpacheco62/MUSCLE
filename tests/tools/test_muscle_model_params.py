# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Unit tests of the parameters and initializers of the registered models in tools/fypp/muscle_meta.py
(``MODEL_PARAMS``, ``MODEL_INITS``, ``model_inputs``; D-018.17e: single source of the model templates
and of the material initializers).

Run: python3 -m unittest tests/tools/test_muscle_model_params.py
"""

import pathlib
import re
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'fypp'))
import muscle_meta as mm  # noqa: E402


def registered_models():
    """Every ``(prefix, model)`` of the model families."""
    return [(p, m) for p, fam in mm.MODEL_FAMILIES.items() for m in fam.models]


def template(prefix, model):
    """Source of the template of a model (``ModelFamily.module`` or ``MODEL_MODULES``)."""
    name = mm.model_module(prefix, model)
    found = list((ROOT / 'src').rglob(name + '.fypp'))
    assert len(found) == 1, (name, found)
    return found[0].read_text()


class TestModelParams(unittest.TestCase):

    def test_every_registered_model_has_parameters(self):
        self.assertEqual(sorted(mm.MODEL_PARAMS), sorted(registered_models()))

    def test_parameters_are_valid(self):
        for key, params in mm.MODEL_PARAMS.items():
            with self.subTest(model=key):
                self.assertEqual(mm.model_params(params), params)
                self.assertTrue(all(isinstance(p, mm.ModelParam) for p in params))

    def test_examples(self):
        self.assertEqual([p.name for p in mm.model_param_list('hard', 'bilinear')], ['y0', 'K'])
        self.assertEqual(mm.model_param_list('kin', 'corotational'), ())
        with self.assertRaises(ValueError):
            mm.model_param_list('hard', 'nope')

    def test_templates_read_the_single_source(self):
        # Each model template takes its parameters from muscle_meta (no literal list left)
        for p, m in registered_models():
            src = template(p, m)
            with self.subTest(model=(p, m)):
                self.assertNotRegex(src, r"#:set PARAMS = \[|MODEL_PARAMS_TYPE\([^)]*\[")
                self.assertIn(f"model_param_list('{p}', '{m}')", src)


class TestModelInits(unittest.TestCase):

    def test_private_parameters_need_an_initializer(self):
        for (p, m), params in mm.MODEL_PARAMS.items():
            with self.subTest(model=(p, m)):
                if any(q.private for q in params):
                    self.assertIsNotNone(mm.model_init(p, m))

    def test_initializers_are_registered_models(self):
        self.assertTrue(set(mm.MODEL_INITS) <= set(registered_models()))

    def test_inputs(self):
        for p, m in registered_models():
            with self.subTest(model=(p, m)):
                inputs = mm.model_inputs(p, m)
                names = [i.name.lower() for i in inputs]
                self.assertEqual(len(set(names)), len(names))
                for i in inputs:
                    self.assertIsInstance(i, mm.ModelInput)
                    self.assertTrue(i.doc.endswith('.'), i)
                    if mm.model_init(p, m) is not None:
                        self.assertIn(i.type, ('real', 'integer'))
                        self.assertIn(i.shape, ('', '(:)'))
        self.assertEqual([i.name for i in mm.model_inputs('elas', 'linear')], ['young', 'poisson'])
        self.assertEqual([i.name for i in mm.model_inputs('hard', 'swift')], ['k', 'e0', 'n'])
        self.assertEqual(len(mm.model_inputs('yield', 'cpb06')), 14)
        self.assertEqual(len(mm.model_inputs('yield', 'cpb06ex2')), 27)
        self.assertEqual(len(mm.model_inputs('yield', 'yld2004')), 19)
        self.assertEqual([i.shape for i in mm.model_inputs('visco', 'maxwell')], ['', '', '(:)', '(:)'])

    def test_public_parameters_without_initializer(self):
        bad = mm.ModelParam('x', 'real', '', 'X.', None, True)
        mm_params = dict(mm.MODEL_PARAMS)
        mm_params[('hard', 'bilinear')] = (bad,)
        saved = mm.MODEL_PARAMS
        try:
            mm.MODEL_PARAMS = mm_params
            with self.assertRaisesRegex(ValueError, 'no initializer'):
                mm.model_inputs('hard', 'bilinear')
        finally:
            mm.MODEL_PARAMS = saved

    def test_templates_use_the_init_header(self):
        for p, m in mm.MODEL_INITS:
            with self.subTest(model=(p, m)):
                src = template(p, m)
                self.assertIn(f"$:MODEL_INIT('{p}', '{m}', k)", src)
                self.assertNotRegex(src, rf'pure subroutine {p}_{m}_init_\$\{{k\}}\$\(')

    def test_signature(self):
        sig = mm.model_init_signature('visco', 'kelvin', 'dp')
        self.assertEqual(sig, '    pure subroutine visco_kelvin_init_dp(prm, bulk, shear, bulk_viscosity, '
                              'shear_viscosity, status)')
        sig = mm.model_init_signature('yield', 'cpb06ex2', 'dp')
        lines = sig.split('\n')
        self.assertGreater(len(lines), 1)
        self.assertTrue(all(len(line) <= 120 for line in lines))
        self.assertTrue(all(line.endswith(' &') for line in lines[:-1]))
        args = re.sub(r'&\s*', '', sig).split('(', 1)[1].rstrip(')').split(',')
        self.assertEqual([a.strip() for a in args],
                         ['prm'] + [i.name for i in mm.model_inputs('yield', 'cpb06ex2')])


if __name__ == '__main__':
    unittest.main()
