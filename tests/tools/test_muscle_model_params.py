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


def type_components(block):
    """Components of a hand-written derived type: ``name``, ``type``, ``shape``, ``default`` (lower case,
    without blanks; continuation lines joined)."""
    text = re.sub(r"&\s*\n\s*&?", "", block)
    out = []
    for line in text.split('\n'):
        found = re.match(r"\s*(real\(\$\{k\}\$\)|integer|type\([^)]*\))(?:\s*,\s*\w+)*\s*::\s*(\w+)"
                         r"(\([^)]*\))?(?:\s*=\s*(.*?))?\s*$", line, re.I)
        if found and not line.lstrip().startswith('!'):
            typ, name, shape, default = found.groups()
            out.append({'name': name.lower(), 'type': typ.replace(' ', ''),
                        'shape': (shape or '').replace(' ', '').lower(),
                        'default': None if default is None else default.replace(' ', '')})
    return out


def expected_type(q):
    """Declared type of a ``ModelParam`` as written in the templates."""
    if q.type == 'real':
        return 'real(${k}$)'
    if q.type == 'integer':
        return 'integer'
    return f'type(ten_{q.type}_${{k}}$)'


def expected_default(q):
    """Default value of a ``ModelParam`` as written in the templates (``None``: no default)."""
    if q.default is None:
        return None
    if q.type == 'real':
        return f'{q.default}_${{k}}$' if re.fullmatch(r'[\d.]+(e[-+]?\d+)?', q.default) else q.default
    if q.type == 'integer':
        return q.default
    zero = f'{q.default}_${{k}}$'
    return f'ten_{q.type}_${{k}}$([{",".join([zero] * 6)}])'


def family_module_source(fam):
    """Source of the module of a family (``muscle_hard_laws``...)."""
    found = list((ROOT / 'src').rglob(fam.family_module + '.fypp'))
    assert len(found) == 1, (fam.family_module, found)
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
                if f"model_param_list('{p}', '{m}')" in src:
                    continue
                # Rewritten by hand (D-020): the type is written in the template and the
                # registered names must be its components, in the registered order.
                block = re.search(rf"type :: {p}_{m}_\$\{{k\}}\$\n(.*?)end type", src, re.S)
                self.assertIsNotNone(block, 'neither model_param_list nor a hand-written type')
                comps = type_components(block.group(1))
                params = mm.model_param_list(p, m)
                self.assertEqual([c['name'] for c in comps], [q.name.lower() for q in params])
                for c, q in zip(comps, params):
                    with self.subTest(component=q.name):
                        self.assertEqual(c['type'], expected_type(q), 'type')
                        self.assertEqual(c['shape'], q.shape.replace(' ', '').lower(), 'shape')
                        self.assertEqual(c['default'], expected_default(q), 'default value')

    def test_hand_written_lists_name_every_registered_model(self):
        # The module of a family lists its models by hand (use, interface, public): a registered model
        # that a list misses would not be reachable through the generics.
        for p, fam in mm.MODEL_FAMILIES.items():
            src = family_module_source(fam)
            for m in fam.models:
                with self.subTest(family=p, model=m):
                    self.assertRegex(src, rf"\b{p}_{m}_\$\{{k\}}\$|\bmuscle_{p}_{m}\b",
                                     'the family module does not mention the model')

    def test_max_visco_branches_literal_matches_the_metadata(self):
        src = (ROOT / 'src' / 'viscoelasticity' / 'muscle_visco_types.fypp').read_text()
        found = re.search(r"parameter\s*(?:,\s*public\s*)?::\s*MAX_VISCO_BRANCHES\s*=\s*(\d+)", src)
        self.assertIsNotNone(found)
        self.assertEqual(int(found.group(1)), mm.MAX_VISCO_BRANCHES)


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


if __name__ == '__main__':
    unittest.main()
