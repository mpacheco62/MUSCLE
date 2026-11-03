# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Sample inputs of every model for the tests of the material initializers (``mat_init_<name>``,
D-018.17a; tests/materials/test_mat_init.fypp loads this file through ``muscle_meta.load_python_file``).

``sample_inputs(mm, prefix, model)``: ``{input name: literal}`` for the inputs of the model
(``mm.model_inputs``), the literals without kind (a list of literals for an input per branch). The values
are those of the first golden set of the model (golden_params.py: component assignments and keyword or
positional initializer calls), the default of a parameter the golden set does not give, and Young
modulus and Poisson ratio for the linear elasticity.
"""

import os
import re

_HERE = os.path.dirname(os.path.abspath(__file__))
_DP_MODES = {'TENSION': '1', 'COMPRESSION': '2', 'SHEAR': '3'}
_NUM = r'-?[0-9.]+(?:e-?[0-9]+)?'


def _golden(mm):
    return mm.load_python_file(os.path.join(_HERE, 'golden_params.py')).GOLDEN_PARAMS


def _literal(text):
    """Literal without kind of a golden value (``500.0_{k}`` -> ``500.0``, ``DP_HARDENING_SHEAR`` -> ``3``,
    ``[1.0_{k}, 2.0_{k}]`` -> ``['1.0', '2.0']``)."""
    text = text.strip()
    if text.startswith('['):
        return [_literal(t) for t in text[1:-1].split(',')]
    if text.startswith('DP_HARDENING_'):
        return _DP_MODES[text[len('DP_HARDENING_'):]]
    return text.replace('_{k}', '')


def sample_inputs(mm, prefix, model):
    """Sample value of every input of a model (see the module docstring)."""
    inputs = mm.model_inputs(prefix, model)
    if (prefix, model) == ('elas', 'linear'):
        return {'young': '210000.0', 'poisson': '0.3'}
    res = {}
    sets = _golden(mm).get((prefix, model), {})
    stmts = ' '.join(next(iter(sets.values()))) if sets else ''
    stmts = stmts.replace('@@', ' ')
    for name, val in re.findall(r'\{v\}%(\w+) = (' + _NUM + r'_\{k\})', stmts):
        res[name] = _literal(val)
    call = re.search(r'call \w+_init\(\{v\},(.*)\)', stmts)
    if call:
        args = re.findall(r'\s*(?:(\w+)\s*=\s*)?(\[[^\]]*\]|[^,]+)', call.group(1))
        args = [(n, v) for n, v in args if v.strip() and not v.strip().startswith('status')]
        for pos, (name, val) in enumerate(args):
            if name == 'status':
                continue
            res[name or inputs[pos].name] = _literal(val)
    defaults = {p.name: p.default for p in mm.model_param_list(prefix, model)}
    out = {}
    for i in inputs:
        val = res.get(i.name, defaults.get(i.name))
        if val is None:
            raise ValueError(f'no sample value for input {i.name!r} of {prefix} {model}')
        out[i.name] = val
    return out


def prop_literal(mm, mat, prop, kind, settings):
    """Fortran real literal (kind ``kind``) of parameter ``prop`` (``mm.material_props``) of the parameter
    array of a material: the sample value of its input (its branch for an input per branch), or for a
    solver setting ``settings[name]`` (a literal without kind)."""
    if prop.component is None:
        val = settings[prop.name]
    else:
        val = sample_inputs(mm, prop.component, getattr(mat, mm.MATERIAL_FIELDS[prop.component]))[prop.name]
        if prop.branch:
            val = val[prop.branch - 1]
    if not any(c in val for c in '.e'):
        val += '.0'
    return f'{val}_{kind}'
