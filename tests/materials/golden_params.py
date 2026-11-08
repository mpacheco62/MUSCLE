# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Parameter sets of the golden cases of the return mapping solvers by name, for the tests of the
declared materials (tests/materials/mat_golden_common.fypp loads this file through
``muscle_meta.load_python_file``).

The sets are those of the historical OO golden tests of the solvers (solvers_golden_common.f90,
test_perzyna_golden.f90, test_time_golden.f90; removed in Phase 7 step 9), plus those of the physical checks.

``GOLDEN_PARAMS[(family prefix, model)][golden name]``: the Fortran statements that set the
parameters of a model to that set, with ``{v}`` = the parameter variable, ``{k}`` = the kind and
``@@`` = a Fortran line break (continuation). An initializer call passes ``status=status`` (a local
integer of the caller) and is followed by its check (``_checked``).
"""
import re

# Long initializer calls, split in Fortran lines with '@@'
_CPB06 = ('call yield_cpb06_init({v}, c11=1.0_{k}, c12=-2.373_{k}, c13=-2.364_{k}, c21=-2.373_{k}, @@'
          'c22=-1.838_{k}, c23=1.196_{k}, c31=-2.364_{k}, c32=1.196_{k}, c33=-2.444_{k}, c44=3.607_{k}, @@'
          'c55=3.607_{k}, c66=3.607_{k}, k=-0.136_{k}, a={a}_{k}, status=status)')
_CPB06EX2 = ('call yield_cpb06ex2_init({v}, 1.283_{k}, 0.208_{k}, -0.421_{k}, 0.208_{k}, 0.535_{k}, @@'
             '-0.791_{k}, -0.421_{k}, -0.791_{k}, 0.746_{k}, 1.156_{k}, 1.418_{k}, 1.329_{k}, 0.030_{k}, @@'
             '1.283_{k}, 0.380_{k}, -0.001_{k}, 0.380_{k}, 2.145_{k}, 0.792_{k}, -0.001_{k}, 0.792_{k}, @@'
             '0.218_{k}, 1.721_{k}, 1.418_{k}, 1.329_{k}, -0.027_{k}, 5.0_{k}, status=status)')
_CPB06_3 = ('call yield_cpb06_init({v}, c11=1.0_{k}, c12=-2.373_{k}, c13=-2.364_{k}, c21=-1.5_{k}, @@'
           'c22=-1.838_{k}, c23=1.196_{k}, c31=-2.364_{k}, c32=0.8_{k}, c33=-2.444_{k}, c44=3.607_{k}, @@'
           'c55=2.9_{k}, c66=4.1_{k}, k=0.25_{k}, a=2.0_{k}, status=status)')
_CPB06EX2_3 = ('call yield_cpb06ex2_init({v}, 0.771_{k}, -0.914_{k}, -1.312_{k}, -0.914_{k}, 0.654_{k}, @@'
               '-1.209_{k}, -1.312_{k}, -1.209_{k}, -0.293_{k}, 1.617_{k}, 1.0_{k}, 1.0_{k}, 0.295_{k}, @@'
               '0.771_{k}, 0.337_{k}, -0.764_{k}, 0.337_{k}, 0.646_{k}, -0.401_{k}, -0.764_{k}, -0.401_{k}, @@'
               '-0.022_{k}, 0.499_{k}, 1.0_{k}, 1.0_{k}, 0.900_{k}, 4.0_{k}, status=status)')
_CPB06_LI = ('call yield_cpb06_init({v}, c11=1.0_{k}, c12=-2.428_{k}, c13=-2.920_{k}, c21=-2.428_{k}, @@'
             'c22=1.652_{k}, c23=-2.236_{k}, c31=-2.920_{k}, c32=-2.236_{k}, c33=1.003_{k}, c44=-3.996_{k}, @@'
             'c55=-3.996_{k}, c66=-3.996_{k}, k=-0.165_{k}, a=2.0_{k}, status=status)')
_YLD2004 = ('call yield_yld2004_init({v}, cp12=-0.069888_{k}, cp13=0.936408_{k}, cp21=0.079143_{k}, @@'
            'cp23=1.003060_{k}, cp31=0.524741_{k}, cp32=1.363180_{k}, cp44=0.954322_{k}, cp55=1.023770_{k}, @@'
            'cp66=1.069060_{k}, cpp12=0.981171_{k}, cpp13=0.476741_{k}, cpp21=0.575316_{k}, @@'
            'cpp23=0.866827_{k}, cpp31=1.145010_{k}, cpp32=-0.079294_{k}, cpp44=1.404620_{k}, @@'
            'cpp55=1.051660_{k}, cpp66=1.147100_{k}, a={a}_{k}, status=status)')
_DP = ('call yield_druckerprager_init({v}, beta_deg={beta}_{k}, K={K}_{k}, @@'
       'hardening_mode=DP_HARDENING_{mode}, status=status)')


def _checked(call):
    """Statements of an initializer call (with ``status=status``) followed by the check of its status."""
    name = re.match(r'call (\w+)\(', call).group(1)
    return [call, f"if (status /= STATUS_OK) error stop '{name} failed'"]


def _sub(text, **kw):
    """``text`` with the given placeholders replaced (the others, ``{v}`` and ``{k}``, are kept)."""
    for key, val in kw.items():
        text = text.replace('{' + key + '}', val)
    return text


GOLDEN_PARAMS = {
    ('hard', 'bilinear'): {'bilinear': ['{v}%y0 = 200.0_{k}', '{v}%K = 1000.0_{k}']},
    ('hard', 'swift'): {'swift': ['{v}%k = 500.0_{k}', '{v}%e0 = 1.0e-3_{k}', '{v}%n = 0.2_{k}'],
                        'swift_small': ['{v}%k = 100.0_{k}', '{v}%e0 = 1.0e-4_{k}', '{v}%n = 0.1_{k}'],
                        # large increments of test_mat_cp_anisotropic (historical OO test)
                        'swift_li': ['{v}%k = 500.0_{k}', '{v}%e0 = 0.01_{k}', '{v}%n = 0.2_{k}']},
    ('hard', 'ludwik'): {'ludwik': ['{v}%sigma0 = 200.0_{k}', '{v}%k = 500.0_{k}', '{v}%n = 0.3_{k}']},
    ('hard', 'vocemod'): {'vocemod': ['{v}%Sy = 272.075784_{k}', '{v}%k = 0.04230655_{k}',
                                      '{v}%q = 297.603923_{k}', '{v}%n = 1.48375781_{k}']},
    ('yield', 'vonmises'): {'vonmises': []},
    ('yield', 'hill48'): {
        'hill48-2': ['{v}%f = 0.25_{k}', '{v}%g = 0.4_{k}', '{v}%h = 0.6_{k}', '{v}%l = 1.3_{k}',
                     '{v}%m = 1.7_{k}', '{v}%n = 1.45_{k}'],
        'hill48-3': ['{v}%f = 0.7_{k}', '{v}%g = 0.2_{k}', '{v}%h = 0.35_{k}', '{v}%l = 2.0_{k}',
                     '{v}%m = 0.9_{k}', '{v}%n = 1.1_{k}'],
        # From the Lankford coefficients r0 = 0.2115, r45 = 1.5769, r90 = 0.6923 (test_mat_cp_anisotropic)
        'hill48-r': ['{v}%f = 0.2115_{k}/(0.6923_{k}*(1.0_{k} + 0.2115_{k}))', '{v}%g = 1.0_{k}/(1.0_{k} + 0.2115_{k})',
                     '{v}%h = 0.2115_{k}/(1.0_{k} + 0.2115_{k})', '{v}%l = 1.25_{k}', '{v}%m = 1.75_{k}',
                     '{v}%n = (0.2115_{k} + 0.6923_{k})*(1.0_{k} + 2.0_{k}*1.5769_{k})/(2.0_{k}*0.6923_{k}*(1.0_{k} + 0.2115_{k}))']},
    ('yield', 'druckerprager'): {
        'dp-1': _checked(_sub(_DP, beta='20.0', K='0.8', mode='TENSION')),
        'dp-2': _checked(_sub(_DP, beta='20.0', K='0.8', mode='COMPRESSION')),
        'dp-3': _checked(_sub(_DP, beta='30.0', K='1.0', mode='SHEAR'))},
    ('yield', 'cpb06'): {'cpb06-2': _checked(_sub(_CPB06, a='2.0')), 'cpb06-3': _checked(_CPB06_3),
                         'cpb06-4': _checked(_sub(_CPB06, a='8.0')),
                         # Shear factors 3.607, 2.9, 4.1 (shears of test_mat_cp_anisotropic)
                         'cpb06-2s': _checked(_sub(_CPB06, a='2.0').replace('c55=3.607_{k}, c66=3.607_{k}',
                                                                          'c55=2.9_{k}, c66=4.1_{k}')),
                         # Ti-6Al-4V at Wp = 48.66 (Tuninetti et al. 2013, Table 2; large increments)
                         'cpb06-li': _checked(_CPB06_LI)},
    ('yield', 'cpb06ex2'): {'cpb06ex2-1': _checked(_CPB06EX2), 'cpb06ex2-3': _checked(_CPB06EX2_3)},
    ('yield', 'yld2004'): {'yld2004-2': _checked(_sub(_YLD2004, a='8.0')), 'yld2004-4': _checked(_sub(_YLD2004, a='6.0'))},
    ('os', 'perzyna'): {'per1': ['{v}%eta = 1.0e-2_{k}', '{v}%n = 1.0_{k}'],
                        'per2': ['{v}%eta = 1.0_{k}', '{v}%n = 0.5_{k}'],
                        'per3': ['{v}%eta = 1.0e-4_{k}', '{v}%n = 5.0_{k}']},
    ('os', 'chaboche'): {'cha1': ['{v}%k = 80.0_{k}', '{v}%n = 2.3_{k}'],
                         'cha2': ['{v}%k = 50.0_{k}', '{v}%n = 1.0_{k}'],
                         'cha3': ['{v}%k = 200.0_{k}', '{v}%n = 0.5_{k}']},
    ('vp', 'jc'): {'jc1': ['{v}%C = 0.1_{k}', '{v}%epdmax = 1.0e-3_{k}'],
                   'jc2': ['{v}%C = 0.014_{k}', '{v}%epdmax = 1.0_{k}'],
                   'jc3': ['{v}%C = 0.05_{k}', '{v}%epdmax = 1.0e-3_{k}', '{v}%epmin = 1.0e-5_{k}']},
    ('vp', 'jcfull'): {'jcfull': ['{v}%C = 0.014_{k}', '{v}%epdmax = 1.0_{k}', '{v}%m = 1.03_{k}',
                                  '{v}%Troom = 293.15_{k}', '{v}%Tmelt = 1793.0_{k}', '{v}%T_current = 700.0_{k}']},
    ('vp', 'mrk'): {'mrk': ['{v}%B01 = 100.020701_{k}', '{v}%B02 = 100.0_{k}', '{v}%epdmax = 0.01_{k}',
                            '{v}%nu1 = 0.01439459_{k}', '{v}%nu2 = 0.14313556_{k}', '{v}%n0 = 3.95295384_{k}',
                            '{v}%D2 = 5.6375e-12_{k}', '{v}%epdmin = 1.0e-5_{k}', '{v}%sig_u = 253.772550_{k}',
                            '{v}%chi1 = 9.1503e-4_{k}', '{v}%chi2 = 0.01942062_{k}']},
    ('vp', 'rk'): {'rk': ['{v}%B0 = 600.0_{k}', '{v}%epdmax = 1.0e7_{k}', '{v}%nu = 0.05_{k}', '{v}%ep0 = 0.01_{k}',
                          '{v}%n0 = 0.3_{k}', '{v}%D2 = 0.25_{k}', '{v}%D1 = 0.2_{k}', '{v}%epdmin = 1.0e-5_{k}',
                          '{v}%sig0 = 200.0_{k}', '{v}%m = 2.0_{k}']},
    ('vp', 'nnl'): {'nnl': ['{v}%sig_a = 334.776175_{k}', '{v}%sig_0 = 5.88318003_{k}', '{v}%KG0 = 0.21906996_{k}',
                            '{v}%n1 = 0.06301105_{k}', '{v}%epd0 = 0.01425422_{k}', '{v}%at = 22.8752986_{k}',
                            '{v}%n0 = 0.97427303_{k}', '{v}%q = 2.0_{k}', '{v}%p = 0.6666667_{k}']},
    # Viscoelastic models of the frozen elasticity (elas_init of the historical solvers_golden_common.f90); {v} is a
    # component of the material, `status` a local integer of the caller
    ('visco', 'maxwell'): {'maxwell': ['call visco_maxwell_init({v}, bulk=20000.0_{k}, shear_inf=1000.0_{k}, @@'
                                       'shear=[2000.0_{k}, 3000.0_{k}], tau=[0.05_{k}, 1.0_{k}], status=status)',
                                       "if (status /= STATUS_OK) error stop 'visco_maxwell_init failed'"]},
    ('visco', 'kelvin'): {'kelvin': ['call visco_kelvin_init({v}, bulk=20000.0_{k}, shear=5000.0_{k}, @@'
                                     'bulk_viscosity=500.0_{k}, shear_viscosity=200.0_{k}, status=status)',
                                     "if (status /= STATUS_OK) error stop 'visco_kelvin_init failed'"]},
    ('vp', 'va'): {'va': ['{v}%B = 374.525856_{k}', '{v}%B1 = 0.00331794_{k}', '{v}%B2 = 0.50280316_{k}',
                          '{v}%n = 0.26412620_{k}', '{v}%m = 0.03766668_{k}', '{v}%sig_u = 195.233128_{k}']},
}

#: Initializers and constants a test imports for the yield criteria of its materials (module, names).
YIELD_USES = {
    'druckerprager': ('muscle_yield_druckerprager',
                      ('yield_druckerprager_init', 'DP_HARDENING_TENSION', 'DP_HARDENING_COMPRESSION',
                       'DP_HARDENING_SHEAR')),
    'cpb06': ('muscle_yield_cpb06', ('yield_cpb06_init',)),
    'cpb06ex2': ('muscle_yield_cpb06ex2', ('yield_cpb06ex2_init',)),
    'yld2004': ('muscle_yield_yld2004', ('yield_yld2004_init',)),
}

#: Initializers a test imports for the viscoelastic models of its materials (module, names).
VISCO_USES = {
    'maxwell': ('muscle_visco_maxwell', ('visco_maxwell_init',)),
    'kelvin': ('muscle_visco_kelvin', ('visco_kelvin_init',)),
}
