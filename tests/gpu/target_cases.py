# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Cases of the generated target tests (testing.md, GPU: every kernel in an OpenMP target region
against the CPU, in every kind), one small test executable per case so that each device image is
small (the JIT of the Intel GPU has memory and time caps, D-019).

Read by the templates (``muscle_meta.load_python_file``) and by CMake
(``python target_cases.py <group> <kinds>``, kinds separated by ``;``), which preprocesses the template
of a group once per case with the fypp define ``TARGET_CASE='<case>'``:

- ``tensor``: ``<family>_<dim>_<kind>``, the free functions of muscle_tensor_functions
  (tests/tensors/test_tensor_target.fypp);
- ``family``: ``<prefix>_<kind>``, the generics of every model of a family
  (tests/gpu/test_family_target.fypp; the kinematics are tested by kinematics-kin_kinds);
- ``material``: ``<material>_<kind>``, ``mat_update`` of materials of the test registry that together
  use every model and every solver (tests/materials/test_mat_kernels_target.fypp);
- ``batch_dev``: ``<material>_<kind>``, the batch API on device memory of the materials of
  ``BATCH_DEV_MATERIALS`` (tests/materials/test_mat_batch_dev.fypp).
"""
import sys

import muscle_meta

#: Free functions of muscle_tensor_functions by family besides `is_approx` and `inverse` (every family):
#: (family, functions); also used by tests/tensors/test_muscle_tensor_functions.fypp.
TENSOR_FUNCTIONS = (('2Osym', ('norm', 'det', 'square')), ('2O', ('norm', 'det', 'transpose')),
                    ('4O3sym', ('norm',)), ('4O2sym', ('norm', 'convert_3sym')), ('4O', ('norm', 'convert_2sym')))

#: Model families with a generated target test.
FAMILY_PREFIXES = ('hard', 'vp', 'os', 'yield', 'elas', 'visco')

#: Materials with a target test and the golden case of each (test_mat_<solver>.fypp): (material,
#: elasticity, overstress or viscoplastic law, hardening, criterion, dt, yield strain a_y, first factor
#: f0 of the strain ramp f a_y, f from f0 to 4). f0 = 2 for jcfull: its barely plastic points (rates
#: near 0, where the rate term of the law is singular) are ill-conditioned, they converge to other
#: iterates with other rounding.
MATERIAL_CASES = (
    ('tm_mrk_vm', 'steel', 'mrk', '-', 'vonmises', '1.0e-2', '1.15e-3', '0.0'),
    ('tm_rk_vm', 'steel', 'rk', '-', 'vonmises', '1.0e-3', '2.47e-3', '0.0'),
    ('tm_nnl_vm', 'steel', 'nnl', '-', 'vonmises', '1.0e-2', '3.57e-4', '0.0'),
    ('tm_va_vm', 'steel', 'va', '-', 'vonmises', '1.0e-2', '8.89e-4', '0.0'),
    ('tm_jcfull_vm_vocemod', 'steel', 'jcfull', 'vocemod', 'vonmises', '1.0e-3', '9.16e-4', '2.0'),
    ('tm_jc_cpb06_ludwik', 'steel', 'jc2', 'ludwik', 'cpb06-2', '1.0e-3', '7.04e-4', '0.0'),
    ('pz_cha_yld2004_ludwik', 'steel', 'cha1', 'ludwik', 'yld2004-4', '1.0e-2', '9.21e-4', '0.0'),
    ('pz_cpb06ex2_bilinear', 'steel', 'per1', 'bilinear', 'cpb06ex2-1', '1.0e-2', '8.94e-4', '0.0'),
    ('cp_dp_ludwik', 'steel', '', 'ludwik', 'dp-3', '0.0', '7.39e-4', '0.0'),
    ('cp_kelvin_hill48_swift', 'kelvin', '', 'swift', 'hill48-2', '0.1', '6.55e-3', '0.0'),
    ('cp2_maxwell_vm_bilinear', 'maxwell', '', 'bilinear', 'vonmises', '0.1', '1.56e-2', '0.0'),
    ('cut_yld2004_ludwik', 'steel', '', 'ludwik', 'yld2004-2', '0.0', '8.97e-4', '0.0'),
)

#: Materials of the test of the batch API on device memory (D-019.7a), one of each input and elasticity:
#: strain (cutting plane, Perzyna), frozen Maxwell elasticity, deformation gradients.
BATCH_DEV_MATERIALS = ('cut_vm_bilinear', 'pz_vm_bilinear', 'cut_maxwell_vm_bilinear', 'kin_ul_pz_vm_bilinear')


def cases(group, kinds):
    """Case names of ``group`` for the real kinds ``kinds``."""
    if group == 'tensor':
        return [f'{fam}_{dim}_{k}' for fam, _ in TENSOR_FUNCTIONS for dim in muscle_meta.DIMS for k in kinds]
    if group == 'family':
        return [f'{p}_{k}' for p in FAMILY_PREFIXES for k in kinds]
    if group == 'batch_dev':
        return [f'{m}_{k}' for m in BATCH_DEV_MATERIALS for k in kinds]
    if group == 'material':
        return [f'{c[0]}_{k}' for c in MATERIAL_CASES for k in kinds]
    raise ValueError(f'unknown group of target tests {group!r}')


def _size(typ, dim):
    """Number of `real` values of an argument or result of type ``typ`` (a scalar or a tensor family)."""
    return 1 if typ in ('real', 'integer', 'logical') else len(muscle_meta.layout(typ, dim).comps)


def family_inputs(prefix, model, dim):
    """Inputs of the point procedure of a model in ``test_family_target``: every argument of its functions
    (family arguments, extra arguments, composed values of a composed model, D-017.6b), in order of first
    use, as (``ModelArg``, first position in the input array, size)."""
    fam = muscle_meta.model_family(prefix)
    args = {}
    for f in fam.funcs:
        if f.models is not None and model not in f.models:
            continue
        comp = [c.arg for c in fam.composition.values if c.arg.name in f.comp] if model in fam.composed else []
        for a in list(fam.args) + list(f.extra) + comp:
            args.setdefault(a.name, a)
    res, pos = [], 1
    for a in args.values():
        res.append((a, pos, _size(a.type, dim)))
        pos += _size(a.type, dim)
    return res


def family_outputs(prefix, model, dim):
    """Outputs of the point procedure of a model in ``test_family_target``: (procedure of
    ``model_procs``, the call (``None`` for the joint evaluation and its flag), result type, first
    position in the output array, size)."""
    fam = muscle_meta.model_family(prefix)
    res, pos = [], 1
    for n in muscle_meta.model_procs(prefix, model):
        if fam.joint and n == fam.joint.flag:
            call, typ, size = None, 'logical', 1
        elif fam.joint and n == fam.joint.name:
            call, typ = None, 'joint'
            size = sum(_size(muscle_meta.model_func(prefix, f).result.type, dim) for f in fam.joint.funcs)
        else:
            f = muscle_meta.model_func(prefix, n)
            names = [a.name for a in fam.args] + [a.name for a in f.extra]
            if model in fam.composed:
                names += list(f.comp)
            call, typ = f'{prefix}_{n}(prm, {", ".join(names)})', f.result.type
            size = _size(typ, dim)
        res.append((n, call, typ, pos, size))
        pos += size
    return res


def family_output_blocks(prefix, model, dim):
    """Blocks of the output array of the point procedure of a model in ``test_family_target``, each
    compared on a device relative to its own largest component (a result of the joint evaluation of the
    criteria split into its value, gradient and Hessian): (name, first, last position)."""
    res = []
    for n, _, typ, s, z in family_outputs(prefix, model, dim):
        if typ == 'joint':
            n2 = _size('2Osym', dim)
            res += [(f'{n} f', s, s), (f'{n} 2Osym', s + 1, s + n2), (f'{n} 4O3sym', s + n2 + 1, s + z - 1)]
        else:
            res.append((n, s, s + z - 1))
    return res


def visco_output_blocks():
    """Blocks of the output array of a viscoelastic point in ``test_family_target`` (``VNOUT``, every
    branch up to ``MAX_VISCO_BRANCHES``): (name, first, last position)."""
    nb = muscle_meta.MAX_VISCO_BRANCHES
    sizes = [('status', 1), ('n_branches', 1), ('k_alg', 1), ('g_alg', 1), ('decay', nb), ('weight', nb),
             ('stress', 6), ('dstress_dstrain', 21), ('stress_dev', 6), ('dstress_dstrain_dev', 21),
             ('strain_e', 6)] + [(f'h({i})', 6) for i in range(1, nb + 1)]
    return _blocks(sizes)


def material_output_blocks(name):
    """Blocks of the packed output of a point in ``test_mat_kernels_target`` (stress, tangent and the
    history array of ``mat_pack_<name>``: plastic strain, equivalent plastic strain and, with a
    viscoelastic elasticity, elastic strain and history tensors): (name, first, last position)."""
    mat = muscle_meta.material(name)
    nvisco = (muscle_meta.material_state_size(mat) - 7)//6
    sizes = [('stress', 6), ('tangent', 36), ('strain_p', 6), ('ep', 1)]
    if nvisco > 0:
        sizes += [('strain_e', 6)] + [(f'h({i})', 6) for i in range(1, nvisco)]
    return _blocks(sizes)


def _blocks(sizes):
    """(name, first, last position) of consecutive blocks of the given (name, size)."""
    res, pos = [], 1
    for n, z in sizes:
        res.append((n, pos, pos + z - 1))
        pos += z
    return res


def material_case(name):
    """Row of ``MATERIAL_CASES`` of a material."""
    return next(c for c in MATERIAL_CASES if c[0] == name)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: target_cases.py tensor|family|material|batch_dev KINDS', file=sys.stderr)
        sys.exit(2)
    print(';'.join(cases(sys.argv[1], sys.argv[2].split(';'))))
