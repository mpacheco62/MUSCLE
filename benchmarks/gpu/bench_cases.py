# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Cases of the batch benchmark GPU vs CPU (benchmarks/gpu/bench_batch_gpu.fypp), one executable per case
``<material>_<kind>`` so that each device image is small (the JIT of the Intel GPU has memory and time
caps, D-019.5d).

Read by the template (``muscle_meta.load_python_file``) and by CMake (``python bench_cases.py <kinds>``,
kinds separated by ``;``).
"""
import sys

#: Materials of the test registry (tools/fypp/muscle_materials.py) whose batch kernels run right on the
#: Intel UHD 770 in sp in the target tests (tools/gpu_skip_tests.sh), one per solver family that IGC
#: compiles right there: (material, time step). The benchmark itself runs right on the UHD 770 in Release
#: only for cut_vm_bilinear and cut_yld2004_ludwik: pz_cpb06ex2_bilinear and tm_jc_cpb06_ludwik only in
#: Debug, pz_cha_yld2004_ludwik not at all. Strain input, no viscoelasticity (the benchmark checks it).
#: The parameters are the sample
#: inputs of tests/materials/props_samples.py and the time steps those of tests/gpu/target_cases.py.
BENCH_MATERIALS = (
    ('cut_vm_bilinear', '0.0'),
    ('cut_yld2004_ludwik', '0.0'),
    ('pz_cha_yld2004_ludwik', '1.0e-2'),
    ('pz_cpb06ex2_bilinear', '1.0e-2'),
    ('tm_jc_cpb06_ludwik', '1.0e-3'),
)


def cases(kinds):
    """Case names for the real kinds ``kinds``."""
    return [f'{m}_{k}' for m, _ in BENCH_MATERIALS for k in kinds]


if __name__ == '__main__':
    print(';'.join(cases(sys.argv[1].split(';'))))
