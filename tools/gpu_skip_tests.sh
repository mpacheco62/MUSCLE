# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Target tests not run on the Intel GPU (sourced by tools/check.sh and tools/kernel_resources.sh):
# GPU_SKIP_TESTS maps ctest name patterns (regex) to the reason.
# Not run on the Intel GPU (PROBLEMAS.md, GPU / compilers): IGC (the JIT of the Intel driver) miscompiles
# the closest point return mappings when more than one work item runs (wrong stresses, non-converged
# statuses or ZE_RESULT_ERROR_DEVICE_LOST; right with one work item or with the IGC optimizations off,
# -ze-opt-disable)
# (also the closest point time solver and the dp of the other solvers in mat_kernels_target: device lost or
# other iterates, and not the same from one run to the next; the sp of perzyna, cutting plane and, among the
# time-dependent materials, only tm_jc_cpb06_ludwik pass)
declare -A GPU_SKIP_TESTS=(
    ['materials-mat_batch(_openmp)?']="IGC miscompiles the closest point kernels with more than one work item"
    ['materials-mat_kernels_target_(tm_(mrk|rk|nnl|va)_vm|tm_jcfull_vm_vocemod|cp_dp_ludwik|cp_kelvin_hill48_swift|cp2_maxwell_vm_bilinear)_(sp|dp)']="IGC: closest point / closest point time kernels"
    ['materials-mat_batch_dev_(pz|kin_ul_pz)_vm_bilinear_sp']="IGC: device lost in the Perzyna batch with more than one work item"
    ['materials-mat_(kernels_target|batch_dev)_[a-z0-9_]+_dp']="IGC with the emulated FP64: device lost (not reproducible run to run)"
)
