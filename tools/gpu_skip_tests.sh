# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
#
# Target tests not run on the Intel GPU (sourced by tools/check.sh and tools/kernel_resources.sh):
# GPU_SKIP_TESTS maps ctest name patterns (regex) to the reason.
# Not run on the Intel GPU (PROBLEMAS.md, GPU / compilers). Most of the closest point, closest point 2,
# closest point time and Perzyna sp kernels excluded before lost the device because IGC reserved too little
# private memory for the run-time-size temporary of COND_INF (D-030, fixed: NORM_INF); they run again. Still out:
# - mat_batch / mat_batch_openmp: IGC miscompiled the closest point kernels with more than one work item
#   (D-019.5c, not checked again after D-030);
# - two sp time-dependent kernels: other iterations than the CPU (and other stresses and tangents in
#   tm_jcfull_vm_vocemod), also with the EU fusion disabled (to investigate);
# - every dp material kernel: the emulated FP64 of the UHD 770 loses the device, not the same from one run to
#   the next.
declare -A GPU_SKIP_TESTS=(
    ['materials-mat_batch(_openmp)?']="IGC miscompiled the closest point kernels with more than one work item (D-019.5c)"
    ['materials-mat_kernels_target_(tm_rk_vm|tm_jcfull_vm_vocemod)_sp']="IGC: other iterations than the CPU (to investigate)"
    ['materials-mat_(kernels_target|batch_dev)_[a-z0-9_]+_dp']="IGC with the emulated FP64: device lost (not reproducible run to run)"
)
