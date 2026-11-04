! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
program muscle_audit_probe
    !! Smallest program with an OpenMP target region: linked with the whole library (target
    !! muscle_audit_probe, src/CMakeLists.txt) so that the offload compiler builds the device image of
    !! every `declare target` procedure, which tools/audit_kernels.sh inspects. Never run.
    implicit none
    integer :: i

    i = 0
    !$omp target map(tofrom: i)
    i = i + 1
    !$omp end target
    if (i /= 1) error stop 1
end program muscle_audit_probe
