! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_status
    !! Module muscle_status
    !! ====================
    !! Status codes of the kernel API (gpu-static-dispatch.md): a kernel procedure never stops; it
    !! returns one of these codes in an `integer, intent(out) :: status` argument and the caller (a
    !! host driver, a historical OO adapter) decides what to do. Shared by every model family.
    !!
    !! - `STATUS_OK`: valid arguments, results computed.
    !! - `STATUS_INVALID_PARAMETER`: a model parameter is out of its domain (an initializer
    !!   `<fam>_<model>_init`); the parameter type keeps its default values.
    !! - `STATUS_INVALID_TIME_STEP`: the time step is out of the domain of the model (e.g. dt < 0,
    !!   or dt = 0 in a Kelvin-Voigt solid, D-017.13b); the results keep their default values.
    !! - `STATUS_INVALID_DEFORMATION`: the deformation gradient is not admissible for a finite strain
    !!   kinematics (D-017.14a): det F <= 0 (or NaN), F not invertible, or an eigenvalue of
    !!   C = F^T F at or below the threshold of the kinematics (`KIN_EIG_MIN_<k>`). The kinematic
    !!   state is still computed with the usual formulas (as the historical OO kinematics, which
    !!   ignore this code), but it is not valid (e.g. R not orthogonal, a logarithmic strain
    !!   component dropped) and must not be used.
    !!
    !! Common convention (D-017.14f): when `status` is not `STATUS_OK` the outputs must not be used.
    !! Each family documents whether they keep their default values or are still computed (to keep
    !! the bits of the historical OO API): the viscoelastic family leaves the defaults, the
    !! kinematics leave the computed (invalid) state. A caller (the solvers of Phase 7) always
    !! checks the status.
    !!
    !! The error codes are distinct from the result codes of the historical return mapping solvers
    !! (`STATUS_UNDEFINED = -1`, `STATUS_ELASTIC_CASE = 0`, `STATUS_CONVERGED = 1`, ... in their
    !! own modules), which are rewritten in Phase 7.
    implicit none
    private

    integer, parameter, public :: STATUS_OK = 0
        !! Valid arguments, results computed.
    integer, parameter, public :: STATUS_INVALID_PARAMETER = 101
        !! A model parameter is out of its domain.
    integer, parameter, public :: STATUS_INVALID_TIME_STEP = 102
        !! The time step is out of the domain of the model.
    integer, parameter, public :: STATUS_INVALID_DEFORMATION = 103
        !! The deformation gradient is not admissible (det F <= 0, singular F or a vanishing
        !! eigenvalue of C).

end module muscle_status
