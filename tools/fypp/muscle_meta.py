# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Python helpers for the MUSCLE FYPP templates.

fypp loads this module with ``-M tools/fypp -m muscle_meta``; ``include/muscle_common.fypp``
re-exports the helpers under short names. Tensor layout metadata is added in Phase 1
(see ``docs/fypp/PLAN.md``).
"""

#: Real kinds known to MUSCLE (Fortran constants ``sp`` and ``dp``).
KINDS = ('sp', 'dp')


def tname(base, kind):
    """Name of the concrete tensor type: ``tname('3D2Osym', 'dp')`` -> ``'ten_3D2Osym_dp'``."""
    if kind not in KINDS:
        raise ValueError(f"unknown kind '{kind}', expected one of {KINDS}")
    return f"ten_{base}_{kind}"
