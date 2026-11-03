# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Merges the records written by a golden test (``--write``) into an existing reference file, for
regenerate.sh: prints the comment header of the old file (every comment line before its first
record: the capture, the format lines and the approved updates and removals) and the new records,
keeping the old line of a key wherever the new values are equal for the golden comparison (a NaN
matches any NaN, D-018.12f: the sign of a NaN depends on the operand order the compiler chooses).

Usage: python3 -I tests/reference/merge_golden.py OLD NEW > MERGED
"""
import sys


def is_nan(bits):
    """Whether the hexadecimal IEEE-754 double ``bits`` is a NaN."""
    v = int(bits, 16)
    return (v & 0x7FF0000000000000) == 0x7FF0000000000000 and (v & 0x000FFFFFFFFFFFFF) != 0


def same(a, b):
    """Whether two record lines have the same values for the golden comparison."""
    va, vb = a.rsplit(' | ', 1)[1].split(), b.rsplit(' | ', 1)[1].split()
    return len(va) == len(vb) and all(x == y or (is_nan(x) and is_nan(y)) for x, y in zip(va, vb))


def key(line):
    """Key of a record line (everything before the count and the values)."""
    return line.rsplit(' | ', 2)[0]


def main(old_path, new_path):
    old_lines = open(old_path).read().splitlines()
    out = []
    for line in old_lines:
        if line and not line.startswith('#'):
            break
        out.append(line)
    old = {key(line): line for line in old_lines if line and not line.startswith('#')}
    for line in open(new_path).read().splitlines():
        if not line or line.startswith('#'):
            continue
        if key(line) in old and same(old[key(line)], line):
            line = old[key(line)]
        out.append(line)
    print('\n'.join(out))


if __name__ == '__main__':
    main(*sys.argv[1:3])
