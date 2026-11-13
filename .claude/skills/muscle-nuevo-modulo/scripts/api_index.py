#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Index of the public API of MUSCLE: for every module of src/, its public names (generics, types, constants) with
the first line of their FORD documentation. Search it before writing any function that is not the physics of the
new model: what already exists must be reused, not rewritten.

Usage: api_index.py REPO [--all] [WORD ...]
  REPO: root of the MUSCLE repository.
  --all: also the specific procedures of each kind (by default only generics, types and constants: call the
         generic, never a specific).
  WORD: optional filters (case-insensitive, all must match the module, the name or its doc line).
Templated names are shown with <k> (kind) and <d> (dimension); names produced by Python expressions in the
templates (lists of specifics) are skipped: their generic is listed.
"""
import re
import sys
from pathlib import Path

PUBLIC = re.compile(r'^\s*public\s*::\s*(.*?)\s*$', re.I)
MODULE = re.compile(r'^\s*module\s+(\w+)\s*$', re.I)


def logical_lines(text):
    """Lines with Fortran continuations joined, each with the index of its first physical line."""
    lines = text.split('\n')
    out, buf, start = [], '', 0
    for i, line in enumerate(lines):
        code = line.rstrip()
        if not buf:
            start = i
        if code.endswith('&'):
            buf += code[:-1] + ' '
            continue
        out.append((start, buf + code.lstrip().lstrip('&') if buf else code))
        buf = ''
    return lines, out


def doc_after(lines, i):
    """First FORD doc line right after line i (the declaration), or ''."""
    for j in range(i + 1, min(i + 4, len(lines))):
        s = lines[j].strip()
        if s.startswith('!!'):
            return s[2:].strip()
        if s and not s.startswith('#'):
            break
    return ''


def normalize(name):
    name = name.replace('${k}$', '<k>').replace('${d}$', '<d>')
    return None if '${' in name or '$' in name or not re.match(r'^[A-Za-z][\w<>]*$', name) else name


def find_doc(lines, name):
    """Doc of the declaration of a public name: generic interface, type or parameter."""
    pat = re.escape(name).replace('<k>', r'\$\{k\}\$').replace('<d>', r'\$\{d\}\$')
    decls = [rf'^\s*interface\s+{pat}\s*$', rf'^\s*type\s*(,[^:]*)?::\s*{pat}\s*$',
             rf'^\s*(pure\s+|elemental\s+)*(function|subroutine)\s+{pat}\s*\(', rf'::\s*{pat}\s*(=|$)']
    for d in decls:
        rx = re.compile(d, re.I)
        for i, line in enumerate(lines):
            if rx.search(line):
                doc = doc_after(lines, i)
                if doc:
                    return doc
    return ''


def index(repo):
    rows = []
    for path in sorted(list(Path(repo, 'src').rglob('*.fypp')) + list(Path(repo, 'src').rglob('*.f90'))):
        text = path.read_text(errors='replace')
        lines, logical = logical_lines(text)
        mod = next((m.group(1) for _, l in logical if (m := MODULE.match(l))), None)
        if not mod:
            continue
        names = []
        for _, l in logical:
            m = PUBLIC.match(l)
            if m:
                raws = m.group(1).split(',')
            else:   # a declaration with the public attribute: `real(dp), parameter, public :: TOL_X = ...`
                m = re.match(r'^[^!]*,\s*public\s*(,[^:]*)?::\s*([\w${}]+)', l, re.I)
                raws = [m.group(2)] if m else []
            for raw in raws:
                n = normalize(raw.strip())
                if n and n not in names:
                    names.append(n)
        for n in names:
            doc = find_doc(lines, n)
            if not doc:     # a generic without doc of its own: that of its first specific
                spec = next((m for m in names if m.startswith(n + '_') and '<k>' in m), None)
                doc = find_doc(lines, spec) if spec else ''
            is_type = any(re.search(rf'^\s*type\s*(,[^:]*)?::\s*{re.escape(n).replace("<k>", r"[$][{]k[}][$]")}\s*$',
                                    l, re.I) for l in lines)
            specific = '<k>' in n and not is_type and not n.isupper() and not n.split('_<')[0].isupper()
            rows.append((path.relative_to(repo).parent.as_posix(), mod, n, doc.replace('${k}$', '<k>'), specific))
    return rows


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    repo, args = sys.argv[1], sys.argv[2:]
    show_all = '--all' in args
    words = [w.lower() for w in args if w != '--all']
    last = None
    for area, mod, name, doc, specific in index(repo):
        if specific and not show_all:
            continue
        if words and not all(w in f'{mod} {name} {doc}'.lower() for w in words):
            continue
        if (area, mod) != last:
            print(f'\n{area}/{mod}')
            last = (area, mod)
        print(f'  {name:42s} {doc[:110]}')


if __name__ == '__main__':
    main()
