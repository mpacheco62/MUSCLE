#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Review aid for a new or changed module of MUSCLE: where its code lives and whether it rewrites what exists.

Usage: check_new_module.py REPO [FILE ...]
  REPO: root of the MUSCLE repository. FILE: sources relative to REPO; by default the .fypp/.f90 files of src/ and
  include/ changed or new in the working tree (git status).

Reports, per file (advisory: every finding is to be justified or fixed, nothing here fails a build):
  1. location: the module is named like its file, and a model module lives in the folder of its family;
  2. procedures outside the contract of a model module (helpers): a helper that does not use the parameters of
     the model is not part of the model: reuse an existing routine (tools: api_index.py) or move it to its area
     (src/math, src/tensors, a support module of the family);
  3. suspicious patterns of code that already exists in the library (trace or deviator by components,
     determinant or eigenvalues written out, magic literals that should be named constants);
  4. blocks repeated between the file and the rest of the library (tools/dup_report.py).
"""
import re
import subprocess
import sys
from pathlib import Path

#: Folder of the model modules of each family prefix (pages/contributing.md, table of the families)
FAMILY_DIRS = {
    'hard': ['src/plasticity/hardening_laws'], 'vp': ['src/plasticity/viscoplastic_laws'],
    'os': ['src/plasticity/viscoplastic_laws'], 'yield': ['src/plasticity/yield_criteria'],
    'elas': ['src/elasticity', 'src/viscoelasticity'], 'visco': ['src/viscoelasticity'],
    'kin': ['src/kinematics'],
}
#: Support modules of a family (not models): no contract to check
SUPPORT = {'muscle_visco_types', 'muscle_kin_finite', 'muscle_yield_anisotropic'}

#: (regex on a code line, what to use instead)
PATTERNS = [
    (r'vals\(1\)\s*\+\s*[\w%]*vals\(2\)\s*\+\s*[\w%]*vals\(3\)',
     'trace or deviator by components: `a .ddot. iden_2O` (trace) and `.dev.` exist in muscle_tensors'),
    (r'(\w+)\(1,\s*1\)\s*\*\s*\1\(2,\s*2\)\s*\*\s*\1\(3,\s*3\)',
     'determinant written out: `det` (muscle_tensor_functions) or the inverses of muscle_math_inverses'),
    (r'\b(acos|atan2)\s*\(|jacobi',
     'eigenvalues by hand: `spectral_decomposition`, `eigenprojections` (muscle_math_eigen, D-032), derivatives '
     'in muscle_math_spectral_derivs'),
    (r'\bmatmul\s*\(.*\bmatmul\s*\(', 'nested matmul: a tensor operator (.ddot., .transform., .tdot.) may exist'),
    (r'\b(sum|maxval|minval|product|any|all|count)\s*\(.*\bdim\s*=',
     'reduction with dim=: in device code a run-time-size temporary overflows the private memory (D-030); '
     'explicit loops, or NORM_INF (muscle_math_inverses)'),
    (r'(?<![\w.])\d+\.\d*e-\d+_', 'small literal threshold: a named constant (muscle_kinds per kind, or a documented '
     'parameter of the module)'),
    (r'\b(epsilon|tiny)\s*\(', 'precision threshold inline: a named constant of muscle_kinds?'),
]

PROC = re.compile(r'^\s*(?:pure\s+|elemental\s+|impure\s+|recursive\s+)*(function|subroutine)\s+([\w${}]+)', re.I)


def changed_files(repo):
    out = subprocess.run(['git', '-C', repo, 'status', '--porcelain'], capture_output=True, text=True).stdout
    files = []
    for line in out.splitlines():
        path = line[3:].split(' -> ')[-1]
        if re.match(r'^(src|include)/.*\.(fypp|f90)$', path):
            files.append(path)
    return files


def code_lines(text):
    for i, line in enumerate(text.split('\n'), 1):
        s = line.strip()
        if not s or s.startswith(('!', '#')):
            continue
        yield i, re.sub(r'!.*$', '', line)


def check(repo, rel):
    path = Path(repo, rel)
    text = path.read_text(errors='replace')
    findings = []
    m = re.search(r'^\s*module\s+(\w+)\s*$', text, re.M | re.I)
    module = m.group(1) if m else None
    if module and module != path.stem:
        findings.append(f'location: module `{module}` in a file named `{path.name}`')
    fam = re.match(r'muscle_([a-z]+)_(\w+)$', path.stem)
    prefix, model = (fam.group(1), fam.group(2)) if fam else (None, None)
    is_model = prefix in FAMILY_DIRS and model != 'laws' and path.stem not in SUPPORT
    if is_model:
        if path.parent.as_posix().replace(str(Path(repo)) + '/', '') not in FAMILY_DIRS[prefix] and \
                rel.rsplit('/', 1)[0] not in FAMILY_DIRS[prefix]:
            findings.append(f'location: a `{prefix}` model outside {FAMILY_DIRS[prefix]}')
        contract = re.compile(rf'^{prefix}_(\w+_)?{model}(_[23]D)?_\$\{{k\}}\$$|^{prefix}_{model}_init(_\$\{{k\}}\$)?$')
        for i, line in enumerate(text.split('\n'), 1):
            p = PROC.match(line)
            if p and not line.strip().lower().startswith('end') and not contract.match(p.group(2)):
                if p.group(2).startswith(f'{prefix}_{model}_'):
                    findings.append(f'model helper (line {i}): `{p.group(2)}`: fine if it uses the parameters of the '
                                    'model; generic algebra in it belongs to src/math or src/tensors')
                else:
                    findings.append(f'HELPER OUTSIDE THE MODEL (line {i}): `{p.group(2)}` is neither a function of the '
                                    'family nor named after the model: reuse an existing routine (api_index.py) or '
                                    'move it to its area (src/math, src/tensors, a support module of the family)')
    for i, code in code_lines(text):
        for rx, advice in PATTERNS:
            if re.search(rx, code, re.I) and 'parameter' not in code.lower():
                findings.append(f'pattern (line {i}): {advice}\n        {code.strip()[:110]}')
    return findings


def duplicates(repo, files):
    tool = Path(repo, 'tools', 'dup_report.py')
    if not tool.exists():
        return []
    out = subprocess.run([sys.executable, str(tool), 'src,include', '400'], cwd=repo, capture_output=True,
                         text=True).stdout.splitlines()[1:]
    names = {Path(f).name for f in files}
    return [l.strip() for l in out if any(n in l for n in names)]


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    repo = sys.argv[1]
    files = sys.argv[2:] or changed_files(repo)
    if not files:
        print('no changed source in src/ or include/ (pass the files explicitly)')
        return
    total = 0
    for rel in files:
        f = check(repo, rel)
        total += len(f)
        print(f'\n== {rel}: {len(f)} finding(s)')
        for x in f:
            print(f'  - {x}')
    dups = duplicates(repo, files)
    print(f'\n== repeated blocks involving these files: {len(dups)}')
    for d in dups:
        print(f'  - {d[:160]}')
    print(f'\n{total} finding(s), {len(dups)} repeated block(s): justify or fix each one in the report.')


if __name__ == '__main__':
    main()
