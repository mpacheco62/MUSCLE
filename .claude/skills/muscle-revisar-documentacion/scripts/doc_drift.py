#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Drift between the documentation and the code of MUSCLE.

    doc_drift.py <repo> [--only undocumented|names|paths|decisions|models]

Checks (each prints its findings; exit status 1 if any):
- undocumented: public modules, types and procedures of src/ (templates) without their FORD `!!` line right after
  the declaration (module, type, function or subroutine statement and its continuation lines);
- names: identifiers in backticks in pages/*.md, README.md, my_project.md and .claude/skills/**/*.md that look like
  code (snake_case with a known prefix, or a MUSCLE module) and appear nowhere in src/, include/, tools/, tests/,
  benchmarks/, cmake/ or CMakeLists.txt (renamed or retired names);
- paths: paths in backticks (with a / and an extension, or a known top directory) that do not exist;
- links: links between FORD pages, `](page/<name>.html)`, whose pages/<name>.md does not exist;
- decisions: references D-0NN (and D-0NN.x) to decisions that docs/fypp/DECISIONES.md does not have;
- models: model modules of each family (muscle_<prefix>_<model>) that README.md, my_project.md and pages/ never
  mention (by module name, model name or its prefix_model).
"""
import argparse
import pathlib
import re
import sys

PREFIXES = ('muscle_', 'hard_', 'vp_', 'os_', 'yield_', 'elas_', 'visco_', 'kin_', 'mat_', 'ten_', 'iden_',
            'fp_', 'eigen', 'spectral_', 'polar_', 'log_', 'lu_', 'cond_', 'norm_', 'STATUS_', 'TOL_', 'KIN_',
            'EPS_', 'NEWTON_', 'JACOBI_', 'FD_', 'QUIET_')
CODE_DIRS = ('src', 'include', 'tools', 'tests', 'benchmarks', 'cmake')
DOC_GLOBS = ('pages/*.md', 'README.md', 'my_project.md', '.claude/skills/**/*.md')
FAMILY_DIRS = {'hard': 'src/plasticity/hardening_laws', 'vp': 'src/plasticity/viscoplastic_laws',
               'os': 'src/plasticity/viscoplastic_laws', 'yield': 'src/plasticity/yield_criteria',
               'elas': 'src/elasticity', 'visco': 'src/viscoelasticity', 'kin': 'src/kinematics'}
NOT_MODELS = {'laws', 'types', 'anisotropic', 'finite', 'frozen'}


def read(p):
    return p.read_text(errors='replace')


def code_corpus(repo):
    words = set()
    for d in CODE_DIRS:
        for p in (repo / d).rglob('*'):
            if p.is_file() and p.suffix in ('.fypp', '.f90', '.F90', '.py', '.sh', '.cmake', '.txt', '.toml'):
                text = re.sub(r'\$\{[^}]*\}\$', '', read(p))     # name_${k}$ -> name_
                words.update(w.lower().rstrip('_') for w in re.findall(r'[A-Za-z_]\w*', text))
    cm = repo / 'CMakeLists.txt'
    if cm.exists():
        words.update(w.lower().rstrip('_') for w in re.findall(r'[A-Za-z_]\w*', read(cm)))
    return words


def doc_files(repo):
    out = []
    for g in DOC_GLOBS:
        out += sorted(repo.glob(g))
    return [p for p in out if p.is_file()]


def backticks(text):
    return re.findall(r'`([^`\n]+)`', text)


def check_undocumented(repo):
    findings = []
    decl = re.compile(r'^\s*(?:(?:pure|elemental|impure|recursive)\s+)*(module|function|subroutine|type)\b'
                      r'(?!\s*(?:procedure|\(|is\b|,\s*extends))\s*(?:,\s*[\w\s(),]+::\s*)?([\w${}]+)', re.I)
    for p in sorted((repo / 'src').rglob('*.fypp')):
        lines = read(p).split('\n')
        public = set()
        for l in lines:
            m = re.match(r'\s*public\s*::\s*(.+)', l, re.I)
            if m:
                public.update(n.strip().lower() for n in re.split(r',', m.group(1)) if n.strip())
        for i, l in enumerate(lines):
            m = decl.match(l)
            if not m or re.match(r'\s*end\b', l, re.I):
                continue
            kind, name = m.group(1).lower(), m.group(2)
            if kind == 'module' and name.lower() == 'procedure':
                continue
            base = re.sub(r'_\$\{\w+\}\$$', '', name).lower()
            if kind != 'module' and not (name.lower() in public or base in public
                                         or any(n.startswith(base) for n in public if '$' in n)):
                continue
            j = i
            while lines[j].rstrip().endswith('&'):
                j += 1
            nxt = next((x for x in lines[j + 1:j + 4] if x.strip() and not x.strip().startswith('#')), '')
            if not nxt.strip().startswith('!!'):
                findings.append(f'{p.relative_to(repo)}:{i + 1}: {kind} {name} without a `!!` line')
    return findings


#: Names that the documentation uses as examples of a name a user would write (not in the code)
EXAMPLES = {'muscle_hard_mylaw', 'hard_mylaw'}
#: Retired names that pages/migration.md documents on purpose (the former API)
RETIRED = {'muscle_kinematics'}


def known(w, corpus):
    """w, or w without the kind/dimension suffixes of a generated specific, is a name of the code."""
    w = w.lower()
    cands = {w, re.sub(r'(_[23]d)?(_[sd]p)$', '', w), re.sub(r'_[sd]p$', '', w)}
    return any(c in corpus for c in cands)


def check_names(repo, corpus):
    findings = []
    for p in doc_files(repo):
        for i, l in enumerate(read(p).split('\n'), 1):
            for t in backticks(l):
                # a name followed by a placeholder (`muscle_mat_<name>`, `name_${k}$`) is a pattern, not a name
                for m in re.finditer(r'[A-Za-z_]\w*', t):
                    w, rest = m.group(0), t[m.end():]
                    if rest.startswith(('<', '$', '{')) or w.endswith('_'):
                        continue
                    if (len(w) < 5 or not w.startswith(PREFIXES) or '_' not in w.strip('_')
                            or w.lower() in EXAMPLES or (p.name == 'migration.md' and w.lower() in RETIRED)):
                        continue
                    if not known(w, corpus):
                        findings.append(f'{p.relative_to(repo)}:{i}: `{w}` not found in the code')
    return sorted(set(findings))


def check_paths(repo):
    findings = []
    top = ('src/', 'tests/', 'tools/', 'include/', 'docs/', 'pages/', 'cmake/', 'benchmarks/', '.claude/', 'Theory/')
    for p in doc_files(repo):
        for i, l in enumerate(read(p).split('\n'), 1):
            for t in backticks(l):
                t = t.strip()
                if ' ' in t or '<' in t or '*' in t or '{' in t or '$' in t:
                    continue
                if not (t.startswith(top) or (('/' in t) and re.search(r'\.\w+$', t))):
                    continue
                if any(e in t for e in EXAMPLES):
                    continue
                # relative to the repository, to src/, or to any directory above the document (a skill root)
                cands = [repo / t, repo / 'src' / t] + [d / t for d in p.parents if repo in d.parents or d == repo]
                if not any(c.exists() for c in cands):
                    findings.append(f'{p.relative_to(repo)}:{i}: path `{t}` does not exist')
    return sorted(set(findings))


def check_links(repo):
    pages = {p.stem for p in (repo / 'pages').glob('*.md')}
    findings = []
    for p in doc_files(repo):
        for i, l in enumerate(read(p).split('\n'), 1):
            for name in re.findall(r'\]\((?:\.\./)?page/([\w-]+)\.html', l):
                if name not in pages:
                    findings.append(f'{p.relative_to(repo)}:{i}: link to page/{name}.html, but pages/{name}.md does not exist')
    return findings


def check_decisions(repo):
    dec = read(repo / 'docs/fypp/DECISIONES.md')
    have = set(re.findall(r'^## (D-\d{3})', dec, re.M))
    findings = []
    files = doc_files(repo) + sorted((repo / 'docs/fypp').glob('*.md')) + sorted((repo / 'src').rglob('*.fypp'))
    for p in files:
        for i, l in enumerate(read(p).split('\n'), 1):
            for d in re.findall(r'\bD-(\d{3})\b', l):
                if f'D-{d}' not in have:
                    findings.append(f'{p.relative_to(repo)}:{i}: D-{d} is not in DECISIONES.md')
    return sorted(set(findings))


def check_models(repo):
    docs = '\n'.join(read(p) for p in doc_files(repo) if '.claude' not in str(p)).lower()
    lines = docs.split('\n')
    findings = []
    for fam, d in FAMILY_DIRS.items():
        for p in sorted((repo / d).glob(f'muscle_{fam}_*.fypp')):
            model = p.stem[len(f'muscle_{fam}_'):]
            if model in NOT_MODELS:
                continue
            words = model.split('_')
            if not (any(k in docs for k in (p.stem, f'{fam}_{model}', model))
                    or any(all(w in l for w in words) for l in lines)):
                findings.append(f'{p.relative_to(repo)}: model `{model}` of family `{fam}` is not mentioned in '
                                'README.md, my_project.md or pages/')
    return findings


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('repo')
    ap.add_argument('--only', choices=('undocumented', 'names', 'paths', 'links', 'decisions', 'models'))
    a = ap.parse_args()
    repo = pathlib.Path(a.repo).resolve()
    checks = {'undocumented': lambda: check_undocumented(repo), 'names': lambda: check_names(repo, code_corpus(repo)),
              'paths': lambda: check_paths(repo), 'links': lambda: check_links(repo),
              'decisions': lambda: check_decisions(repo),
              'models': lambda: check_models(repo)}
    total = 0
    for name, fn in checks.items():
        if a.only and name != a.only:
            continue
        f = fn()
        total += len(f)
        print(f'== {name}: {len(f)} finding(s)')
        for x in f:
            print('  ' + x)
    sys.exit(1 if total else 0)


if __name__ == '__main__':
    main()
