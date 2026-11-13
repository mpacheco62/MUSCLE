# SPDX-License-Identifier: GPL-3.0-or-later
"""List the signatures of the generic operators/assignments of Fortran tensor modules.

Each line is ``<generic> (<type arg1>, <type arg2>) -> <result type>``, independent of the
specific procedure names, so the historical double precision set can be compared with the
generated one (D-015): kind suffixes ``_dp`` and kind parameters ``real64``/``dp`` are
normalized to the historical names.

Usage:
    ops_signatures.py FILE...                 print the signatures (sorted)
    ops_signatures.py --check REF FILE...     fail if a signature of REF is missing
"""
import re
import sys

_GENERIC = re.compile(r'^\s*interface\s+(operator\s*\(\s*\S+?\s*\)|assignment\s*\(\s*=\s*\))', re.I)
_MODPROC = re.compile(r'^\s*module\s+procedure\s+(.+)$', re.I)
_PROC = re.compile(r'^\s*(?:(?:pure|elemental|impure|recursive)\s+)*(function|subroutine)\s+(\w+)\s*'
                   r'\(([^)]*)\)(?:\s*result\s*\(\s*(\w+)\s*\))?', re.I)
_END = re.compile(r'^\s*end\s+(function|subroutine)\b', re.I)
_DECL = re.compile(r'^\s*(type\s*\(\s*\w+\s*\)|real\s*\(\s*\w+\s*\)|logical|integer)\s*(?:,[^:]*)?::\s*(.+)$',
                   re.I)


def _norm(t):
    t = re.sub(r'\s+', '', t).lower()
    t = re.sub(r'_dp\)$', ')', t)
    return t.replace('real(real64)', 'real(dp)')


def _strip(line):
    return line.split('!', 1)[0].rstrip()


def parse(paths):
    generics = {}   # procedure name -> generic
    procs = {}      # procedure name -> (args, {var: type})
    for path in paths:
        with open(path, encoding='utf-8') as f:
            lines = [_strip(l) for l in f]
        generic = None
        cur = None
        for line in lines:
            m = _GENERIC.match(line)
            if m:
                generic = re.sub(r'\s+', '', m.group(1)).lower()
                continue
            if re.match(r'^\s*end\s+interface', line, re.I):
                generic = None
                continue
            if generic:
                m = _MODPROC.match(line)
                if m:
                    for name in m.group(1).split(','):
                        generics[name.strip().lower()] = generic
                continue
            m = _PROC.match(line)
            if m:
                args = [a.strip().lower() for a in m.group(3).split(',') if a.strip()]
                res = (m.group(4) or m.group(2)).lower() if m.group(1).lower() == 'function' else None
                cur = (m.group(2).lower(), args, res, {})
                continue
            if cur:
                if _END.match(line):
                    procs[cur[0]] = cur[1:]
                    cur = None
                    continue
                m = _DECL.match(line)
                if m:
                    for var in m.group(2).split(','):
                        cur[3][var.split('(')[0].split('=')[0].strip().lower()] = _norm(m.group(1))
    sigs = set()
    for name, generic in generics.items():
        if name not in procs:
            raise SystemExit(f'{name}: listed in {generic} but not defined in the given files')
        args, res, types = procs[name]
        sig = f"{generic} ({', '.join(types[a] for a in args)})"
        if res:
            sig += f' -> {types[res]}'
        sigs.add(sig)
    return sigs


def main(argv):
    if argv and argv[0] == '--check':
        with open(argv[1], encoding='utf-8') as f:
            ref = {l.strip() for l in f if l.strip()}
        missing = sorted(ref - parse(argv[2:]))
        for sig in missing:
            print(f'missing: {sig}')
        return 1 if missing else 0
    print('\n'.join(sorted(parse(argv))))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
