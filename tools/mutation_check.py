#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
"""Mutation check of the tests of one area (plan R of docs/fypp/ESTADO_OPTIMO.md).

A reference test may be deleted only if the unit tests that replace it detect the same deliberate
mutations of the code. This script applies, one at a time, mutations of the sources of an area to the
working tree, rebuilds the tests (dp only, Debug), runs the reference tests and the unit tests of the
area, restores the source and reports which group detected (killed) each mutation.

The mutations are sampled with a fixed seed from the executable lines of the sources, round-robin over the
files (so the sample does not depend on who writes the tests), plus the hand-written ones of the area:

- ``aor``: a binary ``+``/``-``/``*``/``/`` replaced by its opposite;
- ``lit``: a real literal changed (0 -> 1, x -> x + 1);
- ``idx``: an integer subscript changed by one (1 -> 2, n -> n - 1);
- ``ror``: a relational operator relaxed or tightened (``<`` <-> ``<=``, ``>`` <-> ``>=``).

The area is described by ``tests/mutations/<area>.toml``::

    sources = ["src/math/muscle_math_operations.fypp"]               # globs, relative to the repository
    reference = "^math-math_reference$"                             # ctest regex of the reference tests
    unit = "^math-"                                                 # ctest regex of the unit tests
    samples = 40                                                    # sampled mutations (default 40; "all")
    seed = 1                                                        # default 1
    [[manual]]                                                      # optional, any number
    file = "src/plasticity/hardening_laws/muscle_hard_swift.fypp"
    old = "res = prm%k*(prm%e0 + ep)**prm%n"                        # must occur exactly once
    new = "res = prm%k*(prm%e0 + ep)**(prm%n + 0.01_${k}$)"

The unit group runs ``ctest -R <unit> -E <reference>|_reference_mutation`` (the self-checks of the reference
comparison are not unit tests). With ``reference = ""`` (the reference of the
area was deleted) only the unit group runs: the report then shows what the unit tests detect. Criterion: every mutation killed by the
reference group is also killed by the unit group (report line ``CRITERION: OK``). A mutation that does not
compile is skipped; one that no group kills is listed (equivalent mutant or a gap of both).

Usage: tools/mutation_check.py AREA [--build DIR] [--jobs N] [--list] [--only ID,...] [--report FILE]
                               [--samples N|all] [--seed S] [--reference REGEX] [--unit REGEX]
  --build: dp-only Debug build directory (default build-mutation; configured if missing). Never run two
           checks on the same working tree at the same time: the mutation is written to the source.
  --list:  print the sampled mutations without building.
  --reference, --unit: override the groups of the area; e.g. to measure whether the legacy tests of an area
           detect anything the others do not: ``--reference '<legacy tests>' --unit '<all tests of the area>'``
           (the unit group excludes the reference group).
"""
import argparse
import glob
import os
import random
import re
import signal
import subprocess
import sys
import time
import tomllib
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

DECL = re.compile(r"^\s*(real|integer|logical|character|type\s*\(|class\s*\(|use\b|implicit|public|private|"
                  r"interface|end\b|contains|module|pure|elemental|function|subroutine|procedure|"
                  r"!|#|@|\$|&\s*!|recursive|impure|import|include|generic|program|return|cycle|exit)",
                  re.IGNORECASE)
REAL_LIT = re.compile(r"(?<![\w.])(\d+\.\d*)(?:([eEdD])([+-]?\d+))?(?=_)")
SUBSCRIPT = re.compile(r"(?<=[\w%)]\()(\d+)(?=\s*[,)])|(?<=,)(\s*)(\d+)(?=\s*\))")
ROR = [("<=", "<"), (">=", ">"), ("<", "<="), (">", ">=")]


@dataclass
class Mutation:
    id: str
    file: str
    line: int          # 1-based; 0 for a manual mutation located by `old`
    old: str
    new: str
    kind: str


def code_part(line):
    """The code of a line without its trailing comment, and the comment."""
    in_str = None
    for i, c in enumerate(line):
        if in_str:
            if c == in_str:
                in_str = None
        elif c in "'\"":
            in_str = c
        elif c == "!":
            return line[:i], line[i:]
    return line, ""


def masked(code):
    """Positions inside ${...}$ or strings, which no mutation may touch."""
    mask = [False]*len(code)
    for m in re.finditer(r"\$\{.*?\}\$|'[^']*'|\"[^\"]*\"", code):
        for i in range(m.start(), m.end()):
            mask[i] = True
    return mask


def is_executable(line):
    code, _ = code_part(line)
    s = code.strip()
    if not s or DECL.match(s):
        return False
    # an assignment, a condition or a continuation of an expression
    return "=" in s or s.startswith(("if", "else if", "do ")) or s.endswith("&") or s.startswith("&")


def candidates(path, lines):
    """Every mutation of every executable line of a source."""
    out = []
    rel = os.path.relpath(path, REPO)
    directive = False      # inside a continued FYPP directive (#:set x = ... &)
    for no, line in enumerate(lines, 1):
        if directive or line.lstrip().startswith("#:"):
            directive = line.rstrip().endswith("&")
            continue
        if not is_executable(line):
            continue
        code, comment = code_part(line)
        mask = masked(code)
        start = code.find("=") + 1 if re.match(r"^\s*[\w%()., ]+=(?!=)", code) else 0

        def add(kind, i, j, rep):
            if any(mask[i:j]) or i < start:
                return
            new = code[:i] + rep + code[j:] + comment
            out.append(Mutation(f"{kind}:{rel}:{no}:{i}", rel, no, line, new, kind))

        for m in re.finditer(r"[+\-*/]", code):
            i, c = m.start(), m.group()
            prev = code[:i].rstrip()
            nxt = code[i + 1:].lstrip()
            if c == "*" and (code[i - 1:i] == "*" or code[i + 1:i + 2] == "*"):
                continue
            if c == "/" and (code[i - 1:i] in "(/" or code[i + 1:i + 2] in "/=)"):
                continue
            if c in "+-":
                if not prev or prev[-1] in "=(,*/+-<>:" or re.search(r"\d[eEdD]$", prev) \
                        or re.search(r"\.(and|or|not|eqv|neqv|lt|le|gt|ge|eq|ne)\.$", prev, re.I):
                    continue        # unary sign or exponent
            if not nxt or nxt[0] == "&":
                continue
            add("aor", i, i + 1, {"+": "-", "-": "+", "*": "/", "/": "*"}[c])
        for m in REAL_LIT.finditer(code):
            v = float(m.group(1)) * (10.0 ** int(m.group(3)) if m.group(3) else 1.0)
            rep = "1.0" if v == 0.0 else repr(v + 1.0)
            add("lit", m.start(), m.end(), rep)
        for m in SUBSCRIPT.finditer(code):
            g = 1 if m.group(1) else 3
            v = int(m.group(g))
            add("idx", m.start(g), m.end(g), str(v + 1 if v == 1 else v - 1) if v > 0 else "1")
        for m in re.finditer(r"<=|>=|(?<![=/<>])<(?!=)|(?<![=/<>-])>(?!=)", code):
            rep = dict(ROR)[m.group()]
            add("ror", m.start(), m.end(), rep)
    return out


def sample(area, cfg):
    files = sorted({f for g in cfg["sources"] for f in glob.glob(str(REPO/g))})
    if not files:
        sys.exit(f"no sources match {cfg['sources']}")
    rng = random.Random(cfg.get("seed", 1))
    pools = []
    for f in files:
        c = candidates(f, Path(f).read_text().splitlines(keepends=True))
        rng.shuffle(c)
        if c:
            pools.append(c)
    picked, n = [], cfg.get("samples", 40)
    if n == "all":
        n = sum(len(p) for p in pools)
    while len(picked) < n and any(pools):
        for p in pools:
            if p and len(picked) < n:
                picked.append(p.pop())
    for i, m in enumerate(cfg.get("manual", []), 1):
        picked.append(Mutation(m.get("id", f"manual{i}"), m["file"], 0, m["old"], m["new"], "manual"))
    return picked


def run(cmd, timeout=None, **kw):
    return subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=timeout, **kw)


def build(bdir, jobs):
    r = run(["cmake", "--build", str(bdir), "--target", "build_all_tests", "-j", str(jobs)])
    return r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]


def ctest(bdir, regex, exclude, jobs):
    cmd = ["ctest", "--test-dir", str(bdir), "-R", regex, "-j", str(jobs), "--timeout", "120"]
    if exclude:
        cmd += ["-E", exclude]
    r = run(cmd)
    if "No tests were found" in r.stdout + r.stderr:
        sys.exit(f"no tests match {regex}")
    failed = re.findall(r"^\s*\d+ - (\S+) \(", r.stdout, re.M)
    return r.returncode != 0, failed


class Applied:
    """Writes a mutation to its source and always restores it."""

    def __init__(self, m):
        self.m = m
        self.path = REPO/m.file
        self.orig = self.path.read_text()

    def __enter__(self):
        if self.m.line:
            lines = self.orig.splitlines(keepends=True)
            if lines[self.m.line - 1] != self.m.old:
                raise RuntimeError(f"{self.m.id}: source changed since the sample")
            lines[self.m.line - 1] = self.m.new
            text = "".join(lines)
        else:
            if self.orig.count(self.m.old) != 1:
                raise RuntimeError(f"{self.m.id}: `old` must occur exactly once in {self.m.file}")
            text = self.orig.replace(self.m.old, self.m.new)
        self.path.write_text(text)
        return self

    def __exit__(self, *exc):
        self.path.write_text(self.orig)
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("area")
    ap.add_argument("--build", default="build-mutation")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only")
    ap.add_argument("--report")
    ap.add_argument("--samples", help='override the samples of the area (a number or "all")')
    ap.add_argument("--seed", type=int, help="override the seed of the area")
    ap.add_argument("--reference", help="override the ctest regex of the reference group")
    ap.add_argument("--unit", help="override the ctest regex of the unit group")
    a = ap.parse_args()
    cfg = tomllib.loads((REPO/"tests"/"mutations"/f"{a.area}.toml").read_text())
    for group in ("reference", "unit"):
        if getattr(a, group) is not None:
            cfg[group] = getattr(a, group)
    if a.samples:
        cfg["samples"] = a.samples if a.samples == "all" else int(a.samples)
    if a.seed is not None:
        cfg["seed"] = a.seed
    muts = sample(a.area, cfg)
    if a.only:
        ids = set(a.only.split(","))
        muts = [m for m in muts if m.id in ids]
    if a.list:
        for m in muts:
            print(f"{m.id}\n  - {m.old.strip()}\n  + {m.new.strip()}")
        return 0

    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))
    bdir = (REPO/a.build).resolve()
    if not (bdir/"CMakeCache.txt").exists():
        r = run(["cmake", "-S", ".", "-B", str(bdir), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Debug",
                 "-DMUSCLE_REAL_KINDS=dp"])
        if r.returncode:
            sys.exit(r.stdout + r.stderr)
    ref, unit = cfg["reference"], cfg["unit"]
    unit_excl = f"{ref}|_reference_mutation" if ref else "_reference_mutation"
    ok, log = build(bdir, a.jobs)
    if not ok:
        sys.exit("baseline build failed\n" + log)
    for name, rx, ex in (("reference", ref, None), ("unit", unit, unit_excl)):
        if not rx:
            continue
        bad, failed = ctest(bdir, rx, ex, a.jobs)
        if bad:
            sys.exit(f"baseline {name} tests fail: {failed}")

    rows, t0 = [], time.time()
    for i, m in enumerate(muts, 1):
        with Applied(m):
            ok, _ = build(bdir, a.jobs)
            if ok:
                kr, fr = ctest(bdir, ref, None, a.jobs) if ref else (False, [])
                ku, fu = ctest(bdir, unit, unit_excl, a.jobs)
        if not ok:
            rows.append((m, None, None, [], []))
            status = "does not compile"
        else:
            rows.append((m, kr, ku, fr, fu))
            status = f"reference {'KILLED' if kr else 'survives'}, unit {'KILLED' if ku else 'survives'}"
        print(f"[{i}/{len(muts)} {time.time() - t0:.0f}s] {m.id}: {status}", flush=True)
    build(bdir, a.jobs)    # leave the build matching the restored sources

    lost = [r for r in rows if r[1] and not r[2]]
    lines = [f"# Mutation check: {a.area}", "",
             f"reference `{ref}`, unit `{unit}` (excluding `{unit_excl}`), {len(muts)} mutations, seed "
             f"{cfg.get('seed', 1)}", "",
             "| mutation | change | reference | unit |", "|---|---|---|---|"]
    for m, kr, ku, fr, fu in rows:
        ch = f"`{m.old.strip()}` → `{m.new.strip()}`".replace("|", "\\|")
        f = lambda k: "—" if k is None else ("killed" if k else "survives")
        lines.append(f"| {m.id} | {ch} | {f(kr)} | {f(ku)} |")
    n_ref = sum(1 for r in rows if r[1])
    n_both = sum(1 for r in rows if r[1] and r[2])
    lines += ["", f"Killed by the reference: {n_ref}; also by the unit tests: {n_both}; "
              f"killed only by the unit tests: {sum(1 for r in rows if r[2] and not r[1])}; "
              f"killed by none: {sum(1 for r in rows if r[1] is False and r[2] is False)}; "
              f"not compiled: {sum(1 for r in rows if r[1] is None)}.",
              "CRITERION: OK" if not lost else
              "CRITERION: FAILED, killed only by the reference: " + ", ".join(r[0].id for r in lost)]
    text = "\n".join(lines) + "\n"
    if a.report:
        Path(a.report).write_text(text)
    print(text)
    return 1 if lost else 0


if __name__ == "__main__":
    sys.exit(main())
