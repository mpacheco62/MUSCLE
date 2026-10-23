#!/usr/bin/env python3
"""Compare CPU API, concrete CPU solver, and the stored 768-state baseline."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reference', required=True)
parser.add_argument('--concrete', required=True)
parser.add_argument('--baseline', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
result = args.output / 'result.json'
result.unlink(missing_ok=True)
report = {'status': 'execution_error', 'rows': 768, 'rtol': 1e-11, 'atol': 1e-10}

def finish(status, code, **details):
    report.update(status=status, **details)
    result.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    raise SystemExit(code)

def execute(label, executable):
    path = args.output / (label + '.csv')
    path.unlink(missing_ok=True)
    try:
        p = subprocess.run([executable, str(path)], capture_output=True, text=True,
                           env={**os.environ, 'OMP_TARGET_OFFLOAD': 'DISABLED'}, timeout=45)
    except subprocess.TimeoutExpired:
        finish('execution_error', 2, failed_executable=label, reason='Timeout after 45 seconds')
    (args.output / (label + '.log')).write_text(p.stdout + p.stderr)
    if p.returncode:
        finish('execution_error', 2, failed_executable=label, exit_code=p.returncode)
    return path

def read(path):
    with path.open() as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        rows = list(reader)
    if len(rows) != 768 or len(fields or []) != 76:
        finish('output_error', 2, reason='Expected 768 states with 76 columns', path=str(path))
    return fields, rows

ref = execute('cpu', args.reference)
concrete = execute('concrete', args.concrete)
fields, baseline = read(args.baseline)
identity = fields[:6]
maximum = {key: 0.0 for key in fields[6:]}
comparisons = {}
for label, path in [('cpu_vs_baseline', ref), ('concrete_vs_baseline', concrete)]:
    names, rows = read(path)
    if names != fields:
        finish('output_error', 2, reason='CSV headers differ', comparison=label)
    failures = []
    exact = True
    for index, (old, new) in enumerate(zip(baseline, rows), 1):
        for key in identity:
            if old[key] != new[key]:
                failures.append({'row': index, 'field': key, 'expected': old[key], 'actual': new[key]})
        for key in fields[6:]:
            try:
                a, b = float(old[key]), float(new[key])
            except (ValueError, TypeError):
                finish('output_error', 2, reason='Malformed numeric field', row=index, field=key)
            d = abs(a - b)
            exact = exact and a == b
            maximum[key] = max(maximum[key], d)
            if not (math.isfinite(a) and math.isfinite(b)) or d > 1e-10 + 1e-11 * abs(a):
                failures.append({'row': index, 'field': key, 'expected': a, 'actual': b})
    comparisons[label] = {'failure_count': len(failures), 'exact_numeric_agreement': exact}
    if failures:
        finish('numerical_error', 3, comparison=label, failure_count=len(failures), first_failures=failures[:10])
# Also compare CPU to concrete directly, rather than only to the stored baseline.
_, cpu_rows = read(ref)
_, concrete_rows = read(concrete)
failures = []
for i, (a, b) in enumerate(zip(cpu_rows, concrete_rows), 1):
    for key in fields:
        if key in identity:
            valid = a[key] == b[key]
        else:
            x, y = float(a[key]), float(b[key])
            valid = math.isfinite(x) and math.isfinite(y) and abs(x-y) <= 1e-10 + 1e-11*abs(x)
        if not valid:
            failures.append({'row': i, 'field': key})
if failures:
    finish('numerical_error', 3, comparison='cpu_vs_concrete', failure_count=len(failures), first_failures=failures[:10])
comparisons['cpu_vs_concrete'] = {'failure_count': 0, 'exact_numeric_agreement': all(
    all(a[k] == b[k] if k in identity else float(a[k]) == float(b[k]) for k in fields)
    for a, b in zip(cpu_rows, concrete_rows))}
finish('passed', 0, comparisons=comparisons, max_absolute_difference_by_field=maximum,
       max_absolute_difference=max(maximum.values()), gpu_solver_validated=False)
