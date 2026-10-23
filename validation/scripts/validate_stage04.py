#!/usr/bin/env python3
"""Reproduce stage 4 builds and regressions; require real GPU execution by default."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=ROOT / 'build/stage04-validation')
parser.add_argument('--cpu-only', action='store_true', help='Validate CPU integration without closing the GPU gate')
args = parser.parse_args()
OUT = args.output.resolve()
OUT.mkdir(parents=True, exist_ok=True)
summary = {'captured_at_utc': datetime.now(timezone.utc).isoformat(), 'profiles': {},
           'cpu_integration_validated': False, 'target_gpu_runtime_validated': False,
           'stage4_complete': False, 'stage5_started': False, 'amd_validated': False}

def run(name, cmd, overrides=None):
    try:
        p = subprocess.run([str(x) for x in cmd], cwd=ROOT, capture_output=True, text=True,
                           env={**os.environ, **(overrides or {})}, timeout=300)
        result = {'exit_code': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}
    except subprocess.TimeoutExpired:
        result = {'exit_code': None, 'stdout': '', 'stderr': 'Timed out after 300 seconds'}
    result['command'] = [str(x) for x in cmd]
    (OUT / (name + '.json')).write_text(json.dumps(result, indent=2) + '\n')
    (OUT / (name + '.log')).write_text(result['stdout'] + result['stderr'])
    return result

profiles = ['native-cpu-release', 'fypp-cpu-release', 'fypp-cpu-debug']
if not args.cpu_only:
    profiles.append('fypp-nvidia-gcc12')
for preset in profiles:
    item = summary['profiles'][preset] = {'status': 'configuration_error'}
    if run(preset + '-configure', ['cmake', '--preset', preset])['exit_code'] != 0:
        print(preset + ': configuration_error', flush=True)
        continue
    built = run(preset + '-build', ['cmake', '--build', '--preset', preset])
    item['build_passed'] = built['exit_code'] == 0
    if not item['build_passed']:
        item['status'] = 'build_error'
        print(preset + ': build_error', flush=True)
        continue
    junit = OUT / (preset + '-ctest.xml')
    junit.unlink(missing_ok=True)
    tests = run(preset + '-ctest', ['ctest', '--preset', preset, '-V', '--output-junit', junit])
    item['tests_passed'] = tests['exit_code'] == 0
    if junit.exists():
        cases = ET.parse(junit).getroot().findall('.//testcase')
        item['test_count'] = len(cases)
        item['failed_tests'] = [case.attrib.get('name') for case in cases
                                if any(case.find(tag) is not None for tag in ('failure', 'error', 'skipped'))]
    expected = 51 if preset == 'native-cpu-release' else 53
    item['tests_passed'] = bool(item['tests_passed'] and item.get('test_count') == expected
                                    and not item.get('failed_tests'))
    item['status'] = 'passed' if item['tests_passed'] else 'test_error'
    comparison = ROOT / 'build' / preset / 'tests/generated/reference-comparison/result.json'
    if preset != 'native-cpu-release':
        if comparison.exists():
            item['cpu_reference_comparison'] = json.loads(comparison.read_text())
            (OUT / (preset + '-comparison.json')).write_text(comparison.read_text())
        item['tests_passed'] = bool(item['tests_passed'] and
            item.get('cpu_reference_comparison', {}).get('status') == 'passed')
    if not item['tests_passed']:
        item['status'] = 'test_error'
    commands = json.loads((ROOT / 'build' / preset / 'compile_commands.json').read_text())
    host = [c['command'] for c in commands if '/muscle.dir/' in c['command']]
    item['cpu_has_no_gpu_flags'] = bool(host) and all('-foffload=nvptx-none' not in c for c in host)
    if preset == 'fypp-nvidia-gcc12':
        device = [c['command'] for c in commands if c['file'].endswith('muscle_kernels_target.f90')]
        required = ('-fopenmp', '-foffload=nvptx-none', '-foffload-options=nvptx-none=-march=sm_80')
        ninja = run(preset + '-commands', ['ninja', '-C', ROOT / 'build' / preset, '-t', 'commands'])
        links = [c for c in ninja['stdout'].splitlines() if ' -o tests/generated/muscle_test_generated_kernels ' in c]
        item['device_flags_verified'] = bool(device) and len(links) == 1 and all(
            all(flag in c for flag in required) for c in device + links)
        item['device_execution_verified'] = item['tests_passed'] and all(
            marker in tests['stdout'] for marker in
            ('calculation_host_points=0', 'PASS_GENERATED_KERNELS_DEVICE'))
        if not item['tests_passed'] and item.get('failed_tests') == ['generated_kernels'] and any(
                token in tests['stdout'] + tests['stderr'] for token in ('libgomp:', 'OFFLOAD_ERROR', 'cuInit')):
            item['status'] = 'offload_error'
    print(preset + ': ' + item['status'], flush=True)

summary['cpu_integration_validated'] = all(
    summary['profiles'][p].get('tests_passed') and summary['profiles'][p].get('cpu_has_no_gpu_flags')
    for p in profiles[:3])
gpu = summary['profiles'].get('fypp-nvidia-gcc12', {})
summary['target_gpu_build_validated'] = bool(gpu.get('build_passed') and gpu.get('device_flags_verified')
                                           and gpu.get('cpu_has_no_gpu_flags'))
summary['target_gpu_runtime_validated'] = bool(gpu.get('device_execution_verified'))
summary['stage4_complete'] = bool(summary['cpu_integration_validated'] and
    summary['target_gpu_build_validated'] and summary['target_gpu_runtime_validated'])
summary['status'] = ('passed' if summary['stage4_complete'] else
                     'cpu_only_passed' if args.cpu_only and summary['cpu_integration_validated'] else 'pending_or_failed')
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({k: v for k, v in summary.items() if k != 'profiles'}, indent=2), flush=True)
raise SystemExit(0 if summary['stage4_complete'] or (args.cpu_only and summary['cpu_integration_validated']) else 2)
