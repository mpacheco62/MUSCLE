#!/usr/bin/env python3
"""Compare unchanged feature kernels; report offload and numerical errors separately."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--reference", required=True)
parser.add_argument("--selected", required=True)
parser.add_argument("--mode", choices=("host", "gpu"), required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
(args.output / "result.json").unlink(missing_ok=True)
report = {"mode": args.mode, "gpu_validated": False}


def finish(status, code, **details):
    report.update(status=status, **details)
    (args.output / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(code)


def execute(label, executable, policy):
    csv_path = args.output / (label + ".csv")
    csv_path.unlink(missing_ok=True)
    command = [executable, str(csv_path)]
    try:
        p = subprocess.run(command, capture_output=True, text=True, timeout=45,
                           env={**os.environ, "OMP_TARGET_OFFLOAD": policy})
        result = {"command": command, "policy": policy, "exit_code": p.returncode,
                  "stdout": p.stdout, "stderr": p.stderr}
    except subprocess.TimeoutExpired:
        result = {"command": command, "policy": policy, "exit_code": None,
                  "stdout": "", "stderr": "Execution timed out after 45 seconds."}
    report[label] = result
    (args.output / (label + ".log")).write_text(result["stdout"] + result["stderr"])
    return result, csv_path


reference, ref_path = execute("reference", args.reference, "DISABLED")
if reference["exit_code"] != 0:
    finish("reference_execution_error", 4)
selected, sel_path = execute("selected", args.selected,
                             "MANDATORY" if args.mode == "gpu" else "DISABLED")
if selected["exit_code"] != 0:
    diagnostic = selected["stdout"] + selected["stderr"]
    if "NUMERICAL_ERROR" in diagnostic:
        finish("numerical_error", 3)
    if any(token in diagnostic for token in ("OFFLOAD_ERROR", "libgomp:", "cuInit", "CUDA_ERROR")):
        finish("offload_error", 2)
    finish("selected_execution_error", 4)
if args.mode == "gpu":
    stdout = selected["stdout"]
    markers = ("entry_initial_device=F", "calculation_host_points=0",
               "reduction_host_visits=0", "PASS_SNAPSHOT_DEVICE")
    if not all(marker in stdout for marker in markers):
        finish("offload_error", 2, reason="Missing proof of device execution in every calculation region.")
try:
    with ref_path.open() as stream:
        ref_rows = list(csv.DictReader(stream))
    with sel_path.open() as stream:
        sel_rows = list(csv.DictReader(stream))
except OSError as error:
    finish("output_error", 4, reason=str(error))
if len(ref_rows) != 257 or len(sel_rows) != 257:
    finish("numerical_error", 3, reason="Unexpected sample count.")

# Field-specific tolerances for these real64 compiler probes, not tolerances
# for the nonlinear constitutive solver. Never widen them after a failed run.
tolerances = {**{f"stress_{i}": {"atol": 1e-10, "rtol": 1e-11} for i in range(1, 7)},
              "vm": {"atol": 1e-10, "rtol": 1e-11},
              "hard": {"atol": 1e-12, "rtol": 1e-11},
              "inverse_residual": {"atol": 1e-12, "rtol": 0.0}}
expected_fields = {"point", "status", *tolerances}
if any(set(row) != expected_fields for row in ref_rows + sel_rows):
    finish("output_error", 4, reason="Unexpected CSV columns.")
maximums = {key: 0.0 for key in tolerances}
differences = []
for index, (ref, selected_row) in enumerate(zip(ref_rows, sel_rows), 1):
    if ref["point"] != str(index) or selected_row["point"] != str(index) or ref["status"] != selected_row["status"]:
        differences.append({"point": index, "field": "point/status", "reference": ref,
                            "selected": selected_row})
    for field, tol in tolerances.items():
        try:
            a, b = float(ref[field]), float(selected_row[field])
        except (TypeError, ValueError):
            finish("output_error", 4, reason=f"Malformed value at point {index}, field {field}.")
        difference = abs(a-b)
        maximums[field] = max(maximums[field], difference)
        if not (math.isfinite(a) and math.isfinite(b)) or difference > tol["atol"] + tol["rtol"]*abs(a):
            differences.append({"point": index, "field": field, "reference": a, "selected": b})
try:
    ref_checksum = float(re.search(r"^checksum=\s*(\S+)", reference["stdout"], re.M)[1])
    sel_checksum = float(re.search(r"^checksum=\s*(\S+)", selected["stdout"], re.M)[1])
except (TypeError, ValueError):
    finish("output_error", 4, reason="Missing checksum.")
if not math.isfinite(sel_checksum) or abs(ref_checksum-sel_checksum) > 1e-8 + 1e-11*abs(ref_checksum):
    differences.append({"field": "checksum", "reference": ref_checksum, "selected": sel_checksum})
report.update(samples=len(ref_rows), tolerances=tolerances, max_absolute_difference=maximums,
              checksum_absolute_difference=abs(ref_checksum-sel_checksum),
              difference_count=len(differences), first_differences=differences[:20])
if differences:
    finish("numerical_error", 3)
report["gpu_validated"] = args.mode == "gpu"
finish("passed", 0)
