#!/usr/bin/env python3
"""Optional evidence of the environment passed to the NVIDIA compiler driver."""
import json
import os
import shutil
import sys

record = {"path": os.environ.get("PATH"), "ptxas_found": shutil.which("ptxas"),
          "compiler": "/usr/bin/gfortran-12", "arguments": sys.argv[2:]}
descriptor = os.open(sys.argv[1], os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o600)
try:
    os.write(descriptor, (json.dumps(record) + "\n").encode())
finally:
    os.close(descriptor)
if record["path"] != "/usr/bin:/bin" or record["ptxas_found"] is not None:
    raise SystemExit("NVIDIA driver environment is not isolated from ptxas.")
