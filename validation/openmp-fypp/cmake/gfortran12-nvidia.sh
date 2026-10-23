#!/bin/sh
# Only the NVIDIA toolchain uses this driver. The restricted PATH is inherited
# by mkoffload/nvptx-as during both compilation and linking, without changing
# the invoking shell or any CPU profile.
if [ -n "${GPU_MUSCLE_DRIVER_AUDIT:-}" ]; then
    driver_directory=$(/usr/bin/dirname -- "$0")
    /usr/bin/env PATH=/usr/bin:/bin /usr/bin/python3 \
        "$driver_directory/audit-driver.py" "$GPU_MUSCLE_DRIVER_AUDIT" "$@" || exit "$?"
fi
exec /usr/bin/env PATH=/usr/bin:/bin /usr/bin/gfortran-12 "$@"
