# CPU/NVIDIA feature probes (stage 3)

See [the validated baseline and reproduction commands](../../docs/gpu/stage03.md).
The probes are a separate CMake project. They do not offload the complete
MUSCLE solver. NVIDIA tests require GPU execution and reject host fallback.
