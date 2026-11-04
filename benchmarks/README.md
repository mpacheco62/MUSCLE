# MUSCLE benchmarks

Performance benchmarks. They are **not** part of the test suite: they are not registered in CTest
and are only built when requested.

```bash
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release -DMUSCLE_BUILD_BENCHMARKS=ON
cmake --build build-rel --target build_all_benchmarks
./build-rel/benchmarks/bench_closest_point_scalability --help
```

Use a Release build: the benchmarks print a warning when the build type is not `Release`. They use materials
of the default test registry (`tools/fypp/muscle_materials.py`): with `-DMUSCLE_MATERIALS=<path>` the registry
must also declare `cp_vm_swift` and `cp_dp_swift`.

| Benchmark | What it measures |
|---|---|
| `bench_closest_point_scalability` | Strong scaling (OpenMP) of the closest point return mapping through the batch API of the declared materials (`muscle_update_cp_vm_swift` / `muscle_update_cp_dp_swift` of the test registry, `dp`): time, speedup, efficiency and Karp-Flatt serial fraction per thread count, solver convergence summary, and a parallel == sequential check (exit code 1 on mismatch). Options: `--points`, `--threads`, `--repeat`, `--yield vm\|dp`, `--load tension\|compression`. |
| `bench_batch_gpu_<material>_<kind>` | Batch API on N integration points, GPU against CPU (Phase 8): `muscle_update_<material>_dev` with tangent on arrays of `omp_target_alloc`, one executable per case of `gpu/bench_cases.py` (small device images). With `MUSCLE_GPU` the arrays stay on the device: kernel time without transfers, device initialization (image load and JIT) and first call apart, host <-> device transfer times of the inputs and outputs. In a CPU build the same API runs the OpenMP batch with each thread count. A sample of points is checked against the point API on the host (`TOL_DEVICE_<k>`, exit code 1 on mismatch). Options: `--points`, `--threads`, `--repeat`, `--sample`. Lines `RESULT`, `TRANSFER`, `CHECK` are for scripts. |

## Batch API: Intel GPU against CPU (2026-11-03)

```bash
# GPU (ifx, Intel GPU): run under memory and time caps, as tools/check.sh does
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-bench-spir64 -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_GPU=SPIR64 \
    -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icx -DMUSCLE_BUILD_BENCHMARKS=ON
cmake --build build-bench-spir64 --target bench_batch_gpu_cut_vm_bilinear_sp
systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0 --quiet \
    env OMP_TARGET_OFFLOAD=MANDATORY OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 \
    timeout -k 30 1800 build-bench-spir64/benchmarks/bench_batch_gpu_cut_vm_bilinear_sp
# CPU (gfortran): the same benchmark in a build without MUSCLE_GPU
cmake -S . -B build-bench-cpu -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_BUILD_BENCHMARKS=ON
```

Hardware and software: Intel Core i9-13900K (8 P + 16 E cores, 32 threads), 32 GB, Intel UHD Graphics 770
(integrated, memory shared with the host), WSL2 (Linux 6.18); ifx 2026.2.0 (`-fiopenmp -fopenmp-targets=spir64`,
Level Zero, intel-opencl-icd 26.05.37020, IGC 2.28.4), gfortran 15.2 (`-O3`, OpenMP). Real kind `sp`; one increment
from the zero state, about 77 % plastic points, all converged; with tangent; minimum of 5 calls (3 for the Debug
GPU rows). Points/s = N / time of one batch call; speedup = GPU points/s over CPU points/s.

| Material | N | GPU (UHD 770) points/s | CPU 1 thread points/s | CPU 32 threads points/s | Speedup vs 1 thread | Speedup vs 32 threads |
|---|---:|---:|---:|---:|---:|---:|
| `cut_vm_bilinear` | 1e4 | 8.9e7 | 3.6e6 | 1.3e6 (*) | 25 | 67 (*) |
| | 1e5 | 2.9e8 | 2.8e6 | 8.7e6 | 102 | 33 |
| | 1e6 | 2.0e8 | 2.3e6 | 5.0e6 | 87 | 41 |
| `cut_yld2004_ludwik` | 1e4 | 2.0e6 | 2.1e5 | 3.7e6 | 9.7 | 0.54 |
| | 1e5 | 2.2e6 | 2.0e5 | 3.3e6 | 11 | 0.67 |
| | 1e6 | 2.2e6 | 2.0e5 | 2.9e6 | 11 | 0.78 |
| `pz_cpb06ex2_bilinear` (GPU Debug) | 1e4 | 4.1e4 | 1.4e5 | 1.7e6 | - | - |
| | 1e6 | 4.0e4 | 8.0e4 | 1.2e6 | - | - |
| `tm_jc_cpb06_ludwik` (GPU Debug) | 1e4 | 1.7e4 | 7.4e4 | 1.5e6 | - | - |
| | 1e6 | 2.8e4 | 1.2e5 | 1.7e6 | - | - |
| `pz_cha_yld2004_ludwik` | 1e4 | does not run | 7.5e4 | 9.8e5 | - | - |
| | 1e6 | does not run | 5.4e4 | 9.1e5 | - | - |

(*) 10^4 points with 32 threads take 7.5 ms in this run and 0.4 ms in another one: below 10^5 points the CPU
timings with all threads are noise (thread start-up, E cores, WSL2).

Device side (Release GPU rows):

- Device initialization (image load and JIT of the program, once per process): 2.3 s (`cut_vm_bilinear`), 19 s
  (`cut_yld2004_ludwik`); 43–49 s for the Debug images. The first batch call costs at most 2x a later one (0.2
  against 0.11 ms with 10^4 points), so the JIT happens at the image load, not at the first call.
- Transfers with 10^6 points (`omp_target_memcpy`, ~52 MB of inputs: strain and history at t_n; ~200 MB of
  outputs: history at t_n+1, stress, tangent, status, iterations): host -> device 9–24 ms (2–6 GB/s), device ->
  host 34–38 ms (5–6 GB/s). For `cut_vm_bilinear` they are 12x the kernel (4.9 ms): the batch only pays off
  with the data resident on the device (FFT solver, D-019.2).
- Check against the point API on the host (1000 points per N): statuses equal; stress, tangent and history
  within 43 eps (`TOL_DEVICE_sp` = 64 eps); up to 0.4 % of the points of the Newton solvers take one iteration more or
  less on the device (sp iterate on the stopping tolerance: FMA, order of the sums), within the same tolerance.
  The CPU batch equals the point API bit for bit.

Caveats:

- UHD 770: integrated GPU (32 EU, no FP64: dp is emulated and not measured here), memory shared with the
  CPU, so "transfers" are copies inside host RAM. A discrete GPU (NVIDIA/AMD, run by the user) will differ.
- Only the batch kernels IGC compiles right run on the GPU (`tools/gpu_skip_tests.sh`). In Release
  (`-O2`) the image of `pz_cha_yld2004_ludwik`, `pz_cpb06ex2_bilinear` and `tm_jc_cpb06_ludwik` loads without
  its batch kernels (IGC drops them silently with the dp emulation on; without it the module fails with "Double
  type is not supported", the sp and dp kernels share the image; with `-ze-opt-disable`: device lost). The Debug
  image (`-O0`) runs `pz_cpb06ex2_bilinear` and `tm_jc_cpb06_ludwik` right but slowly (no speedup computed: the
  Debug GPU rows are only a correctness check); `pz_cha_yld2004_ludwik` does not load in Debug either.
- The CPU rows are gfortran: an ifx CPU build of this benchmark does not link (undefined `omp_target_alloc`
  with `-fiopenmp` and no offload target).
