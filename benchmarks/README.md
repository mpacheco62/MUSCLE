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
