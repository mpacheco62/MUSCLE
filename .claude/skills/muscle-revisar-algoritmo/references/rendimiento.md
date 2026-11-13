# Velocidad en CPU y en GPU

## Qué medir

1. **Costo por llamada** de la rutina aislada: micro-banco con los casos de precisión (mayoría aleatorios), un hilo,
   Release, ns por llamada. Precedente: `benchmarks/math/bench_spectral.fypp`.
2. **Costo del kernel completo** de cada consumidor: es lo que decide. Una rutina 2x más lenta puede pesar 2 % en un
   return mapping, y una 30 % más rápida puede no notarse.
   - Materiales: `bench_batch_gpu_<material>_<kind>` (`benchmarks/gpu/bench_batch_gpu.fypp`, casos en
     `benchmarks/gpu/bench_cases.py`), en CPU con `--threads 1` y en GPU.
   - Cinemática finita: `benchmarks/math/bench_kin_kernel.fypp` (F = R(I + e), puntos elásticos y plásticos).
3. Build: `cmake -S . -B <dir> -G Ninja -DCMAKE_BUILD_TYPE=Release -DMUSCLE_BUILD_BENCHMARKS=ON` y los targets
   `bench_*`. Debug no sirve para tiempos.

## Protocolo de CPU (WSL2 en particular)

- En WSL2 una CPU virtual cae en núcleos P o E del anfitrión: corridas sueltas varían hasta 3.6x. Fija el núcleo
  con `taskset -c <n>`, **intercala** las variantes caso por caso (no todas las de una variante y luego las de
  otra), haz 7 rondas y toma el **mínimo** por caso.
- Incluye siempre un **caso de control** que no use la rutina (p. ej. `cut_vm_bilinear` para algo espectral): debe
  dar 1.00 entre variantes; si no, la medición no vale.
- Reporta relativo a la línea base (hoy = 1.00) y mediana/mínimo ≤ 1.02 como indicador de ruido.
- No corras builds ni GPU mientras mides CPU.

## Protocolo de GPU (Intel UHD 770, ifx SPIR64)

- Build: `-DMUSCLE_GPU=SPIR64 -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icx`, con
  `source /opt/intel/oneapi/setvars.sh` en el mismo shell. Compila con límite de memoria
  (`systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=0 cmake --build ...`).
- Ejecución, siempre con límites y la FP64 emulada (la imagen tiene código dp aunque midas sp):
  ```bash
  systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0 --quiet \
    env OMP_TARGET_OFFLOAD=MANDATORY OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 \
    timeout -k 30 1200 <bench> --points 100000,1000000 --repeat 5
  ```
  Nunca `ulimit -s unlimited`.
- **De a una corrida**: dos procesos en la GPU a la vez dan `DEVICE_LOST` o tiempos falsos. Si un script lanza
  corridas en segundo plano, asegúrate de que un fallo previo (p. ej. el lint) lo detenga.
- El banco imprime `CHECK ... <fallos>`: la última columna es el número de puntos de la muestra que difieren de la
  CPU más allá de `TOL_DEVICE_<k>`. **Un candidato rápido con fallos no es un candidato**: primero averigua por qué
  (`references/kernel_gpu.md`).
- Qué se puede medir: sólo los materiales que corren en la GPU (`tools/gpu_skip_tests.sh`); en Release, IGC entrega
  sin kernels de lote las imágenes de algunos materiales (D-019.8, «kernel ... not found»). Dilo en el informe en
  lugar de extrapolar.
- Los cocientes cambian cuando cambia la base: tras D-031 el kernel anisótropo pasó de 0.459 s a 0.145 s y un
  candidato que costaba 1.7x pasó a 4.5x. Mide siempre la base en la misma corrida.
- NVIDIA/AMD: sólo compilan aquí; no extrapoles tiempos de la GPU Intel a ellas.
