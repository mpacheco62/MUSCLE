# Reproductor del error de IGC con emulación FP64 en EU fusionadas

Reproductores autocontenidos (Fortran + OpenMP offload, sin la librería) de un error de compilación de IGC en la
Intel UHD Graphics 770 (Gen12LP, sin FP64 nativo), encontrado con el prototipo de Jacobi de `PLAN_ESPECTRAL.md`
(D-031). `REPORT.md` es el informe para Intel (en inglés): síntoma, regla exacta, controles, análisis del código
máquina y entorno.

| Archivo | Contenido |
|---|---|
| `mre_igc_scale_fp64emu.F90` | Reproductor original (459 líneas): return mapping con criterio anisótropo; 48 de 256 puntos mal con la entrada natural |
| `mre_igc_scale_fp64emu_patterns.F90` | El mismo con entradas diseñadas por bloque de 16 work items (`-DPATTERN=1..7`, `-DDUMP`) |
| `mre_igc_reduced.F90` | Versión reducida (280 líneas): 128 de 256 mal con un punto plástico por bloque de 16 |
| `REPORT.md` | Informe para Intel |

## Resumen

- **Síntoma:** cuando exactamente uno de dos hilos SIMD8 vecinos (una pareja de EU fusionadas) ejecuta la rama
  divergente, los 8 carriles del otro hilo devuelven un valor erróneo (el contador de iteraciones 50 en lugar de 0;
  la tensión es correcta).
- **Condiciones:** llamadas a funciones de emulación FP64 dentro de la rama (ifx compila `scale()` sobre `real32`
  con operaciones dp), arreglos locales indexados con índices leídos de tablas `PARAMETER`, y `-O1` o más.
- **Causa:** el mismo código máquina es correcto con la fusión de EUs desactivada (`IGC_DisableEuFusion=1` sólo
  agrega `require_disable_eufusion` al `.zeinfo`): es un riesgo de las EU fusionadas que el workaround de IGC
  (`-noMaskWA`, `-fusedCallWA`) no cubre del todo.
- **Versiones:** se reproduce con NEO 26.05 / IGC 2.28.4 y con NEO 26.40.40196.8 / IGC 2.42.6 (lo último publicado
  al 2026-11-12). Issues relacionados: intel/intel-graphics-compiler#397 y #413.
- **Rodeos:** `IGC_EmulationFunctionControl=1` o `IGC_DisableEuFusion=1`. En la librería, `fp_scale`
  (`muscle_math_fp`, D-031) reemplaza `scale()`: los kernels sp ya no llaman a la emulación.

## Cómo correrlo

```bash
source /opt/intel/oneapi/setvars.sh
ifx -O2 -fiopenmp -fopenmp-targets=spir64 -fpp mre_igc_reduced.F90 -o mre_red
OMP_TARGET_OFFLOAD=MANDATORY OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 ./mre_red   # 128 malos
OMP_TARGET_OFFLOAD=DISABLED ./mre_red                                                          # 0 malos
```

En esta máquina, siempre con límite de memoria: `systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0`.
