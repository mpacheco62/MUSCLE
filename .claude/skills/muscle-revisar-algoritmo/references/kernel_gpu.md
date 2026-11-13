# Tamaño del kernel y riesgos de dispositivo

## Medir

- `.claude/skills/muscle-revisar-algoritmo/scripts/igc_kernel_report.sh <ejecutable> [regex]`: corre el ejecutable
  una vez en la GPU con los volcados de IGC y da por kernel:
  - `private`: memoria privada por hilo. **32768 en un kernel SIMD8 = 4 KB por work item, el valor por defecto** que
    IGC usa cuando no puede acotar la pila y el kernel no tiene stack calls: casi siempre un temporal de tamaño
    variable (D-030). Un valor grande calculado (cientos de KB) es normal con stack calls;
  - `spill`: bytes de registros desbordados; sube con la presión de registros (más variables vivas, arreglos
    privados indexados en tiempo de ejecución);
  - `stackcalls`, `calls`: funciones que IGC no metió en línea;
  - `fp64emu`: llamadas a la emulación FP64 (`__igcbuiltin_dp_*`); en un kernel sp deben ser 0;
  - `dyn_alloca`: `alloca` de tamaño dinámico en la imagen y en qué funciones.
  Si el kernel pierde el dispositivo, lo repite con `IGC_StackOverflowDetection=1`.
- `tools/kernel_resources.sh --spir64 <build>` (GRF y spill por material) y `--amdgcn <build>` (VGPR, SGPR, pila por
  árbol de llamadas): para comparar materiales o variantes.
- Comparar variantes: mismo ejecutable con cada variante del worktree; cuenta también las líneas del ISA o del vISA
  si el tamaño importa.

## Riesgos conocidos (revísalos en todo kernel nuevo o cambiado)

| Riesgo | Síntoma | Causa y arreglo |
|---|---|---|
| `scale()` / `exponent()` en un kernel sp | la imagen pide FP64; resultados erróneos en GPU Intel cuando sólo uno de dos hilos SIMD8 de una pareja de EU ejecuta una rama | ifx los compila con dp emulado; IGC compila mal esas llamadas en EU fusionadas (D-031, `docs/fypp/igc_mre/`). Usar `fp_scale`/`fp_exponent` |
| Temporal de tamaño variable (`sum(..., dim=)`, expresión de arrays sobre forma asumida) | `DEVICE_LOST` también con 1 work item; `private` = 32768; «Stack overflow detected» | IGC reserva 4 KB por work item sin stack calls (D-030). Bucles explícitos (`NORM_INF`) |
| Literal o intrínseca dp en un kernel sp | `fp64emu` > 0 | revisar kinds de literales y conversiones |
| Llamadas a la emulación FP64 dentro de ramas divergentes (dp en la GPU Intel) | resultados distintos a CPU que cambian con la entrada de los work items vecinos | bug de IGC (D-031); diagnóstico con las variables de abajo; los kernels dp de material están excluidos (`tools/gpu_skip_tests.sh`) |
| Kernel grande en Release | «kernel ... not found in the program» al cargar | IGC descarta el kernel sin error de build (D-019.8); medir en Debug o reducir el kernel |
| Arreglos privados indexados con valores de tablas `PARAMETER` | más `alloca`, más spill | aceptable si mide bien; con literales el compilador los vuelve registros |
| Desborde intermedio (`theta*theta`, `x**2` de algo cerca de `huge`) | Inf en GPU, aborto con trampas en CPU | reformular sin formar el valor grande (D-028) |

## Diagnóstico (no son soluciones de producción)

Con el mismo ejecutable, una a la vez, y comparando con la corrida normal:

- `IGC_DisableEuFusion=1`: despacha sin fusión de EUs. Si corrige el resultado, es un riesgo de EU fusionadas
  (el código máquina es el mismo: sólo cambia `require_disable_eufusion` en el `.zeinfo`).
- `IGC_EmulationFunctionControl=1`: mete en línea la emulación FP64 (saca las llamadas).
- `IGC_StackOverflowDetection=1`: avisa de desborde de la pila del dispositivo.
- `IGC_ForcePerThreadPrivateMemorySize=<bytes por work item>`: si con más memoria privada pasa, era la pila.
- `LIBOMPTARGET_INFO=4`, `LIBOMPTARGET_DEBUG=1`: lanzamientos y carga de la imagen.
- `IGC_ShaderDumpEnable=1 IGC_DumpToCustomDir=<dir>`: los volcados (`.asm`, `.zeinfo`, `*_optimized.ll`) que lee el
  script.
- Para aislar si un fallo es de un kind: un build con sólo sp exige parches (los alias históricos apuntan a dp,
  D-003); antes de hacerlo, mira `fp64emu` del kernel sp.

Un fallo que cambia al apagar la fusión, al cambiar el ancho SIMD o al reordenar código sin cambiar la aritmética
es del compilador, no del algoritmo: redúcelo a un reproductor sin la librería (precedente: `docs/fypp/igc_mre/`) y
regístralo en `PROBLEMAS.md`.
