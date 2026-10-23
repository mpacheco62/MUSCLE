# Etapa 4: integración Fypp y especialización CPU

La integración CPU/NVIDIA está validada. Fypp 3.3 genera fuentes que se compilan
en la biblioteca real; la API CPU original se conserva. La especialización
concreta CPU comparte los métodos numéricos de Closest_point con la API original.
El solver constitutivo completo en GPU corresponde a la etapa 5.

## Cambios

- `MUSCLE_USE_FYPP=ON` por defecto: los módulos públicos de Elasticidad lineal,
  Swift, Von Mises y Closest_point pasan por Fypp conservando sus nombres.
  El resto de las fuentes se compila como antes. `MUSCLE_USE_FYPP=OFF` conserva
  el build CPU sin depender de Python y omite las nuevas especializaciones.
- Fypp real está incluido sin modificaciones en `tools/fypp`, con su licencia,
  commit release 3.3 y SHA256 de fuente/licencia verificados al configurar.
  Se necesita Python >=3.10; no hay descargas durante el build.
- El nuevo módulo `muscle_solver_closest_point_lsvm` expone `Closest_point`
  con componentes concretos `Elasticity_linear`, `Swift_hardening` y `VonMises`.
  Sus componentes no son allocatables polimórficos. Los métodos iterativos,
  solve y tangentes se comparten mediante `closest_point_methods.inc`.
  Esta especialización sigue siendo CPU: no se declara el solver GPU validado.
- `muscle_kernels_cpu` y `muscle_kernels_target` se generan de una misma plantilla.
  Las tres expresiones Swift y la inversa LU 6×6 se comparten con el código
  original. Los parámetros Swift son un tipo concreto sin asignación dinámica.
- `muscle_target` concentra las opciones Target de compilación y enlace.
  El perfil NVIDIA exige el driver GCC12 con PATH restringido, NVPTX y `sm_80`.
  La biblioteca CPU se compila sin flags NVPTX. Los perfiles y targets separan
  fuentes generadas y módulos `.mod`; no se intercambian módulos entre compiladores.
- Los depfiles Fypp, los INCLUDE Fortran y las dependencias explícitas de targets
  evitan carreras en el escaneo Fortran y regeneran los consumidores afectados.
  Un build sin cambios no realiza trabajo.

## Uso

La API CPU existente no cambia. Para seleccionar la combinación concreta:

```fortran
use muscle_solver_closest_point_lsvm, only: Closest_point_lsvm => Closest_point
! elasticity, hardening y yield son los tipos concretos originales.
type(Closest_point_lsvm) :: solver
! call solver%init(elasticity, hardening, yield)
! call solver%solve(strain, history, status, iterations)
```

Desde la raíz del repositorio MUSCLE:

```sh
cmake --preset fypp-cpu-release
cmake --build --preset fypp-cpu-release
ctest --preset fypp-cpu-release
```

Los presets `fypp-cpu-debug` y `native-cpu-release` verifican Debug y el build
sin generación respectivamente. Para NVIDIA, desde un entorno con GPU accesible:

```sh
cmake --preset fypp-nvidia-gcc12
cmake --build --preset fypp-nvidia-gcc12
ctest --test-dir build/fypp-nvidia-gcc12 --output-on-failure -V
```

La batería reproducible de los cuatro perfiles guarda logs, XML y resumen:

```sh
python3 validation/scripts/validate_stage04.py --output build/stage04-validation
```

Devuelve 0 solamente si se validan los cuatro perfiles y la ejecución de los
nuevos núcleos en GPU real. `--cpu-only` valida CPU y deja `stage4_complete=false`.
Un entorno restringido sin dispositivo puede devolver 2 aunque otro entorno
haya validado la GPU. Los resultados fallidos se conservan por separado.

## Evidencia de validación (2026-10-23)

| Perfil | Resultado |
| --- | --- |
| CPU nativa sin Fypp, GCC11 Release | 51/51 |
| CPU Fypp, GCC11 Release | 53/53 |
| CPU Fypp, GCC11 Debug | 53/53 |
| NVIDIA GCC12, suite completa | 53/53: 52 comprobaciones CPU y 1 prueba de núcleos Target |

La comparación usa 768 estados y 70 campos numéricos: historia comprometida y
candidata, tensión, deformación plástica, tangente y residuos. Los identificadores,
estados e iteraciones coinciden exactamente. Las API original y concreta y la
referencia almacenada coinciden exactamente en todos los campos en los tres
builds Fypp: diferencia máxima **0**. La prueba también exige convergencia y
que `solve` preserve el estado comprometido antes de `commit`.

La nueva prueba Target usa 257 puntos y la API CPU original como referencia
para Swift, sus dos derivadas y LU 6×6, incluyendo matrices con una fila cero.
En GPU registra `omp_num_devices=1`, `calculation_host_points=0`, política
MANDATORY y `PASS_GENERATED_KERNELS_DEVICE`. La máxima diferencia de Swift es
`3.5527136788005009e-15`. Las inversas respetan `atol=1e-12`; sus estados
coinciden exactamente. Las tolerancias se fijaron antes de ejecutar la GPU.

También se verificaron:

- Reconstrucción exacta de las fuentes originales al expandir los INCLUDE
  extraídos de Closest_point y Swift.
- Regeneración de todos los consumidores tras cambiar la fecha de los cuerpos
  compartidos, y build posterior sin trabajo.
- Instalación CPU y consumo externo mediante `find_package`, importación de
  módulos nuevos, enlace y ejecución.
- Rechazo de una copia Fypp con SHA256 alterado y del fallback host forzado.

[Resumen completo](../../validation/evidence/stage04/summary.json),
[prueba GPU](../../validation/evidence/stage04/nvidia-kernels.log),
[comparación CPU](../../validation/evidence/stage04/fypp-cpu-release-comparison.json),
[dependencias](../../validation/evidence/stage04/dependencies.json) y
[controles negativos](../../validation/evidence/stage04/negative-controls.json).

AMD y CPU-001 (Drucker–Prager en tracción) permanecen pendientes. No se modifican
algoritmos para resolver construcción ni se afirma aceleración sin mediciones.
La etapa 5 deberá implementar y validar el recorrido constitutivo completo en
GPU, incluyendo convergencia, historia, tangente y commit/rollback por punto.
