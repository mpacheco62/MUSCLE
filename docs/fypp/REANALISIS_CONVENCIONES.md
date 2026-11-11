# Re-análisis 7: consistencia de convenciones (2026-11-10, rama `fypp/f9-readability`, HEAD `840f486`)

Aspecto 7 de 15 del re-análisis posterior al plan R (nuevo; el ítem 3.6 de `ESTADO_OPTIMO.md` hizo un primer
inventario y corrigió lo interno). Pregunta: ¿se hacen las cosas parecidas de la misma manera, y se cumplen las
convenciones declaradas (`pages/conventions.md`)? Método: chequeos mecánicos sobre las plantillas y el código
generado (`build-maint-Debug`) más un inventario entre familias hecho por un agente de sólo lectura; revisé a
mano los puntos marcados ✔. No se cambió código.

## 1. Convenciones declaradas: cumplimiento (VERIFICADO)

| Convención | Resultado |
|---|---|
| Cabecera SPDX | Todos los archivos (el lint lo comprueba) |
| Literales con kind, sin `real(8)`/`dble`, sin `cpp` en plantillas | Cumple (lint) |
| `STATUS_*` en vez de enteros literales | Cumple: 192 usos, ningún literal |
| Kernels `pure` con `@:DEVICE()` tras `implicit none` | Cumple salvo `muscle_math_spectral_derivs`: 7 procedimientos sin `implicit none` propio (lo heredan del módulo) y `@:DEVICE()` tras la doc; más `linsolvers` y `tensor_io` (host) |
| Doc FORD `!!` en cada procedimiento y argumento | 7 procedimientos y 71 de 1 096 argumentos sin `!!` (sobre todo math: 25 con la doc en la misma línea) |
| Líneas ≤ 132 columnas «también después de expandir» | **No se cumple ni se comprueba**: 127 líneas en plantillas (45 de código, 82 de comentario) y 188 en el Fortran generado de `src` (todas comentarios, así que el compilador no falla) |
| Indentación de 4 espacios | `muscle_math_linsolvers` usa 2 (17 líneas) |
| Comentarios en inglés | `src/CMakeLists.txt:6` en español; la salida de `tools/dup_report.py` en español; los heredados `.f90` tienen comentarios en español (congelados) |

## 2. Inconsistencias entre familias (inventario; ✔ = verificado)

Pública = cambiarla rompe a los usuarios (nombres en `use`, genéricos, palabras clave).

| # | Inconsistencia | Alcance | Pública |
|---|---|---|---|
| I1 ✔ | **Validación de parámetros con 3 estrategias**: el `init` del modelo valida (visco); valida sólo el `mat_init` generado (`elas_linear`, Drucker-Prager, `MODEL_DOMAINS` de JC) y el `init` devuelve siempre `STATUS_OK`; **nadie valida** (hard, vp, os, von Mises, Hill48, CPB06, CPB06ex2, Yld2004) | ~16 modelos sin validación | Comportamiento |
| I2 | **Cuatro canales de error**: `status` entero (modelos, solvers, materiales); `ok` + `cond` (inversa de tensores); `OK_FLAG` + `KAPPA` (`M33INV`, `lu_solve`); `info` (`gauss_solve`, host) | tensores, math | Sí |
| I3 | **Nombres de math en tres estilos**: MAYÚSCULAS (`M33INV`…, `COND_INF`, `LU<n>SOLVE`), camelCase (`dEigenvalues_dTensor`), snake_case (`lu_solve`, `eigenvals`); sin prefijo de familia | ~20 | Sí |
| I4 | **Sufijo de dimensión**: yield es sólo 3D y lleva `_3D_` (20 específicos, salvo `yield_joint_evaluation_*`); kin y visco, también sólo 3D, no | ~25 | Específicos |
| I5 | **Nombres de derivadas**: `hard_dstress_dep`, `vp_dflow_depd`, `os_dflow_dsigma` nombran la variable; `yield_dstress_eq` no. `vp_newton_dstress_dep` no está en la lista de genéricos de la cabecera de `muscle_vp_laws` | 3 | Sí |
| I6 | **Archivo/módulo frente a prefijo**: `os` vive en `muscle_vp_*` y `muscle_vp_overstress_laws`; `elas` en `muscle_elasticity_*`; `elas_frozen` en `muscle_visco_frozen`; `muscle_yield_criteria` (las demás `*_laws`); `muscle_anisotropic_yield` sin prefijo. Es lo que confundió al recién llegado (aspecto 1) | 6 módulos | Sí (`use`) |
| I7 | **Constantes**: unas por kind en `muscle_kinds`, otras locales (`KIN_EIG_MIN_`, `MAXWELL_GAMMA_SERIES_`…); `DGAMMA_MIN` duplicada en dos cuerpos de solver; `1.0e-12` mágico en MRK y VA; prefijo `DP_` de Drucker-Prager que se lee como el kind `dp` | ~15 | Mixta |
| I8 | **Firmas de `*_init`**: `visco_state_init(state, n, status, strain_e)` no termina en `status`; `prm` `intent(inout)` en Drucker-Prager (documentado) | 2 | Sí |
| I9 ✔ | **Nombres de argumentos**: la tensión es `sigma` en yield y `stress` en kin/solvers/materiales; la cabecera de `muscle_kin_laws` documenta `kin_geometric_tangent(prm, state, sigma)` y los específicos usan `stress`; `st` frente a `state`; `dlambda` frente a `dgamma` | ~55 | Por palabra clave |
| I10 | Nombres de resultado (`res` salvo `result(elas)` en `visco_freeze_*`), helpers sin prefijo en `kin_finite` | ~20 | Mayormente interna |
| I11 | Estilo de doc en math (Markdown `##`, `@note`, doc en la línea de la declaración); referencias bibliográficas sólo en algunos modelos | ~45 | Interna |
| I12 | **Tests**: nombres en tres patrones (`test_<fam>_kinds`, `test_<fam>_laws_kinds`, `test_muscle_tensor_<fam>_kinds`) y 13 plantillas nuevas con el prefijo `test_muscle_` de los heredados; dos estructuras (`record()` con código por chequeo frente a `passed` + `stop 1`, los 19 `test_mat_*`); tolerancias locales y helpers `check_*` que duplican `muscle_test_utils` (7 en yield, 6 en math) | ~30 archivos | Interna |
| I13 | Herramientas: sólo `mutation_check.py` usa `argparse`; funciones CMake `add_muscle_*` frente a `muscle_add_*` | 6 + 8 | Interna (CMake la usan colaboradores) |

**Consistente** (según el inventario): tipos de parámetros `<prefix>_<model>_<k>` y específicos
`<prefix>_<func>_<model>[_dim]_<k>` en las 7 familias; `prm` primer argumento de toda función de modelo; todo
`*_init` de modelo termina en `status` (salvo I8); firma de `SOLVE` idéntica en los 5 cuerpos; cada módulo de
familia lista sus genéricos y tiene el párrafo «cómo agregar»; los 33 módulos de modelo tienen «Kernel API»;
la API generada de materiales sigue su tabla.

## 3. Conclusión

El **núcleo de modelos es consistente** (nombres, tipos, firmas, documentación): lo que la Fase 9 escribió a mano
quedó uniforme. Las inconsistencias se concentran en tres lugares: **math** (nombres heredados en tres estilos,
canales de error propios, doc con otro formato), **los nombres de archivo/módulo de algunas familias** (I6, que
ya costó fricción en la prueba del recién llegado) y **la validación de parámetros** (I1: la mayoría de los
modelos no valida nada; es también robustez, aspecto 10). Las convenciones de forma más incumplidas son el
límite de 132 columnas (nadie lo comprueba) y la doc de argumentos en math.

## 4. Propuestas (no aplicadas)

**Internas, sin decisión de API (≈ 0.5 día; snapshot idéntico salvo comentarios):**
- V1. **Hecho (2026-11-10, `0989d80`)**: regla `line-length` sobre toda fuente (comentarios incluidos) y las 128 líneas corregidas; el generado mantiene su regla sin comentarios. Antes: Lint de 132 columnas sobre el generado de `src` (y corregir las 127 + 188 líneas).
- V2. **Hecho (2026-11-11)** (y `gauss_solve` con `A`, `b` `intent(in)`, decisión del usuario). `muscle_math_spectral_derivs` y `linsolvers` al estilo de los modelos (`implicit none`, `@:DEVICE()`,
  `!!` en la línea siguiente, 4 espacios); comentarios en inglés en `src/CMakeLists.txt` y `dup_report.py`.
- V3. **Hecho (2026-11-11).** Cabeceras: `kin_laws` (`stress`), `vp_laws` (`vp_newton_dstress_dep`); `DGAMMA_MIN` en un solo lugar;
  `1.0e-12` de MRK/VA como constante por kind.
- V4. **Hecho de otra forma (2026-11-10, D-025)**: nombres `test_muscle_<x>[_<tema>]` espejo de `src`. Propuesta original: Tests: unificar nombres de los `*_kinds` y quitar el prefijo `test_muscle_` de las plantillas nuevas
  (cambia nombres de ctest); helpers locales a `muscle_test_utils`. Puede ir junto con Q2 (retiro de heredados).

- V9. **Hecho (2026-11-11, D-026).** Tarea acordada con el usuario (2026-11-10):** indentar las directivas de control de FYPP (`#:if`, `#:for`,
  `#:set`, `#:assert`, sus cierres) al nivel del código que controlan; `#:include` y `#:def` siguen en la columna 0.
  FYPP 3.1 lo acepta (probado) y borra la línea de la directiva, así que el generado no cambia (`gen_snapshot.sh`
  idéntico). Hoy 1 467 directivas en columna 0 en `src`. Commit sólo de formato, decisión en `DECISIONES.md` y regla
  en `conventions.md`/`fypp-templates.md`; con un script que indente según la línea de código siguiente y revisión
  a mano de los casos raros (directivas dentro de listas de continuación).

**Públicas (requieren decisión; alias de compatibilidad durante una versión):**
- V5. **Descartado por ahora** (el usuario descartó el punto de robustez, 2026-11-10). I1: validar en el `init` de cada modelo y que `mat_init` lo reuse (cambia comportamiento: parámetros
  hoy aceptados devolverían `STATUS_INVALID_PARAMETER`). **Recomendado**, junto con el aspecto 10.
- V6. **Hecho (2026-11-11, D-027, sin alias).** I6: renombrar módulos de `os`/`elas`/`yield` a `muscle_<prefix>_*` con módulos alias. Beneficio en
  comprensión; costo en usuarios con `use`. Alternativa barata: la tabla familia → archivo → módulo (C1).
- V7. I3/I2: nombres de math en snake_case con alias en mayúsculas; un solo par `ok`/`cond`.
- V8. I4/I5/I9: sufijo de dimensión, nombres de derivadas y `sigma`/`stress`: **no recomiendo** cambiarlos
  (romperían llamadas por palabra clave y específicos con poco beneficio); documentarlos en `conventions.md`.
