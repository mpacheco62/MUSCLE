# Análisis: sobreingeniería (2026-11-04, rama `fypp/f9-readability`)

Aspecto 3 de 4 (1: `ANALISIS_COMPRENSION.md`, 2: `ANALISIS_REPETICION.md`; queda mantenibilidad).
Criterio: un mecanismo es **sobreingeniería** si su costo (código, conceptos, tiempo de CI, pasos para
quien lo toca) supera lo que aporta hoy, aunque haya tenido sentido cuando se creó. No lo es lo que sirve
a un objetivo declarado del proyecto (sp/dp, 2D/3D, GPU), aunque sea grande. Método: medir consumidores de
cada mecanismo (grep, análisis estático de `muscle_meta.py`, lectura de campos); los veredictos son mi
juicio sobre esas cifras. VERIFICADO = medido; ESTIMADO = mío.

## 1. Veredicto por mecanismo

| # | Mecanismo | Tamaño | Evidencia (VERIFICADO) | Veredicto |
|---|---|---|---|---|
| 1 | **Registro de funciones de familia en Python** (`MODEL_FAMILIES` + `ModelFamily/ModelFunc/ModelArg/ModelJoint/ModelComposition`) | 222 líneas de datos + 5 clases + ayudantes | Nació para generar declaraciones (retirado en la Fase 9). Hoy producción lee sólo `prefix`, `module`, `family_module`, `models`. Campos leídos únicamente en tests: `ModelFamily.args/fd_scale/fd_lower/support`; `ModelFunc.result/deriv/extra/default/outs/opts`; `ModelJoint.outs/optional`; `ModelParam.default`. `ModelArg`: ningún campo se lee por atributo. `Material.os`: no se lee. | **Sobra en gran parte.** Lo que queda lo usan `tests/gpu/target_cases.py` y `test_family_target.fypp` para generar tests de familia, duplicando en Python las firmas que ya están en Fortran. |
| 2 | `MODEL_PARAMS` (188 líneas) + `MODEL_INITS` (57) | 245 líneas | Alimentan la tabla de props de `mat_init`; el test las compara con los tipos escritos a mano (nombre, orden, tipo, defecto) | **Se justifica** (es el dato real de cada modelo) |
| 3 | **Generación de materiales desde un registro** | 4 plantillas, 1 147 líneas → 80 módulos, 148 886 líneas | La librería (`muscle`) compila en ~5 s Debug / ~9 s Release; el registro cubre cada modelo por solver | **Se justifica**: Fortran no puede declarar un módulo por entrada |
| 4 | **Puntos de entrada por material** (`mat_update`, `muscle_point_`, `muscle_update_`, `_dev`, `mat_pack/unpack`, `mat_rotate_hsv`, `mat_step`, `mat_init`) | 9 nombres | Cada uno aparece en 1–18 archivos de tests/benchmarks y en 1–4 páginas de docs | **Se justifica** (ninguno sin consumidor) |
| 5 | **Macros de `include/`** | 12 macros, 460 líneas | `@:DEVICE` 401 usos; el resto 2–6 usos cada una | **Se justifica**; tamaño pequeño |
| 6 | **API histórica de tipos (`*_tbp`, `class`)** | 16 bindings en 6 familias | Todos con ≥ 2 usos en tests/docs/benchmarks (`init` 888, `is_approx` 164, `inv` 87) | **Se justifica** por compatibilidad; `class` está confinado ahí |
| 7 | **Infraestructura GPU y de verificación** | 1 459 líneas (`check.sh` 399, `audit_kernels` 193, `kernel_resources` 190, `MuscleGpu`/`MuscleFypp` 363, snapshots y métricas ~230) | Tres builds GPU y dos auditorías en cada `check.sh`; lista de exclusiones por bugs de IGC (4 patrones) | **Se justifica** por el objetivo multi-vendor. **Costo recurrente:** `check.sh` completo tarda del orden de decenas de minutos (no cronometrado) y no tiene modo intermedio entre `--lint-only` y todo |
| 8 | **Tests y goldens** | tests 29 395 líneas (≈ 2 : 1 frente a `src`), referencias 8 839 líneas en 36 archivos, infraestructura de goldens 1 937 líneas | La política «todo material del registro tiene golden» obliga a 3 pasos extra por material (lo mostró la prueba del recién llegado) | **Dudoso**: la rigurosidad es valiosa, pero agrega fricción y no hay forma de saltarla |
| 9 | **Proceso documental** | `docs/fypp` 2 089 líneas (+ 3 análisis), 384 referencias `D-0xx` en `src/` | DECISIONES 1 204 líneas, PLAN 709 | **Dudoso**: útil para el autor, pesado para otro lector (ver análisis 1) |
| 10 | **Capacidades sin consumidor interno** | 2D: tensores + `elas_linear` 2D + `visco_frozen` 2D; UDTIO 180 líneas; esqueletos de UMAT/USERMAT/umat41 | Ningún material usa 2D (solo 3D); UDTIO está desactivado por defecto con gfortran < 16.3; los esqueletos FE están rotulados «no compilado ni probado» | **No es sobreingeniería**: responden a objetivos declarados. Son **carga a vigilar**: código mantenido sin uso propio |
| 11 | **Soporte sp** (tolerancias por kind, pisos de redondeo, `TOL_*_k`) | disperso | Requisito explícito del usuario | **Se justifica** |

## 2. Lo central

1. **El único exceso claro es el registro de funciones de familia (#1).** Es la parte del diseño de las
   Fases 6–8 que sobrevivió a la retirada de las macros de contrato: guarda en Python información (argumentos,
   resultados, derivadas por defecto, salidas opcionales) que ya está escrita en Fortran y que sólo leen
   los generadores de tests. Cualquier cambio de firma exige editarla además del Fortran, y nada comprueba que
   coincidan salvo la compilación de esos tests.
2. **El resto del andamiaje se justifica con consumidores reales.** Ninguna macro, función de `muscle_meta.py`
   ni punto de entrada resultó sin uso (0 funciones muertas). El código de producción de `muscle_meta.py` son
   1 586 líneas alcanzables desde plantillas/cmake/tools; sólo 85 sirven únicamente a tests.
3. **La complejidad pesa más en el proceso que en el código**: política de goldens por material (#8),
   `check.sh` sin modo intermedio (#7) y volumen documental (#9).
4. **Vigilar (#10):** lo que no se usa dentro del proyecto sólo está protegido por sus tests; si no hay
   planes de usar 2D en materiales, ese es el candidato a recibir menos inversión (no a borrarse).

## 3. Propuestas ordenadas

| Prioridad | Cambio | Ahorro / efecto (ESTIMADO) | Riesgo |
|---|---|---|---|
| 1 | Reducir el registro de familias (#1) a lo que produce nombres (prefijo, módulo, modelos) y reescribir `target_cases.py`/`test_family_target` con una tabla pequeña o a mano | ~250–300 líneas de Python y 4–5 clases menos; una fuente menos que mantener | Bajo (afecta sólo tests y datos; el generado de `src` no cambia) |
| 2 | `check.sh --quick` (lint + build CPU + ctest sin GPU) y dejar el completo para antes de integrar | Ciclo de desarrollo mucho más rápido; sin pérdida de cobertura | Muy bajo |
| 3 | Escape controlado para la política de goldens (p. ej. `golden=False` explícito y listado en un informe) o un script que genere el golden de un material nuevo | Un material nuevo deja de exigir 3 ediciones manuales | Bajo; debe quedar visible y revisable |
| 4 | Recortar referencias `D-0xx` en `src/` (análisis 1) | Menos ruido | Nulo |
| 5 | Decidir qué hacer con la carga sin uso (2D en materiales, UDTIO, esqueletos FE) | Evita mantener sin dueño | Decisión del usuario |

No recomiendo tocar los puntos 2–6 y 11 de la tabla: tienen consumidores y objetivos que los sostienen.

## 4. Pendiente de decidir por el usuario
- Si se acomete la reducción del registro de familias (prioridad 1).
- Si se agrega `check.sh --quick` (prioridad 2).
- Si se acepta una vía de escape para los goldens de materiales nuevos (prioridad 3) o la política se mantiene.
- Qué se hace con las capacidades sin consumidor interno (prioridad 5).
