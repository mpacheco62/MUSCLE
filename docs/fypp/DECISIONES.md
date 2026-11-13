# Registro de decisiones (migración FYPP)

Formato: cada decisión tiene estado (**Propuesta**, **Aceptada**, **Reemplazada por D-xxx**),
contexto, decisión, alternativas descartadas y consecuencias. Las decisiones aceptadas son
vinculantes para Claude; las propuestas deben confirmarse con el usuario antes de implementar
algo que dependa de ellas.

---

## D-001 — FYPP como generador único
**Estado:** Aceptada (2026-10-28, implementada en la Fase 0)
**Contexto:** Se necesita generar variantes por kind y dimensión. El proyecto ya usa `-cpp`.
**Decisión:** Usar FYPP (Python, mismo enfoque que fortran-lang/stdlib). Los `.f90`
generados viven en el build, no en git. `cpp` sólo para flags de compilación (`USE_MKL`,
`ENABLE_UDTIO`) mientras existan.
**Alternativas:** macros `cpp` + `include` (limitadas, sin bucles); PDTs (ver D-002);
generar con script Python propio (reinventar FYPP).
**Consecuencias:** FYPP pasa a ser dependencia de build.

## D-002 — Un tipo por kind con sufijo `_sp`/`_dp`
**Estado:** Aceptada (2026-10-29, implementada en el piloto de la Fase 2)
**Contexto:** Fortran no tiene genéricos de tipo; los PDT (`type :: t(k)`) tienen bugs en
gfortran y soporte incierto en compiladores GPU.
**Decisión:** Tipos concretos `ten_3D2Osym_sp` / `ten_3D2Osym_dp`, procedimientos con sufijo,
genéricos que agrupan ambos.
**Consecuencias:** Más símbolos; resolución de genéricos por kind; sin coste en tiempo de ejecución.

## D-003 — Nombres históricos = doble precisión
**Estado:** Aceptada (2026-10-29, implementada en el piloto de la Fase 2; mecanismo en D-011)
**Decisión:** `muscle_tensors` re-exporta `ten_3D2Osym => ten_3D2Osym_dp` (y así para
todos los tipos). Tests y código de usuario existentes compilan sin cambios.

## D-004 — Conservar layouts; unificar 2D/3D con metadatos
**Estado:** Propuesta
**Contexto:** Los tipos 2D tienen órdenes de almacenamiento propios (p.ej. `ten_2D4O3sym`
con 10 componentes en orden `xxxx,yyyy,zzzz,xyxy,xxyy,yyzz,zzxy,xxzz,yyxy,xxxy`).
**Decisión:** Los layouts se describen en `tools/fypp/muscle_meta.py` copiados del código
actual; las plantillas generan código para cualquier layout. Cambiar un layout es una
decisión aparte.

## D-005 — Objetivo GPU
**Estado:** Aceptada (2026-10-25, decisión del usuario)
**Decisión:** OpenMP offload estándar. Primero gfortran (`-foffload=nvptx-none`,
`amdgcn-amdhsa`), después ifx (`-fopenmp-targets=spir64`). Objetivo: NVIDIA, AMD e Intel
con la misma fuente. La macro `@:DEVICE()` emite `!$omp declare target`.
**Descartadas:** OpenACC/nvfortran y CUDA Fortran (atados a un vendor).
**Consecuencias:** sólo características de OpenMP soportadas por los tres backends; el
prototipo de la Fase 2 debe validar tipos derivados con operadores sobrecargados en
offload con gfortran. gfortran 15.2 de la máquina tiene offload configurado, pero faltan
los paquetes de backend (`gcc-15-offload-nvptx`, `gcc-15-offload-amdgcn`).

## D-005b — Todo en simple y doble precisión
**Estado:** Aceptada (2026-10-25, decisión del usuario)
**Decisión:** tensores, math, modelos y solvers se generan para `sp` y `dp`.

## D-006 — Estrategia de despacho estático
**Estado:** Propuesta
**Decisión:** (a) Unión etiquetada por familia de modelo (`integer :: id` + parámetros en
array fijo) con `select case` **generado** desde un registro FYPP; (b) monomorfización
opcional de solvers para combinaciones declaradas. La API OO actual queda como fachada CPU.
**Alternativas:** mantener `class`/`deferred` (no compila para GPU); monomorfizar todo
(explosión combinatoria elasticidad × fluencia × endurecimiento × kind × dim).
**Nota (2026-10-25, experimento con gfortran 15.2 -O3):** con `id` en tiempo de ejecución,
gfortran inlinea el dispatcher y los modelos y hace *loop unswitching*: el `select case` se
evalúa una vez fuera del bucle y genera un bucle vectorizado por modelo (más tamaño de código).
Con `id` constante (instancia monomorfizada) el `select case` desaparece. Ambas cosas
dependen de heurísticas de inlining, así que **el `select case` se coloca en el nivel más
externo posible** (una vez por material/bloque de puntos, antes del Newton), y por debajo
todo es monomórfico. Verificar en Fase 7 con `-fopt-info-optimized` y en el build de offload.
**Límite a la monomorfización (2026-10-25):** monomorfizar todas las combinaciones explota
(p.ej. solver viscoplástico: 6 criterios × 6 leyes vp × 4 endurecimientos × 2 kinds × 2 dims
= 576 instancias). Regla: se monomorfiza sólo el eje caro (criterio de fluencia) junto con
kind y dim (inevitables); endurecimiento y leyes viscoplásticas (escalares) van por
`select case`. Opción CMake `MUSCLE_MODELS` para compilar sólo los modelos que se usan, y CI
que compile la matriz completa. Riesgos a vigilar: tiempo de compilación, tamaño de la
imagen de offload (JIT de PTX/SPIR-V al primer lanzamiento), nombres ≤ 63 caracteres.
**En GPU (2026-10-28):** un kernel con todos los criterios en un `select case` reserva los
registros de la rama más cara (p.ej. Yld2004) para todos los hilos, lo que baja la
ocupación incluso para von Mises y puede provocar *spills* a memoria local. Por eso el
`select case` del criterio se hace **en el host**, fuera de la región `target`, y lanza un
kernel monomorfizado por criterio (cada uno con su propio consumo de registros). Dentro
del kernel sólo quedan `select case` de leyes escalares (endurecimiento, vp). Los
parámetros de material viven en un array `mats(nmat)` mapeado explícitamente una vez
(`target enter data map(to: ...)`) y se pasan por referencia: nunca copia local ni
`firstprivate` por hilo. Los puntos se agrupan por material para evitar divergencia.
Medir registros/spills por kernel en Fase 8 (`ptxas -v` en NVIDIA; metadatos `.amdhsa_*`
en AMD).

## D-007 — Cómo se escribe un modelo (y qué pasa con la API OO)
**Estado:** Aceptada (2026-10-28) en cuanto al contrato de modelo y al retiro de la OO. Sin fachada
OO: los pocos usuarios externos se reescriben contra la nueva API. Las uniones (punto 3)
quedan como opcionales, ver D-008.
**Contexto:** Hoy la OO resuelve dos cosas con el mismo mecanismo (`type, abstract` +
`deferred`): (1) *cómo se escribe* un modelo (extender `Base_yield_critera` e implementar
`stress_eq`) y (2) *cómo lo llama* el solver (despacho dinámico por vtable). (1) es cómodo;
(2) impide el offload portable y la inlining entre solver, modelo y tensores.
**Decisión propuesta:** separar ambas cosas.
1. **Contrato de modelo** (para quien escribe un modelo): un módulo con un `type` de
   parámetros y procedimientos de nombre fijo (`stress_eq`, opcionales `dstress_eq`,
   `ddstress_eq`), escritos con nombres cortos `real(wp)`, `type(sym2)`, `type(sym4)`.
   Una macro FYPP (`#:call yield_model(...)`) instancia el cuerpo para cada kind (y dim),
   inyecta los `use ... only: sym2 => ten_3D2Osym_dp, ...`, añade `@:DEVICE()` y genera
   las derivadas numéricas que falten. El autor no ve sufijos ni `${k}$`.
2. **Registro**: una línea por modelo (`YIELD_MODELS = [..., 'hill48']`).
3. **Despacho generado**: `type :: yield_criterion_${k}$` (unión etiquetada con un
   componente por modelo + `id`) y genéricos `stress_eq(y, s)`, `dstress_eq(y, s)` con
   `select case` generado. Los solvers se escriben una vez contra la unión; las
   combinaciones calientes se pueden monomorfizar (D-006).
4. **API de usuario**: tipos concretos con constructores y procedimientos genéricos
   libres: `y = hill48(F=..., G=...)`, `f = stress_eq(y, s)`,
   `call return_map(solver, y, h, strain, data)`. **Corrección (2026-10-25):** no se
   puede usar la sintaxis `obj%metodo()` con `type`: Fortran exige que el argumento
   *passed-object* de un type-bound procedure sea `class(...)` (gfortran: "Non-polymorphic
   passed-object dummy argument"), y los tipos `sequence`/`bind(c)` no admiten type-bound
   procedures. En kernels se usan genéricos libres; `obj%metodo()` sólo en la fachada CPU.
5. **Tipos abstractos actuales** (`Base_*`): se retiran como mecanismo de extensión. Sólo
   se mantiene una fachada de compatibilidad si hay código externo (UMAT) que los usa,
   marcada como obsoleta, y se elimina al final (Fase 9).
**Alternativas:** (a) mantener OO completa en paralelo a la versión estática — doble
mantenimiento y dos caminos de resultados; (b) eliminar todo rastro de "objetos" — más
simple internamente pero peor ergonomía para el usuario.
**Consecuencias:** añadir un modelo exige recompilar la librería (una línea de registro);
no se puede extender desde fuera sin recompilar. Coste en ejecución: un `select case` por
llamada, uniforme por material (sin divergencia en GPU), o cero si se monomorfiza.

## D-008 — Materiales declarados + API por lotes residente en GPU
**Estado:** Aceptada (2026-10-28, decisión del usuario)
**Contexto:** el caso de uso principal es un esquema FFT → modelo constitutivo → iFFT con
todos los datos residentes en la GPU (sin transferencias host↔device salvo escalares de
control). Hay pocos usuarios externos de la API OO y se reescribirán.
**Decisión:**
1. **Kernel por material (prioridad):** el usuario declara combinaciones en un registro
   (`MATERIALS`, archivo Python pasado con `-DMUSCLE_MATERIALS=...`): criterio, endurecimiento,
   ley vp, elasticidad y solver. FYPP genera un módulo por material reutilizando un único
   cuerpo de solver; la variación entre materiales es sólo de `use ..., only: yield_f => stress_eq`
   (renombre). Sin `select case` en el kernel; `NSTATE` conocido en compilación.
2. **API por lotes sobre datos ya en GPU:** `muscle_update_<mat>(params, eps, state, sig, status)`
   con arrays en layout SoA (`eps(npts, 6)`); dentro del kernel cada hilo carga a tensores
   locales. La librería no hace `map(tofrom)`: usa `map(present: ...)`/direcciones de
   dispositivo. Multifase = un kernel por material sobre la lista de índices de su fase.
3. **Adaptadores UMAT** (LS-DYNA/ANSYS) generados por material sobre la misma API.
4. **Uniones etiquetadas (D-006/D-007 punto 3):** opcionales y posteriores, sólo si se
   necesita elegir el modelo en tiempo de ejecución.
5. **Sin fachada OO.**
**Condición previa:** firmas unificadas por familia (fluencia, endurecimiento, vp,
elasticidad) antes de la Fase 7.
**Consecuencias:** cambiar la combinación de modelos exige recompilar (los parámetros no);
compilación lineal en nº de materiales; CI con un registro de prueba que cubre cada modelo
al menos una vez. La FFT (cuFFT/rocFFT/oneMKL) queda fuera de MUSCLE, pero la API debe
aceptar punteros de dispositivo.

## D-009 — Marcadores de línea y preprocesado de los `.f90` generados
**Estado:** Aceptada (2026-10-28, Fase 0)
**Contexto:** con `--line-numbering`, los formatos `cpp` y `gfortran5` de fypp producen el
aviso "linemarker ignored due to incorrect nesting" en el paso de preprocesado separado que
CMake usa para Fortran (los marcadores de `#:include` llevan flags de anidamiento).
**Decisión:** fypp siempre con `--line-numbering --line-marker-format=std` (`#line`), en
todas las configuraciones (el comando no depende de la configuración, así los generadores
multi-config comparten el `.f90`; errores y cobertura apuntan al `.fypp`). Los `.f90`
generados se compilan con `Fortran_PREPROCESS ON`. Las plantillas no usan directivas `cpp`:
los interruptores de build (`USE_MKL`, `ENABLE_UDTIO`) llegan a fypp vía
`MUSCLE_FYPP_DEFINES` y se prueban con `#:if defined(...)` (el lint lo comprueba).
**Pendiente:** verificar `#line` con ifx/IntelLLVM (`-fpp`) cuando esté instalado. *(2026-10-30:
ifx 2026.2 acepta los marcadores en todos los módulos generados; falta comprobar que los errores
apunten al `.fypp`.)*

## D-010 — Constantes por kind de la Fase 1
**Estado:** Aceptada (2026-10-28, Fase 1)
**Decisión:**
- `TOL_APPROX_sp = 1e-5` (~84 ulp) y `TOL_APPROX_dp = 1e-12` (histórico) en `muscle_kinds`:
  tolerancia relativa por defecto de `.approx.` y de `check_approx`. `EPS_ABS_<k> = 1e-30`
  (cota inferior de la escala, normal en ambos kinds).
- `ADJ_KAPPA_MAX_sp = ADJ_KAPPA_MAX_dp = 1e3`: medido con 30k matrices aleatorias 3x3/4x4 por
  kind, el adjunto en sp degrada igual que en dp en unidades de κ·ε.
- `QUIET_NAN` (sin argumentos) sólo puede ser genérico de un específico: es el `dp`; las
  plantillas llaman `QUIET_NAN_<k>`.
**Resuelto en D-011:** `is_approx` usa `1.0e-12_dp` (antes literal real32 ≈9.9999998e-13).
**Resuelto (2026-10-30, decisión del usuario):** tolerancia por defecto de `4O2sym`: `1e-8` en dp
(se mantiene) y `1e-4` en sp (10× la base sp), a validar midiendo errores reales de 4O2sym en sp;
se pasa con `tol_default` de `TENSOR_COMMON_PROCEDURES` (D-012).

## D-011 — Patrón de las plantillas de tipos (piloto 2Osym)
*(Enmendada por D-012: código común en macros.)*
**Estado:** Aceptada (2026-10-29, punto de control de la Fase 2 con el usuario; incluye los
cambios de comportamiento listados abajo)
**Decisión:** lo que copiarán los demás tipos (Fase 3) y operadores (Fase 4):
1. **Un módulo por familia** (`muscle_tensor_2osym`) con todas las dims × kinds de esa familia,
   generado desde `layout(familia, dim)`. Lo que es genuinamente distinto por dimensión y
   cuyos bits importan (`square`, `det`) se escribe explícito en `#:if d == ...`.
2. **Operaciones como procedimientos libres con `type`** marcados `@:DEVICE()`; genéricos
   públicos del módulo de familia (`norm`, `det`, `square`, `is_approx`, `inverse`, operadores).
   Los genéricos con nombre de función **no** se re-exportan en `muscle_tensors` (evita chocar
   con variables `norm`/`det` del usuario); los kernels los importan del módulo de familia
   (`use muscle_tensor_2osym, only: norm, det`). **Decisión del usuario (2026-10-29):** esto
   ahora, y al terminar las familias (Fase 3–4) se crea `muscle_tensor_functions`, que reúne las
   funciones libres de todas las familias en un solo `use` (fuera de `muscle_tensors`). Hecho
   2026-10-30: `src/tensors/muscle_tensor_functions.f90`.
3. **Type-bound procedures = envoltorios `*_tbp` con `class`** (Fortran exige `class` en el
   passed-object), sólo para la API histórica; `init` y accesores (`%xx()`) sólo existen así.
   Ningún operador ni kernel llama a un procedimiento con `class` (con `class` gfortran pasa
   un contenedor con puntero a vtab). `tools/check.sh` lo verifica (reglas `oop`, `device`).
4. **Nombres específicos:** operación + tipos + kind al final, una vez:
   `ddot_3D2Osym_3D2Osym_dp`, `mul_3D2Osym_real_dp`; listas con la macro `specifics()`.
5. **Compatibilidad:** `muscle_tensor_3d2osym` / `muscle_tensor_2d2osym` (`.f90`) re-exportan
   `ten_3D2Osym => ten_3D2Osym_dp` y los operadores; se borran en la Fase 4, cuando los
   operadores usen el módulo de familia. `muscle_tensors` (`.fypp`) exporta además `_sp/_dp`.
   (Hecho en la Fase 4, commit dcf55c1: la fachada da los alias dp directamente.)
6. **Salida formateada (UDTIO)** en `muscle_tensor_io` (host-only), generada del layout; opción
   CMake `MUSCLE_UDTIO=AUTO|ON|OFF` y `tools/check.sh` compila el camino `ON`.
**Cambios de comportamiento (verificados bit a bit, 200k entradas/operación, -O0/-O2/-O3):**
- `is_approx`/`.approx.`: umbral `1.0e-12_dp` y `EPS_ABS = 1.0e-30_dp` (antes literales real32);
  sólo cambia el veredicto en el borde exacto del umbral (aprobado por el usuario).
- Inversa 3D: la diagonal ya no pasa por `0.5*(x + x)`: idéntica salvo |x_ii| > huge/2, donde
  antes daba Inf y ahora el valor correcto (como ya hacía la 2D).
- UDTIO 3D modo lista: orden de almacenamiento (`… xy, yz, xz`); antes `… xy, xz, yz`.
- (Fase 3, decisión del usuario 2026-10-29) El umbral `1.0e-12_dp` aplica a todas las familias.
  UDTIO de `ten_3D2O` en modo matriz: las tres filas en `ES` (el original imprimía las filas 2–3
  con `F`, errata). `4O3sym` gana salida formateada (no tenía): lista = componentes en orden de
  almacenamiento, matriz = matriz de Voigt simétrica \(C_{IJ}\) (6x6 / 4x4, componentes
  tensoriales); igual que el resto, sólo con `ENABLE_UDTIO` (`MUSCLE_UDTIO`) por el bug de gfortran.
- (Fase 3, tanda 2, 2026-10-30) 4O2sym: umbral por defecto `1.0e-8_dp` (antes literal real32 `1E-8`
  ≈9.99999994e-9) y `EPS_ABS_dp` (antes `1e-30` real32): sólo cambia el veredicto en el borde, como
  el `1e-12` de arriba. 4O2sym gana salida formateada (lista + matriz 6x6/4x4, decisión del usuario).
  4O conserva su salida histórica en cualquier modo (`IO_MATRIX_ONLY`); única diferencia: el 3D
  ahora rellena `iomsg` al fallar (el 2D ya lo hacía). API nueva: genéricos libres `norm`,
  `is_approx`, `convert_3sym` (4O2sym) y `norm`, `is_approx`, `convert_2sym` (4O), más las
  variantes `sp`. Explícito por dimensión (bits): orden de suma de `norm` en 4O2sym (por filas en
  3D, por columnas en 2D, como los originales).
Resto de operaciones: idénticas bit a bit. Benchmark `closest_point` Release: sin cambio.

**Evidencia (2026-11-11, experimento del usuario):** ¿sirve un procedimiento ligado al tipo en una región `target` si
se llama desde una variable `type(t)` del mismo tipo? Programa aislado con función libre, método del mismo tipo (también
`non_overridable`), método heredado y sobrescrito de un tipo extendido, y despacho dinámico real (`class(t)` decidido en
ejecución), dentro de `!$omp target teams distribute parallel do`:

| Caso | ifx → Intel UHD 770, sp | ifx, dp emulado | gfortran 15.2 → NVPTX / AMDGCN (sólo build: no hay GPU) |
|---|---|---|---|
| Función libre con `type(t)` | funciona | funciona | enlaza |
| Método del mismo tipo (y `non_overridable`) | funciona | funciona (programa mínimo) | **no enlaza con `-O0`**; enlaza con `-O2` |
| Heredado / sobrescrito en tipo extendido | funciona | el JIT no termina (> 240 s) | no enlaza con `-O0` |
| Despacho dinámico | funciona | el JIT no termina | no enlaza con `-O0` |

Causa con gfortran: Fortran exige `class(t)` en el argumento del método aunque se llame con `type(t)`, y gfortran le
pasa un contenedor con puntero a la tabla de métodos: «variable '__vtab_min_m_T' has been referenced in offloaded code
but hasn't been marked to be included in the offloaded code». La tabla la genera el compilador y no se puede marcar
`declare target`; con `-O2` la referencia desaparece por optimización (dependería del nivel de optimización). ifx sí lo
soporta. **La regla se mantiene** (los kernels sólo llaman funciones libres; los `%metodo()` son API histórica de host)
mientras el proyecto apunte a NVIDIA y AMD con gfortran; revisarla si gfortran lo corrige o si se usara sólo ifx.
A raíz de esto, `is_identity` (antes privada de `muscle_kin_finite`) pasó a ser función libre de `2Osym` y `2O`
(exportada por `muscle_tensor_functions`, comparación exacta) con su envoltorio `%is_identity()`.

## D-012 — Enmienda a D-011: código común de las plantillas de tipo
**Estado:** Aceptada (2026-10-30; refactor pedido por el usuario como prioridad antes de la
segunda tanda de la Fase 3)
**Contexto:** tras migrar 2Osym, 2O y 4O3sym, ~145 líneas eran idénticas carácter a carácter
entre las tres plantillas, y la inversa de 4O3sym necesitaba `#:if d` sólo porque `M66INV` no
devolvía el número de condición y `M44INV` sí.
**Decisión:**
1. `include/muscle_tensor_common.fypp` con macros que generan lo común a toda familia:
   `TENSOR_COMMON_INTERFACES(F)`, `TENSOR_COMMON_PROCEDURES(d, F, k)` (assign, is_approx,
   approx, sum, sub, subU, mul ×2, div, inv), `TENSOR_COMMON_TBP(...)` (init_tbp, envoltorios
   `%op()`, is_approx_tbp) y `LOAD_MATRIX` (matriz desde `Layout.matrix()`). Se llaman con `$:`.
   En cada plantilla queda sólo lo específico de la familia (init2/set, fórmulas, escritura de la
   inversa).
2. Todas las inversas de `muscle_math_inverses` aceptan `KAPPA` opcional (`M55INV`, `M66INV`,
   `M99INV` lo calculan con `COND_INF` tras el LU, sólo si se pide); las plantillas pasan `cond`
   directamente y llaman siempre a la rutina específica por kind (`M66INV_${k}$`).
3. Los tests comparten helpers en `tests/common/`: `muscle_test_helpers` (generador, referencias
   dp, `record`), `check_scaled` en `muscle_test_utils` y macros de test por familia en
   `muscle_test_tensor_cases.fypp` (ver `.claude/rules/testing.md`). Inventario dinámico de checks:
   ninguno se pierde; mutaciones detectadas 74/78 (antes 73/78).
**Verificación:** resultados bit a bit idénticos a `eaee906` (20 000 casos × 12 tipos, todas las
operaciones, librería -O0/-O3 × programa -O0/-O2/-O3); mismos símbolos exportados.
**Consecuencia:** una familia nueva (4O2sym, 4O) sólo escribe lo específico; cambiar un criterio
común (p.ej. `is_approx`) se hace en un solo lugar.

## D-013 — Umbral adjunta/LU en precisión simple y tolerancia sp de 4O2sym
**Estado:** Aceptada (2026-10-30, decisión del usuario tras la tanda 2 de la Fase 3)
**Contexto:** el `test-builder` de 4O2sym midió que con `ADJ_KAPPA_MAX_sp = 1e3` la inversa por
adjunta en sp pierde precisión frente a LU para κ∞ entre 1e2 y 1e3 (contradiciendo la medición de
la Fase 1, D-010). Medición de confirmación (matrices aleatorias U·diag·Vᵀ, 2000 por banda de κ):
| κ∞ | 4x4 adjunta p99 / LU p99 | 3x3 adjunta p99 / LU p99 |
|---|---|---|
| 30–100 | 1.3e-6 / 6.6e-7 | 9.1e-7 / 6.8e-7 |
| 100–300 | 5.9e-6 / 1.7e-6 (3.5×) | 3.7e-6 / 1.7e-6 (2.1×) |
| 300–1000 | 3.5e-5 / 4.6e-6 (7.6×) | 1.8e-5 / 5.9e-6 (3.0×) |
**Decisión:**
1. `ADJ_KAPPA_MAX_sp = 1e2` (dp sigue en 1e3: sin cambios en dp). Afecta sólo a `M33INV_sp`/`M44INV_sp`.
2. La tolerancia por defecto de 4O2sym en sp se mantiene en 1e-4: adecuada hasta κ∞ ≈ 1e3
   (elasticidad hasta ν = 0.499, tangentes con endurecimiento moderado); para inversas peor
   condicionadas se documenta pasar `tol ≈ 2·κ∞·ε` explícito.

## D-014 — Identidades (`iden_*`) en simple y doble precisión
**Estado:** Aceptada (2026-10-30, decisiones del usuario antes de migrar las identidades)
**Contexto:** las identidades no dependen de la dimensión (valen para 2D y 3D) y las no escaladas
(`iden_2O`, `iden_4O3T`, `iden_4O4T`) no guardan ningún número.
**Decisión:**
1. Las identidades sin escalar siguen siendo **un único tipo sin kind** (`iden_2O`, ...): sirven
   para sp y dp; el kind del resultado lo fija el escalar (`2.0_sp * I` → `iden_2OS_sp`).
2. Las escaladas tienen una versión por kind (`iden_2OS_sp`, `iden_2OS_dp`, ...) y el nombre
   histórico es el dp (`iden_2OS => iden_2OS_dp`, D-003).
3. Las operaciones entre identidades sin escalar, que no tienen escalar que fije el kind:
   `I + I`, `I - I`, `-I` devuelven **`iden_2OS_dp`** (doble precisión), y `I * I` devuelve
   `iden_2O`, como hoy.
4. Salida formateada: la de `iden_2O`/`iden_2OS` se conserva (movida a `muscle_tensor_io`) y las de
   4º orden (`iden_4O3T(S)`, `iden_4O4T(S)`) **ganan** salida formateada; todo sólo con `ENABLE_UDTIO`.
5. (Tras migrar, decisión del usuario) Tolerancia de `.approx.` de las identidades escaladas:
   `1.0e-10_dp` (antes el literal real32 `1e-10` ≈ 1.00000001335e-10: sólo cambia el veredicto en el
   borde, como D-011) y `1.0e-5_sp`. El mensaje de error de la salida de `iden_2O`/`iden_2OS`, que
   decía `print_ten_3D2Osym` (errata), nombra ahora su propia rutina.

## D-015 — Organización y kinds de los operadores entre tipos (Fase 4)
**Estado:** Aceptada (2026-10-30, decisión del usuario al planificar la Fase 4)
**Contexto:** 6 módulos `.f90` (`muscle_tensor_ops_*`) con 222 procedimientos dp (119 en 3D, 103
en 2D: todos tienen su par en la otra dimensión salvo las 16 sumas/restas de dos identidades de 4º
orden, que devuelven siempre `ten_3D4O3sym`); casi todos combinan tensores con identidades.
**Decisión:**
1. Se conservan **6 módulos, uno por operación**, convertidos a `.fypp` y generados desde
   una tabla de pares por operador en `muscle_meta.py` × `DIMS` × `KINDS`. Alternativas descartadas:
   un único módulo (diff enorme, sin tandas paralelas) y operadores dentro del módulo de la familia
   del primer operando (ciclos de dependencias entre familias).
2. **Sin mezcla de kinds:** cada operador existe en sp y en dp con ambos operandos del mismo kind;
   `sp ⊗ dp` no compila (el usuario convierte explícitamente). Las identidades sin escalar se
   combinan con ambos kinds y las escaladas sólo con su kind (D-014).
3. El conjunto de interfaces dp generado debe contener al actual (script de inventario en `tools/`,
   comprobado por `check.sh`).
4. (Paso 2, decisión del usuario 2026-10-30) Se amplía D-014.3 a los módulos de operadores: **toda
   operación entre identidades sin escalar devuelve dp**, también con resultado tensorial
   (`iden_4O3T ± iden_4O4T` y espejos → `ten_3D4O3sym_dp`, 4 procedimientos sólo en dp; una versión
   sp tendría los mismos argumentos y haría ambiguo el genérico). En sp: `1.0_sp*iden_4O3T() + iden_4O4T()`.
   En la tabla: `OP_KINDLESS_KIND` de `muscle_meta.py`.
5. Las 16 sumas/restas de dos identidades de 4º orden existen sólo en 3D (como en el original).
6. (Paso 3, decisión del usuario 2026-10-30) Cambio de comportamiento aprobado: la diagonal de
   `a .tdotsym. b` en 3D se escribe `a_I*b_I` como en 2D, en vez de `0.5*(a_I*b_I + b_I*a_I)`.
   Mismos bits salvo cuando `a_I*b_I` > huge/2, donde la forma antigua daba `Inf`.
   Las plantillas son reescrituras completas (p.ej. 1886 → ~80 líneas): git no las reconoce como
   renombre y `git log --follow` se corta en el `.f90` original (su historia sigue en git).
7. Diferencia 2D/3D heredada que **se conserva** (cambiarla movería ~1 ulp en datos normales):
   en `.transform.` de `2Osym` la componente 33 en 2D es `(F33*F33)*S33`, como el original, y no
   `(F33*S33)*F33` como daría la regla general (`op_congruence`, regla de producto único).
8. (Tanda B, decisión del usuario 2026-10-30) **Cambio de comportamiento aprobado:** `.ddot.` se
   escribe con una sola regla en todos los kernels y dimensiones: suma sobre el índice Voigt
   contraído en orden ascendente, con el peso de corte sobre el factor del operando derecho,
   `x*(2*y)` (`op_weighted`). Antes: `4O3sym:2Osym` ponderaba `A*(2b)` en 3D y `2*(A*b)` en 2D;
   `4O2sym:2Osym` usaba `(2a)*b`; `4O3sym:4O3sym` era en 3D una suma en el orden de sympy y en 2D un
   `matmul` 4x4. Sólo `4O3sym:4O3sym` cambia bits con datos normales (≤ 1 ulp de la escala Σ|términos|; por
   componente hasta 32 ulp / 4e-15 relativo donde hay cancelación, dentro del 1e-14 de `testing.md`;
   golden actualizado sólo en esos 6 registros); las demás formas dan los mismos bits salvo
   subnormales o desborde de `2*x`. La validación pasa a ser contra la contracción de los tensores
   completos (3x3 / 3x3x3x3), con tolerancia por kind.
9. FMA: con `-mfma`/`-march=native` el compilador puede fusionar productos y sumas y cambiar el
   redondeo frente al golden. El proyecto no usa esos flags; la equivalencia bit a bit con el golden
   se garantiza sin FMA.

## D-016 — Math sin closures: alcance de la Fase 5
**Estado:** Aceptada (2026-10-31, decisiones del usuario al planificar la Fase 5)
**Contexto:** las diferencias finitas de `muscle_math_derivatives` reciben la función como
`procedure(...)`; las llaman 8 sitios de modelos y solvers con funciones internas que capturan
`self` (trampolines, no aptas para GPU), que sólo desaparecen al migrar los modelos (Fase 6). El
paso relativo dp es 1e-7 (escalar, 2O), 1e-6 (escalar de 2Osym) y 1e-4 (segunda derivada y
2Osym→2Osym), mínimo 1e-40 (1e-35 en 2Osym→2Osym).
**Decisión:**
1. **Convivencia hasta la Fase 6:** la Fase 5 añade macros FYPP de diferencias finitas que llaman
   directamente a la función concreta (sp/dp, `pure`, `@:DEVICE()`); la API con `procedure(...)`
   se conserva (sp/dp), marcada host-only, y se borra en la Fase 6 cuando los modelos usen las
   macros.
2. **Pasos por kind:** dp conserva exactamente los valores actuales (mismos bits); sp usa pasos
   derivados de ε_sp: ~ε^(1/3) (≈5e-3) en la primera derivada central, ~ε^(1/4) (≈2e-2) en la
   segunda, mínimo absoluto 1e-30. Constantes por kind en `muscle_kinds` (como D-010).
3. **Se conservan** las `derivative_escalar_tensor3x3*` de `muscle_math_operations` (migradas a
   sp/dp como API host-only) y la copia sin registrar `tests/test_muscle_math_operations.f90`.
4. (Decisión del usuario 2026-10-31, tras el paso 0) **Bug corregido en el paso 2:** con argumento
   nulo el paso mínimo absoluto 1e-40 (1e-35 en 2Osym→2Osym) hace que f(±paso) se redondee a f(0)
   y la derivada numérica sale 0 sin serlo (p.ej. `derivative2O` en el tensor cero). Regla elegida
   por el usuario (paso 2, 2026-10-31): `h = rel` si el argumento es exactamente 0, si no el paso
   histórico `h = max(|x|·rel, h_min)` (|x| = `abs` o la norma del tensor; `h_min` el mínimo
   histórico por kind: 1e-40 dp, 1e-35 en el jacobiano, 1e-30 sp), de modo que sólo cambia x = 0.
   **Revisión (cierre de la Fase 5, decisión del usuario 2026-10-31):** en x = 0 el paso de las
   primeras derivadas (`FD_DERIV`, `FD_GRAD`, `FD_JAC`) es `FD_STEP_ZERO_<k> = sqrt(ε_k)` (≈1.5e-8 dp,
   ≈3.5e-4 sp) en vez del paso relativo. Motivo: el paso en cero debe ser grande para que el
   redondeo no lo borre (error ~ε|f(0)|/h) y chico para no salir de la zona suave (en los return
   mappings, `derivative(wrapper, strain)` con `strain = 0` en el primer incremento: con el paso
   relativo 1e-4 la perturbación podía cruzar la fluencia y dar un tangente elastoplástico en cero;
   el histórico 1e-35 funcionaba sólo por coincidencia, porque σ(0) = 0 y es lineal). √ε da ~1e-8 de
   precisión y una perturbación de deformación despreciable. La segunda derivada (`FD_HESS`)
   conserva el paso relativo ~ε^(1/4) en x = 0 (con √ε el redondeo ε/h² sería de orden 1). En sp
   √ε ≈ 3.5e-4 puede acercarse a una deformación de fluencia: en la Fase 6 los return mappings deben
   pasar un paso explícito ligado al material (p.ej. una fracción de la deformación de fluencia) o
   usar el tangente consistente analítico, dejando la diferencia finita como respaldo. Se descartó `h = max(|x|·rel, rel)` porque
   cambiaba todo |x| < 1, incluidos los tangentes de los modelos (deformaciones ~1e-3, hasta 6e-3
   relativo), y un suelo `rel·√ε` porque en x = 0 deja un paso demasiado chico (error ~10%).
   **Cambio de comportamiento aprobado**, sólo para argumentos nulos; los registros golden
   afectados se actualizan documentados. Efecto en los modelos: las derivadas numéricas por defecto
   de `muscle_vp_base` (`ep`, `epd`) y `muscle_vp_overstress_base` (`ep`), evaluadas típicamente
   con `ep = 0`/`epd = 0` en el primer incremento plástico, pasan de 0 a la derivada real: cambia la
   respuesta de los modelos viscoplásticos en ese primer paso (ningún test de modelo lo fijaba).
   Limitación: un x no nulo pero diminuto conserva el paso histórico (y su imprecisión).
   Se extiende a `derivative_escalar_tensor3x3sym` (`muscle_math_operations_host`), cuyo paso es
   relativo a |f(x)|: pasa a ser relativo a |x| con el mínimo por kind (decisión del usuario,
   2026-10-31).
5. (Paso 1, decisión del usuario 2026-10-31) `eigenvals`: el umbral histórico `EPS=1e-10` era un
   literal real32 (1.0000000134e-10); pasa a `1.0e-10_dp` (**cambio aprobado**, sólo la franja
   (1e-10, 1.0000000134e-10]). En sp el umbral es 1e-9: es una guarda contra 0/0, no un compromiso
   de precisión (medido: error ≤ 2e-7 frente a ~3e-5 con 4e-5).
6. (Paso 1, decisión del usuario 2026-10-31) Derivadas espectrales: los umbrales `EPS_TOL` y
   `EPS_TOL_SECOND` (antes absolutos, en unidades de tensión²) pasan a ser **relativos a
   max(λ1², λ2², λ3²) en sp y dp** (no depende del orden de los autovalores), con los mismos valores por kind. **Cambio de comportamiento
   aprobado en dp** (igual para tensiones de orden 1; distinto para escalas lejanas de 1, donde el
   umbral absoluto era incorrecto: en sp daba |H| hasta 32 con autovalores repetidos y escala
   1000). `REPEATED_TOL` ya era relativo (1e-9 dp / 1e-4 sp). Ningún golden cambió (casos de
   orden 1).
7. Pasos de diferencias finitas por kind en `muscle_kinds` (D-016.2), compartidos por la API
   host-only y las macros del paso 2; los umbrales propios de un algoritmo (eigenvals, espectral)
   quedan locales en su plantilla. Paso adelantado sp ≈ √ε_sp (3e-4), con el mismo criterio que
   el 1e-8 ≈ √ε_dp histórico.
8. Patrón host-only: las rutinas que reciben `procedure(...)` viven en módulos `*_host`
   (`muscle_math_operations_host`, `muscle_math_derivatives_host`), re-exportados por el módulo de
   kernel, como `muscle_tensor_io`. Se borran en la Fase 6 (D-016.1).
9. Bits macro ↔ API host: en dp las macros de diferencias finitas dan los mismos bits que la API
   host sólo con las mismas opciones de compilación; si el compilador integra (inline) y vectoriza
   la función derivada, cambian los últimos bits (medido 4e-11 relativo con gfortran -O3). En la
   Fase 6, al pasar los modelos de `derivative(wrapper, …)` a las macros, la equivalencia se
   comprueba con tolerancia (no bit a bit).
10. (Paso 3, 2026-10-31; recomendación del revisor) Solver lineal
    pequeño `muscle_math_lu_solve`: `LU<n>SOLVE` (n = 2..9; el 1x1 es una división dentro de `lu_solve`; sp/dp, `pure` + `@:DEVICE()`) con la LU de
    `muscle_math_lu_inverse` compartida con `M<n>INV` (mismos bits verificados). En kernels se llama
    a `LU<n>SOLVE` con n fijo por la plantilla; el frente `lu_solve` (`select case` por tamaño) es
    una comodidad para host y tests (en dispositivo la rama 9x9 fijaría los registros, D-006).
    `KAPPA` opcional es κ∞ exacto (misma definición que en `M<n>INV`, sólo se paga si se pide); un
    estimador O(n²) sólo si el perfil de la Fase 7 lo pide. Contrato: pivote nulo o NaN → `OK_FLAG =
    .false.`, `x = 0`, `KAPPA = huge`; un tamaño inválido se reporta igual que una matriz singular
    (quien llama puede comprobar `size(A,1) <= LU_SOLVE_MAX_N`). `muscle_math_linsolvers` (Gauss,
    Gauss-Seidel, `mklsolve` con `sgesv`/`dgesv`) queda host-only en sp/dp.
11. (Paso 4, 2026-10-31) `spectral_decomposition` escala las filas de T − λI por una potencia de 2
    antes de los productos vectoriales (corrige en sp la pérdida de normalización con cortes
    ~1e-11). El escalado es exacto: dp da los mismos bits salvo cuando esas filas rozan el
    underflow/overflow (≲ 1e-154 o ≳ 1e154), donde el resultado pasa a ser correcto (cambio de
    comportamiento dp sólo en ese borde). (Nota 2026-11-12: desde D-032, K escala el desviador, no las
    filas.)

## D-017 — Fase 6: modelos con macro de contrato y adaptadores OO hasta la Fase 7
**Estado:** Aceptada (2026-10-31, decisiones del usuario al planificar la Fase 6)
**Contexto:** ~7800 líneas de modelos y solvers con despacho dinámico (`class(Base_*)`); los
solvers de return mapping se reescriben en la Fase 7 (materiales declarados, D-008).
**Decisión:**
1. **Adaptadores OO hasta la Fase 7:** cada modelo gana un `type` concreto de parámetros y
   funciones libres `pure` sp/dp (aptas para GPU). Los `Base_*` y las extensiones OO actuales
   quedan como adaptadores finos host-only que llaman a esas funciones, para que solvers, tests y
   UMAT sigan funcionando sin cambios; se borran en la Fase 7 al reescribir los solvers.
2. **Macro de contrato (D-007.1 confirmada):** un modelo se escribe en un archivo corto con
   nombres cortos (`real(wp)`, `type(sym2)`, ...); la macro de su familia (p.ej.
   `#:call hardening_model(...)`) instancia el cuerpo para cada kind (y dimensión cuando aplique),
   inyecta los `use`, los nombres `_<k>` y `@:DEVICE()`, y genera con las macros FD de la Fase 5
   las derivadas que el modelo no dé. Piloto: hardening bilinear (como el piloto de la Fase 2).
3. Golden dp de cada familia capturado contra el código OO original antes de tocarla
   (testing.md); el comportamiento dp no cambia salvo decisión registrada.
4. (Decisiones del usuario 2026-10-31, tras el golden de hardening) Se conservan tal cual, con mismos
   bits en dp, y se documentan: (a) Ludwik: `dstress_dep`/`ddstress_ddep` usan `ep + 1e-10` (no
   `ep`) como regularización en ep = 0 para n < 1 (no son las derivadas exactas de `stress`; en sp
   el offset se escala por kind); (b) Swift con `e0 = 0` (Hollomon): en ep = 0 las derivadas son
   ±Inf; se documenta que hace falta `e0 > 0` (validación opcional en el constructor del adaptador
   OO host).
   Implementación (familia hardening, 2026-10-31; decisión del usuario tras la revisión): Ludwik
   usa `LUDWIK_EP_OFFSET_<k>` = 1e-10 en **ambos** kinds (misma regularización en ep = 0, de modo
   que en la primera fluencia sp y dp se comportan igual y se pueden comparar en los solvers; la
   perturbación relativa de las derivadas es (1-n)·δ/ep: en sp ≤ ~1e-6 en ep = 1e-4 y nula por
   redondeo desde ep ≈ 1.7e-3). Se descartó 1e-14 en sp (derivadas exactas desde ep ≈ 2.4e-7, pero
   un módulo de endurecimiento en ep = 0 ~1600 veces mayor que en dp). Swift no valida `e0 > 0`
   porque tests existentes construyen con `e0 = 0` (sólo se documenta: con e0 = 0 las derivadas en
   ep = 0 son infinitas o NaN).
5. (Piloto hardening bilinear, decisiones del usuario 2026-10-31)
   (a) **Escala por familia en las diferencias finitas:** las macros FD aceptan una escala opcional,
       `h = rel·max(|x|, escala)`; hardening usa escala 1 (escala natural de la deformación
       plástica), de modo que con ep chico pero no nulo el paso no colapsa (σ_y(0) ≠ 0). Los demás
       usos no cambian. Se recomienda dar derivadas analíticas.
   (b) **Genéricos con prefijo de familia:** `hard_stress`, `hard_dstress_dep`, `hard_ddstress_ddep`
       (y análogos `vp_*`, `yield_*`, ...), iguales en todos los modelos de la familia (contrato
       D-007), sin choques entre familias ni con variables locales. Viven en el módulo de familia
       (`muscle_hard_laws`), generado desde el registro `HARDENING_MODELS` de `muscle_meta`.
   (c) **Adaptador OO en el mismo módulo del modelo**, generado por la macro con `oo_type` (un
       `*_host` aparte formaba un ciclo y los usuarios esperan el tipo OO en el módulo histórico):
       mismo patrón que los `*_tbp` de los tensores (D-011); se retira en la Fase 7 quitando
       `oo_type` (desviación aceptada de la letra de D-017.1).
   (d) `FD_DERIV2` (segunda derivada escalar centrada, paso `FD_STEP_CEN2`) se añade a
       `include/muscle_fd.fypp`.
6. (Revisión del piloto, decisiones del usuario 2026-10-31) **Enmienda a D-017.2 / D-007.1:**
   (a) **Patrón invertido:** el autor escribe su plantilla `.fypp` con `#:for k` (y `#:for d`
       cuando aplique) y sus funciones como Fortran normal: errores y cobertura apuntan a sus líneas
       (D-009), puede hacer `use` y escribir código dependiente del kind o la dimensión. Macros
       chicas, comunes a todas las familias, generan sólo lo repetitivo: el `type` de parámetros, la
       cabecera de cada función (`pure`, `@:DEVICE()`, declaraciones), las derivadas por diferencias
       finitas que el autor no dé, el adaptador OO (macro aparte, retirada en la Fase 7) y el registro
       en el módulo de familia. Motivo: con el cuerpo como texto dentro de una macro grande no se
       puede depender de `k`/`d`, no se puede hacer `use` y todas las líneas colapsan a la del
       `#:call`. La ley bilinear se reescribe con el patrón invertido.
   (b) **Composición desacoplada:** las funciones de una familia que hoy contiene otra ley (p.ej.
       las viscoplásticas, que contienen una de endurecimiento) reciben los valores ya evaluados
       (σ_y y sus derivadas), no la ley; quien compone (material declarado, Fase 7) llama primero a
       `hard_*` y luego a `vp_*`. Sin combinatoria, apto para GPU; los adaptadores OO componen en host
       como hoy.
7. (Golden viscoplásticas, decisión del usuario 2026-10-31) Los bugs encontrados en las leyes
   viscoplásticas **sólo se documentan** (FORD del modelo y esta lista) y se migran con los mismos
   bits en dp; no se corrigen en la Fase 6: RK `dstress_dep` omite ∂σ/∂ep; RK y JCfull dan
   `flow_stress = 0` con epd ≤ 0; MRK y RK dan Inf en epd = epdmax y NaN por encima; MRK salta en
   epd = 0⁺ y su `dstress_dep` es 0 con epd ≤ 1e-12 y finita donde σ es NaN/Inf; con
   F = 1 − χ1·L ≤ 0 la flujo de MRK es plana (σ_u) pero sus derivadas usan F = 1e-10 y, con χ2 > 1,
   son grandes; NNL `dstress_dep`
   ignora los recortes de `flow_stress` (suelos de ep distintos, X ≤ 0.999, `eps_log`); VA deriva en
   el suelo donde σ es plana; Perzyna divide por cero si Y(ep) = 0; literales sin kind en MRK/RK.
   Además (bug de gfortran 15, no de MUSCLE): reasignar un `class(...), allocatable` ya asignado con
   otro tipo dinámico (p.ej. `law%hard_law = otra_ley`) no realoca y desborda el heap; hay que
   desasignar antes. (Nota: los literales de RK escritos con kind dan los mismos bits; ese punto
   desaparece al migrar.)
8. (Piloto viscoplástico JC/RK/Perzyna, decisiones del usuario 2026-10-31) Dos familias en
   `muscle_meta`: `vp` (tensión de flujo σ(ep, epd); módulo `muscle_vp_laws`) y `os` (overstress,
   g(σ_eq, ep); módulo `muscle_vp_overstress_laws`; Perzyna y Chaboche). API de kernel de `vp`:
   parciales `vp_dflow_dep` (∂σ/∂ep) y `vp_dflow_depd` (∂σ/∂epd) **y** la total histórica
   `vp_dstress_dep(…, dt)`, generada por defecto como parcial_ep + parcial_epd/dt (RK escribe la
   suya: su `dt` va dentro de cada término y sacarlo cambiaría bits). Las leyes compuestas reciben
   `sy`/`dsy_dep` (D-017.6b); `MODEL_FD` no admite diferencias en ep con σ_y congelado. RK: las
   parciales son coherentes con su bug documentado (`vp_dflow_dep_rk = 0`). JCfull: ∂σ/∂ep analítica
   nueva en la API de kernel (con test contra diferencias finitas) y su adaptador OO sigue heredando
   el numérico de la base (golden sin cambios). Las familias `vp`/`os` mantienen el paso FD actual
   (mismos bits que el numérico de las bases).
9. (Golden yield, decisión del usuario 2026-10-31) **Se corrige** el NaN (0/0) de las derivadas
   analíticas de los criterios con tensión equivalente nula (tensor cero e hidrostáticos: von Mises,
   Hill48, CPB06, CPB06ex2, Yld2004): en f = 0 el gradiente y la hessiana devuelven 0 (0 pertenece
   al subdiferencial de un criterio convexo y homogéneo en ese punto; Drucker-Prager ya tenía
   guardas y no cambia). **Cambio de comportamiento aprobado en dp** sólo en esos estados (golden
   de esos registros actualizado). El resto de anomalías del golden de yield **sólo se documenta**
   (como D-017.7): umbral absoluto 1e-10 de `eigenvals` que descarta cortantes chicos sobre una
   hidrostática (CPB06, CPB06ex2, Yld2004); Yld2004 no exactamente insensible a la presión
   (`linear_transform` sobre σ, f ≈ 7e-17·p con AA2090); Drucker-Prager sin `init` (K = 0) da NaN/Inf
   en las derivadas numéricas; en el ápice de Drucker-Prager la derivada analítica (subgradiente) y
   la numérica difieren; la hessiana de CPB06 con k = ±1 y un autovalor nulo es la elegida en
   `curvature`.
10. (Piloto yield von Mises/Hill48/CPB06, 2026-10-31) Diseño de la familia `yield`: argumento
    `sigma` (`2Osym`) y funciones `stress_eq`, `dstress_eq` (→ `2Osym`), `ddstress_eq` (→ `4O3sym`);
    nombres específicos con dimensión `<fam>_<func>_<model>_<d>_<k>` (amplía D-011.4; los modelos sin
    dimensión no la llevan); `MODEL_FD` sobre `2Osym` con `FD_GRAD`/`FD_HESS` (resuelve el límite (c)
    anotado en PLAN); adaptadores OO en 3D (`OO_DIM`); por ahora sólo 3D (los cuerpos usan posiciones
    3D; 2D cuando haya test 2D vs 3D). **Decisiones del usuario:**
    (a) **Evaluación conjunta** en todos los criterios: subrutina
        `yield_stress_eq_and_derivatives(prm, sigma, f, df [, ddf])` con **`ddf` opcional** (si no se
        pide no se calcula: sirve a closest point y a cutting plane sin repetir la descomposición) y
        bandera constante `yield_joint_evaluation(prm)`; por defecto llama a las tres funciones;
        Yld2004 escribe la suya.
    (b) **Parámetros derivados privados:** la maquinaria admite marcar parámetros como privados; los
        derivados (B y L de CPB06, análogos en CPB06ex2, Yld2004 y Drucker-Prager) son privados y el
        inicializador `yield_<m>_init` es la única forma de construir el tipo.
    (c) **CPB06 reescalado por max|λ|** (df homogénea de grado 0, d²f de grado −1): correcto en todo el
        rango en sp (antes Φ subnormal → error o Inf con |σ| ≲ 3e-6) y dp. Implementado con la potencia
        de 2 de max|λ| (escalado exacto): ningún registro analítico dp cambia; sólo 1 registro de la
        derivada numérica de la base (1.8e-12 relativo, 1 ulp de f amplificado por la diferencia).
    (d) **Hessiana de von Mises en forma estable** (9/(4·seq))·(n⊗n), n = dev/seq: sin
        under/overflow (antes seq³ desbordaba en sp con seq ≲ 1.5e-13 o ≳ 7e12). **Cambio de
        comportamiento aprobado en dp**: 17 registros, ≤ 1.7e-16 relativo (redondeo).
11. (Drucker-Prager, CPB06ex2, Yld2004; decisiones del usuario 2026-11-01)
    (a) **Reescalado siempre (D-017.10c) aceptando el redondeo:** con exponentes no enteros o
        grandes (CPB06ex2 a = 5, 12; Yld2004 a = 6) `pow` no es exactamente homogéneo y f cambia en
        ≤ 12 ulp (≤ 2.4e-15 relativo): 77 registros analíticos dp (máx. 1.02e-14 relativo en una hessiana, justo por encima
        del 1e-14 de testing.md: excepción aprobada) y 89 de la derivada numérica de la base (máx.
        9.3e-9). CPB06 (a = 2, 8) no cambia.
    (b) **Umbral sp del eje hidrostático de Drucker-Prager:** 1e-5 relativo (≈ 100·ε_sp; cubre el
        ruido de redondeo q/|σ| ~ 8e-8 y limita la zona donde el gradiente pierde la parte desviadora,
        relevante en estados muy confinados; revisado por el usuario tras la revisión: la regla
        tol^0.442 daba 1e-4); dp sin cambios (1e-9); absoluto sp 1e-30 (equivale a |σ| = 0). Son
        constantes propias del modelo, locales en su módulo.
    (c) El código común de los criterios anisótropos (transformación lineal, pull-backs, términos
        principales y curvatura de la familia CPB) vive en `muscle_anisotropic_yield` (no es un modelo
        registrado).
    (d) **Hessiana de Drucker-Prager en forma estable** (cocientes en vez de q³, q⁴; como D-017.10d):
        antes, en sp, q⁴ desbordaba con q ≳ 4.3e9 y la hessiana salía finita pero incorrecta (10–50%).
        **Cambio de comportamiento aprobado en dp** (redondeo): 94 registros, ≤ 6.7e-16 relativo.
        `stress_eq` y el gradiente conservan r³/q³ (mismos bits): en sp valen para 2e-13 ≲ q ≲ 7e12
        (documentado).
12. (Golden elasticity, decisiones del usuario 2026-11-01)
    (a) **Un solo camino para el desviador:** la plantilla 2D/3D usa la fórmula analítica del
        desviador (la del `stress_dev` 3D actual) en todos los caminos (2D y `_transient`). **Cambio
        de comportamiento aprobado en dp** sólo por redondeo (~1e-16 relativo a la tensión) en los
        registros 2D y `_transient` del desviador; los 3D no transitorios no cambian. Ampliación
        (usuario, 2026-11-01): con ν = 0.5 esos caminos pasan de NaN a su valor finito y correcto (el
        desviador analítico sólo usa μ), igual que el `stress_dev` 3D desde siempre (17 registros);
        la tensión completa sigue siendo NaN (D-017.12b).
    (b) **ν = 0.5 sólo se documenta:** `set_parameters` divide por cero (λ = +Inf; tensión normal NaN
        y tangente con Inf); la elasticidad lineal exige ν < 0.5 (incompresible exacto no soportado).
        También se documentan: las rutas `_transient` y `stress_and_update_3D` ignoran `dstrain`,
        `dtime` e `internal`; no existen `stress_and_update_2D` ni `dstress_dstrain_dev` 2D en la API OO
        (en la API de kernel sí existe el desviador del tangente 2D).
    (c) **Un solo camino también para el desviador del tangente:** la variante `_transient` deja de
        heredar el proyector de la base y usa la fórmula analítica (4/3·μ, −2/3·μ). Cambio aprobado en
        dp en 61 registros: 43 por redondeo (≤ 60 ulp, ≤ 7.6e-15 relativo); 9 con ν = 0.4999 que
        cambian 3.4e-13 relativo (excepción aprobada al 1e-14 de testing.md: el proyector perdía
        precisión por cancelación, error ∝ (λ/μ)·ε; el analítico es el correcto); 9 con ν = 0.5 de NaN
        a finito (como la ampliación de (a)).
    (d) Diseño de la familia `elas` (piloto, 2026-11-01): prefijo `elas`, argumento `strain` (`2Osym`),
        **dimensiones 3D y 2D en la API de kernel** (una plantilla; embebido 2D/3D testeado); funciones
        `stress`, `stress_dev`, `dstress_dstrain`, `dstress_dstrain_dev` (esta última también en 2D,
        nueva en kernel); módulo de familia `muscle_elasticity_laws`; las variantes `_transient` y
        `stress_and_update_3D` quedan fuera de la API de kernel (sólo en el adaptador OO: hoy ignoran
        `dstrain`/`dtime`/`internal`; un modelo con estado irá a la familia viscoelástica). Extensiones
        de `muscle_meta`: bindings OO con dimensión (4-tuplas), `oo_class` (argumentos `class(...)` en
        el adaptador) y parámetros tensoriales de dimensión fija (`'3D4O3sym'`). Parámetros derivados
        privados con `elas_linear_init`.
    (e) (Decisión del usuario tras la revisión, 2026-11-01) **Caché del tangente sólo en el adaptador
        OO:** el tipo de kernel `elas_linear_<k>` guarda sólo λ, μ y e1 y el tangente se ensambla al
        vuelo (mismos bits que el caché); los tangentes cacheados y los flags `calc_tan_3D/2D` quedan
        en el tipo OO histórico. La maquinaria admite componentes sólo-OO en el adaptador.
13. (Golden viscoelasticity, decisiones del usuario 2026-11-01)
    (a) **γ de Maxwell:** γ = (1 − e^{−x})/x, x = dt/τ: el cambio de la serie (3 términos) a la fórmula
        directa pasa de x < 1e-5 a **x < 1e-4 en dp** (umbral óptimo ≈ (24ε)^{1/4}: la fórmula directa
        cancelaba ~3.7e-12 en x = 1e-5, por encima del 1e-12 de testing.md, con un salto en el cambio);
        en sp el umbral se elige por la misma regla. **Cambio aprobado en dp** sólo para
        1e-5 ≤ x < 1e-4; golden de esos registros actualizado.
    (b) **dt ≤ 0 con status y OO unificado:** la API de kernel devuelve un código `status` (constantes
        `STATUS_*`) en vez de abortar: dt < 0 es error en Maxwell y Kelvin; dt = 0 es válido en Maxwell
        (respuesta instantánea) y error en Kelvin. Los adaptadores OO validan igual: dt < 0 también
        aborta en Maxwell en host, como Kelvin. **Cambio de comportamiento aprobado** (Maxwell con
        dt < 0: de Inf/NaN a error).
    (c) **Sólo se documentan:** Maxwell no valida parámetros negativos ni ν = 0.5 (Kelvin sí); Kelvin
        pierde dígitos con dt muy chico (diferencia de términos ∝ ζ/dt); `stress_2D` de
        `Elasticity_frozen` ignora yz/xz de la tensión histórica; las `_transient` ignoran
        dstrain/dtime/internal; `unpack_from_fea` exige `init` previo; los modelos y el estado son
        sólo 3D. Las ramas ya son un array fijo (`MAX_VISCO_BRANCHES = 10`): el pendiente del PLAN
        (`MAXWELL_MAX_BRANCHES`) ya está resuelto.
    (d) Diseño de la familia `visco` (piloto, 2026-11-01): primera familia con estado. Tipos por kind
        de tamaño fijo en `muscle_visco_types`: `visco_state_<k>` (n, strain_e, h(6, MAX)) y
        `visco_coefs_<k>` (dt, k_alg, g_alg, decay, weight); los nombres históricos son alias dp.
        API de kernel (sólo 3D): `visco_step_coefficients(prm, dt, coefs, status)` una vez por
        material y paso (único punto que valida); `visco_freeze(prm, coefs, state_n)` → `elas_frozen_<k>`
        por punto; `visco_update(prm, coefs, state_n, strain_e, state_np1)` con `state_n` intent(in)
        y `state_np1` intent(out) (las ramas inactivas se copian de `state_n`). El genérico se llama
        `visco_step_coefficients` (decisión del usuario) para no chocar con el tipo histórico
        `Visco_coefficients` (Fortran no distingue mayúsculas). La **elasticidad congelada** es el
        modelo `frozen` de la familia `elas` (un solver escrito contra `elas_*` sirve también para
        visco), con **parámetros por punto** (k_alg, g_alg, stress_hist públicos, sin init): en la Fase 7
        el `prm` elástico de un material visco es local del hilo, no un dato mapeado una vez. Parámetros
        de Maxwell/Kelvin privados porque `init` valida (amplía D-017.10b de "derivados" a "validados").
        **Status común** en `src/common/muscle_status.f90`: `STATUS_OK = 0`,
        `STATUS_INVALID_PARAMETER = 101`, `STATUS_INVALID_TIME_STEP = 102`; con error las salidas quedan
        con sus valores por defecto. γ de Maxwell en sp: serie para x < 1e-2. Extensiones de
        `muscle_meta`: `ModelFunc.outs` (subrutinas), tipos `'type:<name>'`, `MODEL_MODULES`, `support`
        y defaults de array/tensor.
14. (Golden kinematics, decisiones del usuario 2026-11-01)
    (a) **Status en la API de kernel; OO igual:** la API de kernel de las cinemáticas finitas valida
        y devuelve `STATUS_*` (nuevo código de deformación inválida en `muscle_status`) cuando
        det F ≤ 0, cuando la inversa de F falla (`ok` de `inverse`) o cuando un autovalor de C queda bajo
        el umbral 1e-12 (antes: resultados numéricos erróneos sin aviso, p.ej. R no ortogonal o ln λ
        perdido con F casi singular). Los adaptadores OO conservan el comportamiento histórico
        (mismos bits).
    (b) **Sólo se documentan:** la pérdida de precisión de la rotación polar R cerca de C ≈ I
        (|RᵀR − I| ≈ 3.3e-11 con F ≈ I + 1e-3, por los proyectores de Cayley-Hamilton con autovalores
        casi repetidos); `stress_spatial`/`stress_material` se ignoran en los tangentes (sin término
        geométrico ni de Jaumann); `to_spatial_tangent` de la logarítmica material es R·C·R/J sin la
        derivada de ln C; `update_incremental` ignora dt y es aditiva; `Corotational%update` sin R
        conserva el R anterior sin usarlo; `update_F_base` deja `strain` desactualizado. Las funciones
        privadas sin uso `finite_get_E_log` y `finite_get_b` se conservan.
    (c) Diseño de la familia `kin` (piloto, 2026-11-01): sin parámetros físicos: tipos de parámetros
        vacíos `kin_<modelo>_<k>` que sólo eligen la cinemática en compilación (mismo contrato `prm`
        que las demás familias; warnings de `prm` sin usar aceptados). Estado `kin_state_<k>` de tamaño
        fijo (F, F_inv, R, J, C, b, b_inv, strain; por defecto la configuración de referencia).
        `kin_update(prm, F, state, status)` es la única función que valida (STATUS_INVALID_DEFORMATION
        = 103: inversa fallida, J ≤ 0 o NaN, autovalor de C ≤ `KIN_EIG_MIN_<k>`: 1e-12 dp, **1e-4 sp**);
        medidas, conversiones PK1/PK2/Cauchy/Kirchhoff (en minúsculas, `kin_pk2_to_cauchy`),
        push-forward/pull-back de tangentes y rotaciones corrotacionales como funciones comunes del
        estado. Módulo de familia `muscle_kin_laws` (la fachada OO `muscle_kinematics` sigue exportando
        sólo lo histórico hasta la Fase 7). Small strain sólo escribe `strain`.
    (d) **Corrotacional:** función propia de la familia `kin_corot_update(prm, strain, state, status
        [, R])` sobre el mismo `kin_state`; sin R → R = I y `has_rotation` falso, como hoy.
    (e) **Logarítmica espacial:** el status se decide sólo con los autovalores de C (los de b son los
        mismos); se reutiliza la descomposición si no cambia bits.
    (f) **Convención común de error (todas las familias):** con status ≠ `STATUS_OK` las salidas no
        deben usarse (pueden quedar calculadas o en sus valores por defecto, según convenga a cada
        familia para conservar bits; visco las deja por defecto y kin calculadas). El solver de la
        Fase 7 comprueba siempre el status. Documentado en `muscle_status`.

## D-018 — Fase 7: alcance y diseño de los materiales declarados
**Estado:** Aceptada (2026-11-01, decisiones del usuario)
**Contexto:** plan de la Fase 7 (D-008). Hay cinco solvers OO, sólo 3D y sólo dp, y ninguno tiene golden
(`closest_point`, `closest_point_perzyna`, `closest_point_time`, `closest_point_2`, `cutting_plane`).
**Decisión:**
1. **Se portan los cinco solvers.**
   - `closest_point_2` lleva LU propio en vez de LAPACK (D-016).
   - `cutting_plane` se porta como algoritmo explícito.
   - Ambos se limpian de `print`/`write` (no hay salida en el kernel).
2. **`closest_point_time` se porta tal cual,** con los mismos bits en dp. Sus bugs se documentan en
   PROBLEMAS.md: `strain_p` sin inicializar en el tangente elástico, chequeo elástico con epd = 0 y
   `strain_p` sin omega.
3. **Tolerancia del Newton relativa en ambos kinds** (|f| < tol·σ_y), con valor por defecto por kind y
   configurable en el material. Valores por defecto: `sp` 1e-5 (≈ 84 ε) y `dp` 1e-10 (hoy 1e-5
   absoluto ≈ 3e-8 relativo con σ_y ≈ 300 MPa). El solver OO `closest_point` recibe el mismo
   criterio, para que OO y kernel den los mismos bits hasta el retiro.
   - Es un cambio aprobado en dp: cambian iteraciones y bits.
   - El golden se captura antes con el código original y después se actualizan sólo los registros
     afectados («Updated (D-018.3, approved)»).
4. **Alcance:** todos los modelos, incluidas las cinemáticas finitas (TL/UL/logarítmicas/corrotacional)
   y la elasticidad visco congelada.
5. **Cuerpo único de solver:**
   - Llama a los genéricos de familia `<fam>_<func>` de la Fase 6, que se resuelven en compilación por
     el tipo de parámetros.
   - Por material sólo cambian los tipos de los componentes.
   - **Enmienda D-008.1:** no hay renombre en `use`.
   - El cuerpo es un `.fypp` incluido por la plantilla de material con `k` fijado (D-009, D-017.6).
6. **Tangente:** analítico en el kernel. Las diferencias finitas quedan para tests/host, como
   verificación y sin closures.
7. **API por lotes en SoA:** `eps(npts,6)`, `sig(npts,6)` y `hsv(npts,NSTATE)`.
   - La tensión va en `sig`, no en el estado.
   - Para el multifase se pasa la lista de índices `idx(:)`.
8. **Retiro del OO:** los `Base_*`, adaptadores OO, la API `procedure(...)`, `*_host` de closures y
   `Plastic_material_history` se retiran al final de la Fase 7, cuando todos los tests de solvers usen
   la API nueva contra el golden.
9. **Registro de materiales:** archivo separado `tools/fypp/muscle_materials.py`.
   - Trae por defecto un registro de prueba que cubre cada modelo.
   - `-DMUSCLE_MATERIALS=<ruta>` lo reemplaza.
   - Lo valida `muscle_meta`.
10. **Anomalías del golden de solvers (2026-11-01):**
    - (a) `closest_point`: con iteraciones agotadas reporta `min(iters, iter_nw)`, como perzyna (cambio
      aprobado en los registros no convergidos).
    - (b) `closest_point_time`: el chequeo elástico pasa a f = q − σ_flujo ≤ tol, sin división. Hoy
      σ_flujo(ep, 0) = 0 da NaN. Es un cambio aprobado que afecta sólo a los registros hoy NaN o no
      convergidos; el resto de `time` queda tal cual (D-018.2).
    - (c) `closest_point_2` sin dγ ≥ 0 (ep negativo con incremento grande): documentado, port fiel.
    - (d) `closest_point_time` + RK: el tangente analítico difiere ~8 % del numérico (NNL 3.5e-4):
      documentado, sin investigar.
    - Las salidas en pasos no convergidos (tensión 0 en `time`/`cp2`) quedan cubiertas por D-017.14f.
11. **Piloto `closest_point` (2026-11-02):**
    - (a) El residuo de deformación plástica del Newton usa la misma tolerancia por kind (adimensional),
      en vez de 1e-5 absoluto. Es un cambio aprobado en dp, incluido en los registros
      «Updated (D-018.3/D-018.10a)».
    - (b) El denominador de |f|/σ_y es `max(σ_y, SY_MIN)`, para que un ablandamiento hasta cero no
      bloquee la convergencia.
    - (c) `TOL_ELASTIC` por kind en `muscle_kinds`: `dp` 1e-8 (histórico) y `sp` 1e-5. Se mantiene
      «iters = 0 ⇔ elástico».
    - (d) El hsv del material no lleva la tensión: `MAT_<NAME>_NSTATE = 7` (strain_p 6 + ep 1), con
      pack/unpack propio del material (D-018.7). Con elasticidad visco congelada son
      7 + 6 + 6·`visco_branches` (D-018.14).
    - (e) Los tests históricos OO (hill48, cpb06, cpb06ex2, yld2004) pasan a `iters ≤ 7` por D-018.3.
    - (f) Nuevo status `STATUS_NOT_CONVERGED = 104`. Hay un genérico `mat_update` común a todos los
      módulos de material (los específicos difieren en el tipo). `golden_finish(partial)` sirve para
      tests que comparan un subconjunto de un golden.
    - (g) Materiales sólo 3D: la familia yield es sólo 3D. El 2D queda fuera de la Fase 7.
12. **`perzyna` y `time` (2026-11-02):**
    - (a) Los residuos de Perzyna se escalan a tensión: C:R_εp y R_λ/(Δt·∂g/∂σ_eq) (error de
      sobretensión, independiente de η) se comparan con tol·max(σ_y, SY_MIN). R_λ·n:C:n se descartó:
      con η chico el redondeo de σ/σ_y − 1 se amplifica por dt/η·n:C:n y no converge.
      C:R_εp tiene un piso de redondeo: max(tol·max(σ_y, SY_MIN), 8·ε_k·max|C_ij|·|εp|)
      (`RES_FLOOR_ROUNDINGS`). Así sp converge con incrementos grandes (ep ≈ 20× la deformación de
      fluencia). En dp no cambia ningún registro. Así se corrige la precisión de sp
      (con residuos en deformación y tol 1e-5 el error llegaba a 5e-3). Cambio aprobado en dp, aplicado
      también al OO.
    - (b) `jcfull`: el OO usa el `dstress_dep` analítico del kernel en lugar del de diferencias finitas
      (D-017.8), con los mismos bits OO/kernel. Se actualizan los registros jcfull del golden de
      `time` y sus `dstress_dep` en `vp_laws_dp.txt` (18 registros). En sp no se comprueba jcfull: el
      paso 1 no converge y después difiere ~6e-4, porque el flujo es 0 a tasa nula (D-017.7).
    - (c) Firma uniforme: todos los `mat_update` reciben `dt` después de `strain`, y `closest_point`
      lo ignora.
    - (d) El chequeo elástico de `time` (D-018.10b) usa `TOL_ELASTIC·max(σ_flujo, SY_MIN)`.
      `SY_MIN_<k>` se define por kind en `muscle_kinds`.
    - (e) El pack/unpack de hsv es común (`muscle_plastic_state_hsv.fypp`, hoy `src/materials/muscle_material_kernel.fypp`). `MaterialSolver` declara
      las familias compuestas (`composes`) y si el solver usa `dt`.
    - (f) Los goldens comparan NaN = NaN sin mirar signo ni payload. El signo del NaN depende del
      orden de operandos que elige el compilador.
13. **`closest_point_2` y `cutting_plane` (2026-11-02):**
    - (a) `cp2`, en OO y en kernel, resuelve el 7x7 con `LU7SOLVE` en lugar de LAPACK `dgesv` (D-018.1,
      D-016.10). Con LU singular deja el iterado sin cambiar y devuelve `STATUS_NOT_CONVERGED`.
    - (b) `cp2` usa la tolerancia relativa por kind, tanto en f como en el residuo de deformación
      plástica (como D-018.11a). Su chequeo elástico pasa del histórico −1e-16 a
      `TOL_ELASTIC_k·max(σ_y, SY_MIN_k)`, igual que `closest_point`. Cambio aprobado en dp.
    - (c) El tangente de `cp2` es el tangente consistente de `closest_point` (complemento de Schur
      compartido en `muscle_closest_point_tangent.fypp`). La guarda es `|dγ| ≤ 1e-12`, para que el
      estado con dγ < 0 (D-018.10c) reciba el tangente consistente. En `closest_point` no cambia nada,
      porque allí dγ ≥ 0.
    - (d) `cutting_plane` devuelve el tangente continuo elastoplástico C − (C:n⊗n:C)/(n:C:n + h')
      (el original no tenía tangente), verificado contra diferencias finitas con la tolerancia que
      corresponde a la diferencia entre continuo y algorítmico. Esa diferencia llega a 0.38 de la
      componente mayor, así que el chequeo por diferencias finitas (tol 1.2) sólo detecta errores
      groseros. Sólo el kernel tiene tangente: el `cutting_plane` OO no tiene salida de tangente.
    - (e) `MaterialSolver.iter_nw` es el límite de iteraciones por solver: 200 por defecto, 100 en
      `cutting_plane` (valor histórico).
14. **Materiales con elasticidad visco congelada (2026-11-02):**
    - (a) El registro tiene `elas='frozen'`, `visco='maxwell'|'kelvin'` y `visco_branches` (máximo de
      ramas por material, ≤ `MAX_VISCO_BRANCHES`, que ahora vive en `muscle_meta.py`).
    - (b) El hsv es 7 + 6 (strain_e) + 6·`visco_branches`. El estado es `visco_plastic_state_<k>`
      (composición de `plastic` + `visco`). Al hacer unpack, `n_branches = visco_branches`; los
      kernels usan `prm%n`.
    - (c) Los cuerpos de solver reciben `elas` como argumento (`SOLVE(mat, elas, ...)`), sin copiar
      el material por punto.
    - (d) `visco_step_coefficients` se calcula **una vez por lote/paso** (D-017.13d), no por punto,
      y el kernel de punto recibe los coeficientes. Ahí se comprueban también el número de ramas
      (`STATUS_INVALID_PARAMETER`) y dt.
    - (e) Orden por punto: `visco_freeze` → `SOLVE` → `visco_update(strain − strain_p)`, igual que
      el golden. `visco_update` se llama también si no converge (D-017.14f).
    - (f) Queda pendiente para la Fase 8, después de perfilar: el kernel sigue usando
      `visco_state` con `h(6, MAX_VISCO_BRANCHES)`; `visco_branches` sólo reduce el hsv.
15. **Materiales con cinemáticas finitas (2026-11-02):**
    - (a) El kernel recibe F_n y F_{n+1} (DFGRD0/1). La firma de pequeñas deformaciones no cambia.
      El cuerpo del solver trabaja igual sobre `ks%strain` (E, e, ½ln C o ½ln b), y `strain_p`/ep se
      guardan en esa medida.
    - (b) El término geométrico/de tasa se agrega ahora: función de familia
      `kin_geometric_tangent` y campo del material `tangent_rate = 'none'|'jaumann'`, verificado
      contra diferencias finitas de σ(F).
    - (c) Se corrige el tangente espacial logarítmico material para incluir la derivada de ln C.
      Cambio aprobado: se actualizan los registros de kin afectados y el adaptador OO da los mismos
      bits.
    - (d) Corrotacional: R_n se guarda en el hsv (+9 reales sólo en esos materiales).
    - (e) Referencias: no hay composición histórica de solver con cinemática, así que se usan
      oráculos físicos y un golden dp del kernel nuevo como regresión:
      - con F → I coincide con pequeñas deformaciones;
      - es invariante ante rotación rígida;
      - tracción uniaxial log + von Mises contra la solución analítica;
      - tangente contra diferencias finitas.
    - (f) **Enmiendas tras el análisis con diferencias finitas (2026-11-02):**
      - Se corrigen los tangentes espaciales de **todas** las cinemáticas, como cambio aprobado igual
        que (c):
        - UL: le faltaban ∂e/∂ε y σ⊗I (34 % de error);
        - log espacial: le faltaba d ln b/db;
        - log material: le faltaban la derivada de ln C y el giro de R.
      - `kin_geometric_tangent` = (1/J)∂τ/∂ε a tensión constitutiva fija:
        - TL: el término de σ;
        - UL: σ⊗I;
        - log espacial y pequeñas deformaciones: 0;
        - log material: el término de dR.
      - `tangent_rate='jaumann'` = `to_spatial_tangent` + `kin_geometric_tangent` = (1/J)∂(Jσ)/∂ε con
        δF = sym(δε)·F (tasa de Jaumann de Kirchhoff, Abaqus/ANSYS). `'none'` = `to_spatial_tangent`.
      - ~~Corrotacional: medida total ln U de la polar de F_{n+1}.~~ Reemplazado por la «Corrotacional
        rotada por el FE» (más abajo). **Se descarta (d):** no se guarda R_n.
      - Medidas espaciales (UL, log espacial): `strain_p` se rota con ΔR = R_{n+1}·R_nᵀ en cada paso
        (formulación objetiva). R_n sale de la polar de F_n, que llega como argumento, así que no hay
        hsv extra. El oráculo de rotación rígida se aplica a todas las medidas.
      - Criterios anisótropos con medidas espaciales (UL, log espacial): la anisotropía se evalúa en
        el marco rotado, f(Rᵀσ R), con R de la polar de F_{n+1}, y el gradiente y la hessiana se
        llevan de vuelta con R. Así la formulación es objetiva con cualquier criterio. Los criterios
        isótropos (von Mises, Drucker-Prager) no rotan, para conservar los bits. Exige elasticidad
        isótropa, y el validador lo comprueba. En el registro de
        `kin_geometric_tangent` de UL, σ⊗I va en el término geométrico; `to_spatial_tangent` sólo
        lleva ∂e/∂ε.
      - **Cinemáticas en sp (2026-11-02):** la polar y el logaritmo espectrales fallan en sp con
        deformaciones chicas: R no es ortogonal (|RᵀR − I| = 7e-4 con 4e-4 de deformación) y ln C pierde
        la parte desviadora cuando las brechas entre autovalores son ≲ 1e-2. Se usan algoritmos
        robustos **sólo en sp**, y dp conserva sus bits:
        - R por iteración de Newton/Higham con un número fijo de iteraciones;
        - ln C / ln b y sus derivadas con una descomposición espectral de Jacobi (6 barridos fijos);
          la serie (`log_ratio`) sólo se usa en el tangente cuando los autovalores están cerca.
      - **Corrotacional rotada por el FE (2026-11-02, reemplaza «medida total ln U»):** como en
        Abaqus y ANSYS, el FE entrega la tensión y la deformación ya rotadas al marco corrotacional.
        - El material corrotacional es el kernel de pequeñas deformaciones alimentado con la
          deformación rotada, más el tangente de Jaumann.
        - `strain_p` del hsv la rota el adaptador UMAT con `DROT` (paso 8): Abaqus no rota `STATEV`.
          Falta confirmar cómo trata ANSYS sus `ustatev`.
        - Se elimina el ln U total calculado con la polar de F en el kernel.
      - Cobertura de solvers con cinemáticas: `closest_point`, `perzyna`, `time` y `cp2`. El validador
        rechaza las cinemáticas con `cutting_plane`, porque su tangente continuo no es consistente
        (D-018.13d), y también con elasticidad visco.
16. **API por lotes (2026-11-02):**
    - (a) La API es `muscle_update_<name>(mat, eps | F_n, F_np1, dt, hsv_n, hsv_np1, sig, status, iters
      [, ddsdde][, idx])`, en SoA, con arrays `contiguous` y `npts = size(sig,1)`.
      - F va por filas en el orden `ten_3D2O%vals`, column-major, igual que `reshape(DFGRD,[9])`.
      - `ddsdde(p,i,j) = ∂σ_i/∂γ_j`, con corte ingenieril en j = 4..6, igual que DDSDDE de UMAT.
      - Corrotacional recibe `eps`.
      - Con `idx`, el resto de las filas no se toca, y `idx` no puede tener índices repetidos.
    - (b) Los opcionales (`ddsdde`, `idx`) se resuelven en el host, con 4 bucles privados por kind
      generados desde una sola fuente. `iters` es obligatorio. Ningún opcional entra en la región de
      dispositivo (decisión del usuario).
    - (c) Macro `BATCH_LOOP`:
      - en CPU expande a `!$omp parallel do`;
      - con `MUSCLE_GPU` expande a `!$omp target teams distribute parallel do map(present, alloc: mat,
        arrays...)`, más `map(to: coefs)` en los materiales visco.
      - `mat` también tiene que estar presente: quien llama la mapea una vez (decisión del usuario,
        coherente con D-008.2).
      - gfortran 15 rechaza `map(present: x)`; la forma OpenMP 5.1 `present, alloc` funciona.
    - (d) Visco: si `mat_step` falla, cada punto a calcular recibe ese status con iters = 0 y no se
      escribe nada más. Se hace dentro del bucle porque status/iters viven en el dispositivo.
    - (e) No hay chequeos de tamaño en tiempo de ejecución: el contrato de extensiones, rango de `idx`
      y contigüidad queda documentado en la API.
    - (f) Pendiente para la Fase 8:
      - camino `is_device_ptr` para memoria de dispositivo pura (`omp_target_alloc`, cuFFT), que no
        está en la tabla `present`;
      - verificar el modificador `present` en nvhpc. OpenMP 6.0 depreca `alloc` en favor de `storage`.
17. **API genérica para integración con programas FE (2026-11-02; reemplaza D-008.3 «adaptadores
    UMAT generados por material»):** la librería no trae adaptadores. Cada usuario escribe el suyo
    para cualquier programa FE, con estas piezas generadas por material:
    - (a) `mat_init_<name>(mat, props(:), status)`: el orden de parámetros se genera por material,
      con `MAT_<NAME>_NPROPS` y una tabla del orden en la documentación. Llama a los `*_init`.
    - (b) Una API de punto con arrays planos, `muscle_point_<name>(mat, strain(6) | F_n(9), F_np1(9),
      dt, hsv_n, hsv_np1, stress(6), status, iters[, ddsdde(6,6)])`, sin tipos derivados de tensores.
      Usa las mismas convenciones que el lote (D-018.16a).
    - (c) Un helper para rotar las partes tensoriales del hsv con una rotación incremental (DROT u
      otra), para los FE que no rotan sus variables de estado (D-018.15f, corrotacional).
    - (d) Ejemplos en `pages/integration.md`: esqueletos de UMAT (Abaqus), USERMAT (ANSYS) y umat41
      (LS-DYNA), sin compilar ni probar contra el FE. La duda de los `ustatev` de ANSYS queda
      anotada en el ejemplo.
    - (e) Parámetros (decisión del usuario, 2026-11-02):
      - `muscle_meta` es la **fuente única**: los `PARAMS` de cada modelo y, cuando difieren, las
        entradas de su init. Las plantillas los leen de ahí, con los mismos bits.
      - Orden de props:
        - Maxwell: bulk, shear_inf y luego `visco_branches` pares (shear, tau);
        - `hardening_mode` va como real redondeado con `nint`;
        - `T_current` de jcfull cuenta como prop;
        - `tol` e `iter_nw` van al final.
    - (f) `mat_rotate_hsv_<name>` sólo se genera para materiales con entrada de deformación (pequeñas
      deformaciones y corrotacional). Los que reciben F ya son objetivos (ΔR).
    - (g) La API de punto de los materiales visco recibe `coefs` de `mat_step`, igual que `mat_update`.
      Las notas de LS-DYNA en los ejemplos (`IHYPER`, `IBULK`/`IG`, sin rotación en `umat41`) están
      sin verificar.
    - (h) `mat_init_<name>` (2026-11-02):
      - Con menos props que `MAT_<NAME>_NPROPS` devuelve `STATUS_INVALID_PARAMETER`; los sobrantes se
        ignoran.
      - `tol` ≤ 0 o `iter_nw` ≤ 0 también son error, sin valor por defecto (decisión del usuario).
      - Valida además lo que puede comprobar barato: el `hardening_mode` de DP fuera de 1..3, ν fuera
        de (−1, 0.5) y E ≤ 0. Los init de los modelos y el OO no cambian.
      - Para los materiales visco congelados, los componentes van en el orden del tipo del material:
        el modelo visco ocupa el lugar de elas.
18. **Retiro del OO (paso 9, 2026-11-02; detalla D-018.8):**
    - (a) Los casos del golden de solvers sin material en el registro (cp 12/22, perzyna 7/16) se
      cubren añadiendo unos 15 materiales al registro de prueba. No se borra ningún registro de
      solver.
    - (b) Se borra la API de math con `procedure(...)` (`derivative*`, módulos `*_host`) y quedan las
      macros de `include/muscle_fd.fypp`. Se quitan las exenciones host-only del lint que ya no
      hagan falta. Los goldens de derivadas se generan con las macros, y si el compilador integra la
      función y cambian bits se documenta una tolerancia (D-016.9).
    - (c) No hay capa host de conveniencia. Se publica en `pages/` una tabla de equivalencias de API
      vieja a nueva, con lo que se pierde:
      - commit/rollback transaccional (el FE guarda `hsv_n`);
      - elegir el modelo en tiempo de ejecución (queda la unión etiquetada opcional);
      - modelos propios por extensión de clase;
      - `set_parameters`/`to_string`.
    - (d) Se borra la fachada `muscle_kinematics`.
    - (e) Se mantienen `muscle_tensors`, `muscle_tensor_functions`, los `*_tbp` de tensores y
      `muscle_tensor_io` (UDTIO).
    - (f) Goldens: se conservan los `.txt`, y cada generador pasa a la API kernel con las mismas
      claves. Sólo se borran los registros sin equivalente en kernel, con una línea
      «Removed (D-018.8)»:
      - elas: `to_string`, flags de caché y defaults de `Base`;
      - yield: `base numeric`;
      - vp: defaults de diferencias finitas de las bases;
      - visco: commit/rollback.
      `plastic_history_dp.txt` se retira entero, porque `pack/unpack` ya se compara en
      `plastic_state_kinds`. En cada commit de migración `regenerate.sh` tiene que dejar los `.txt`
      sin diferencias, salvo las bajas declaradas.
    - (g) Las pruebas físicas que sólo existen en OO se portan al kernel antes de borrar:
      - valores r y cortante de los criterios anisótropos;
      - incremento grande;
      - relajación exacta y η→0 en Perzyna;
      - radial cerrado;
      - relajación analítica y `time_step_limits` en visco;
      - roundtrips de cinemáticas;
      - Chaboche.
      La prueba `joint` (extensión de clase) se borra.
    - (h) Golden de `time`: el tangente numérico del OO volvía a resolver arrancando en caliente desde la
      iteración convergida (dλ, ω), y el kernel no lo hace. Los registros `tangent_num` pasan a
      diferencias finitas desde t_n; cambian 80 de 136, ≤ 1.5e-7 relativo salvo el caso 14, que no
      converge. Cambio aprobado: «Updated (D-018.18h, approved)».
    - (i) Antes de borrar los solvers OO se portan a materiales los tests de acople que los usaban:
      visco Maxwell/Kelvin + closest_point y el tangente mixto de Chaboche.
    - (j) Golden de yield: se borran los 15 registros `constructor` de cpb06 (6), cpb06ex2 (5) y
      yld2004 (4), porque leen parámetros derivados que en el kernel son privados (D-017.10b). Ya pasaba
      lo mismo con Drucker-Prager, y los registros analíticos siguen detectando un init erróneo.
    - (k) Golden de elas: además de lo listado en (f), se borran los 413 registros de las sobrecargas
      transient y de `stress_and_update`, porque el kernel no tiene API equivalente y los registros
      no transient que quedan ya cubren esos valores. En total se quitan 962 de 1311.
    - (l) Golden de cinemáticas: se borran los ~247 registros de comportamiento sólo OO:
      - push/pull identidad de pequeñas deformaciones;
      - corrotacional «sin R» y copia polimórfica;
      - `update_F_base`;
      - defaults.
      Además se quita `has_rotation` de los 20 registros de estado corrotacional. El resto sale del
      kernel, bit a bit.
    - (m) Math (E4): se borran los 20 registros de diferencias hacia adelante de matrices 3x3
      (`derivative_escalar_tensor3x3`, `_3x3sym2`) de `math_operations_dp.txt`, porque no hay macro
      equivalente y ningún kernel las usa. También se borra `FD_STEP_FWD_<k>`, que queda sin usos.
      Los goldens de derivadas se reproducen bit a bit con las macros (`-fno-inline`, D-016.9).
19. **Bugs del smoke test (2026-11-03):**
    - (a) **dp, logarítmicas cerca de F = I.** (**Reemplazada por D-032**, 2026-11-12: `eigenprojections` usa
      `spectral_decomposition`, precisa (pocos ε) en todas las brechas, y la franja `JACOBI_GAP_LO/HI` se retiró.) La forma
      cerrada (autovalores + proyectores de
      Cayley-Hamilton) toma las raíces como repetidas cuando |D_a| ≤ 1e-11·λmax². Con brechas
      < ~3e-6 queda P = I/3 y ln C pierde la parte desviadora (σ = 0 en corte puro).
      - Corrección **híbrida** (decisión del usuario): en dp se usa Jacobi sólo cuando la brecha
        relativa mínima entre autovalores cae en una franja intermedia; fuera de ella sigue la forma
        cerrada. La franja propuesta es [~1e-12, ~1e-3); se mide y se fija como constante en
        `muscle_kinds`.
      - Las brechas de redondeo (~1e-16, casos rotados casi repetidos) quedan en la forma cerrada,
        donde I/3 es exacto.
      - sp sigue siempre con Jacobi.
      - Los registros afectados se actualizan como cambio aprobado.
    - (b) **Newton de `time` en el piso de tasa de JC.** `vp_jc_rate` da derivada 0 en el piso, el
      primer paso se pasa y el recorte con ω deja una convergencia lineal (124–200 iteraciones). En
      el piso se usa la **derivada por la derecha**. Cambian sólo los registros de `time` cuya
      historia tocaba el recorte, más sus tangentes.
    - (c) **JC con flujo negativo en el piso:** `mat_init_<name>` devuelve `STATUS_INVALID_PARAMETER`
      si C·|ln(epmin/epdmax)| ≥ 1 (regla de dominio en `muscle_meta`). Los init de modelo no cambian.
    - (d) Aplicación (2026-11-03):
      - Franja de Jacobi en dp: `KIN_JACOBI_GAP_LO/HI_dp` = [1e-12, 1e-2), medida contra una
        referencia en cuádruple precisión.
        - Con la forma cerrada, el error de ln C es ≈ eps/gap², y por debajo de ~1.8e-6 pierde todo
          el desviador.
        - Jacobi da error ≈ eps/gap y RᵀR − I ≤ 1.4e-15.
        - Saltos en los bordes: σ ≤ 1.1e-11 y tangente ≤ 2.6e-11 relativos.
        - sp: [0, huge), es decir, siempre Jacobi.
      - JC: la ley pública conserva su derivada verdadera, 0 en la región plana. El solver `time`
        usa la derivada por la derecha sólo para la dirección de Newton (decisión del usuario).
        - Iteraciones: jc1 124–129 → 7–8, jc3 130–178 → 7.
        - jc2 sube de 4–9 a 6–11 iteraciones; aceptado.

## D-019 — Fase 8: build GPU
**Estado:** Aceptada (2026-11-03, decisiones del usuario)
**Contexto:** el único hardware disponible es una Intel UHD 770 (ifx/spir64): sp nativo, dp emulado,
que falla con kernels grandes. No hay GPU NVIDIA ni AMD. Los backends de offload de gfortran 15
(`gcc-15-offload-nvptx`/`amdgcn`) ya están instalados.
**Decisión:**
1. **Verificación:**
   - En la Intel se ejecutan los tests target vs CPU en sp, y en dp sólo en los kernels chicos que
     la emulación soporta.
   - Para nvptx y amdgcn se exige compilar y enlazar.
   - La ejecución en NVIDIA/AMD la hace el usuario en otra máquina.
2. El prototipo FFT → MUSCLE → iFFT pasa a la Fase 9.
3. `tools/check.sh` compila **siempre** los tres targets (nvptx y amdgcn con enlace; ifx spir64 con
   build y ctest de los tests target en la GPU, con `systemd-run` y límites de memoria y tiempo).
4. **Paso 1 (2026-11-03):**
   - (a) Flags (`cmake/MuscleGpu.cmake`, PUBLIC en `muscle`):
     - NVPTX: `-fopenmp -foffload=nvptx-none`.
     - AMDGCN: `-foffload=amdgcn-amdhsa -foffload-options=amdgcn-amdhsa=-march=gfx90a`, con
       `MUSCLE_GPU_ARCH` en caché.
     - gfortran enlaza también `-lgomp` con `--no-as-needed`.
     - SPIR64: `-fiopenmp -fopenmp-targets=spir64`.
   - (b) `check.sh` en los builds GPU construye sólo la librería y los tests target, con NVPTX y AMDGCN
     en paralelo (decisión del usuario). El enlace completo de todos los tests queda en
     `check.sh --gpu-full`.
   - (c) El `matmul` vector×matriz de `2Osym : 4O2sym` se escribe como sumas explícitas **siempre**,
     porque ifx genera código de dispositivo roto. Es un cambio aprobado en dp: ≤ 2.2e-16 relativo a la
     componente mayor de cada resultado (no ≤ 1 ulp por componente: en una componente chica de 2D llega
     a 144 ulp de esa componente). Cambian 6 registros del golden de contracciones.
   - (d) fypp corre con `--line-numbering-mode=nocontlines`, porque ifx rechaza un `#line` entre las
     mitades de un literal partido. Enmienda D-009; el código CPU no cambia.
   - (e) Pendientes del paso 2:
     - ifx con fallback en host da NaN en las regiones target de `kin_kinds`;
     - los `mat_batch*` necesitan datos ya mapeados (`map(present)`);
     - `math_kernels_target` no carga con dp emulado (`__isinfd`/`__isnand` sin resolver).
   - (f) ifx 2026.2 da un error interno al compilar 5 tests `tensors-types-*_kinds`: `call p%init`
     sobre un tipo exportado con dos nombres por `muscle_tensors`. Documentado en PROBLEMAS.md.
5. **Paso 2 (2026-11-03):**
   - (a) Tolerancia en dispositivo `TOL_DEVICE_<k>` = 64 ε por punto, multiplicada por la amplificación
     del punto cuando el kernel la tiene:
     - κ∞ en LU;
     - κ(C) en cinemáticas;
     - 1/h en diferencias finitas.
     En host fallback se exigen los mismos bits. Errores medidos en la Intel:
     - eigen 1.7/2.0 ε;
     - LU 0.1 ε·κ∞;
     - cinemáticas ≤ 7.4/19 ε·κ(C).
   - (b) `exponent()` se reemplaza por `fp_exponent`, leído de los bits IEEE. Da los mismos bits en CPU.
     ifx implementa `exponent` con `for_exponent8_v`, que no existe con dp emulado. Un lint lo prohíbe.
   - (c) `mat_batch*` se excluyen en la GPU Intel y se documentan (decisión del usuario): IGC compila mal
     closest_point/closest_point_2 con más de un work item. NVIDIA/AMD los corre el usuario.
   - (d) `mat_kernels_target` se divide en un test por material, porque la imagen de 12 materiales
     supera los 4 GB / 300 s del JIT.
   - (e) Se agregan tests target directos, generados por fypp, para las funciones libres de tensores y
     los genéricos de familia (decisión del usuario).
   - (f) Aplicación de (d) y (e):
     - `tests/gpu/target_cases.py` define los casos y genera una imagen chica por test, con label
       `target` en ctest y target `build_target_tests`.
     - Tests target generados: 20 de tensores (5 familias × 2 dims × 2 kinds), 12 de genéricos de
       familia y 24 de materiales (12 materiales × kind). En CPU son 172 tests, y en host fallback los
       target dan los mismos bits.
     - En la GPU Intel corren 38 tests, dentro de 64 ε. Se excluyen por IGC:
       - los materiales closest_point, cp2 y time en sp (device lost u otros iterados);
       - todos los materiales en dp (fallos no reproducibles con dp emulado).
       (Nota 2026-11-12: el device lost de la mayoría de los sp era el temporal de `COND_INF` (D-030); con D-030 y
       D-031 vuelven a pasar 8 de ellos, `tools/gpu_skip_tests.sh`.)
       NVIDIA/AMD sólo compilan y enlazan.
6. **Pasos 3–4 (2026-11-03):**
   - (a) `eigenvals` con raíz repetida: si `n2 <= n1/huge` va al límite sin dividir. Mismos bits; hay
     test con trampas de división por cero, overflow e invalid.
   - (b) `tools/audit_kernels.sh` busca símbolos de runtime host en las imágenes de dispositivo:
     PTX/GCN con `-save-temps`, y `llvm-nm` en spir64. `check.sh` lo corre siempre, incluido un build
     Release dp nvptx (estricto). Reemplazos sin cambio de bits:
     - `reshape` por `VEC9_TO_MAT33`/`MAT33_TO_VEC9`;
     - el frente de `lu_solve` y los resultados anidados de `kin_*_tangent` por locales (sin
       `malloc`/`free`);
     - `leadz` por un bucle (lint `fortran-rt`).
   - (c) Permitidos:
     - en Debug, `_gfortran_runtime_error_at` (`-fcheck`) y `_gfortran_matmul_*` (sin inlining en
       -O0);
     - en SPIR64, `for_scale4_v/8_v` (`scale()`): su reemplazo provoca device lost en IGC.
       **Reemplazado por D-030/D-031** (2026-11-12): el device lost venía del temporal de `COND_INF`; `scale()` se
       reemplazó por `fp_scale` y la auditoría ya no lo permite.
7. **Paso 5 (2026-11-03, decisiones del usuario):**
   - (a) Hay una variante `muscle_update_<name>_dev` para memoria de dispositivo pura
     (`omp_target_alloc`, cuFFT). Usa el mismo cuerpo de lote con `has_device_addr`/`is_device_ptr` en
     lugar de `map(present)`, y la API actual no cambia.
   - (b) Cada material tiene su propio tipo de estado visco, con NB = `visco_branches`, en lugar de
     `h(6, MAX_VISCO_BRANCHES)`. Mismos bits.
   - (c) También entran:
     - TL/UL sin la polar, cuando el kernel no usa R (mismos bits en las salidas usadas, o se reporta);
     - medir registros/spills donde la cadena lo permita (amdgcn: metadatos `.vgpr_count`; spir64:
       dumps de IGC). En nvptx sin `ptxas` de CUDA no se puede.
     - `present`/`storage` (OpenMP 6): sólo se documenta mientras los compiladores no acepten
       `storage`.
   - (d) Etapa A aplicada:
     - `muscle_update_<name>_dev` usa `has_device_addr` en los arrays y `map(present, alloc:)` sólo
       en `mat`, sin `is_device_ptr` (ifx lo depreca en arrays que no son `c_ptr`). En CPU llama a los
       bucles normales.
     - El nombre más largo de material baja de 48 a 45 caracteres.
     - En la Intel corren sólo los casos `_dev` de cutting plane en sp. Perzyna sp y todo dp quedan
       excluidos (IGC, dp emulado).
   - (e) Etapa B aplicada. Cada material visco tiene:
     - `mat_state_<name>_<k>` (`plastic` + `visco`);
     - `mat_vstate_<name>_<k>` (`strain_e`, `h(6, visco_branches)`).
     La familia visco agrega `visco_freeze_branches`/`visco_update_branches` (`h(6, nb)` de forma
     explícita), y los `freeze`/`update` de Maxwell y Kelvin pasan a llamarlas, con los mismos bits.
     Se borran `muscle_material_state` y `visco_plastic_state_<k>`: el estado ya no tiene
     `n_branches`. Esto reemplaza D-018.14b.
   - (f) Etapa C aplicada, sólo TL: UL y log espacial necesitan R para la rotación ΔR de `strain_p`,
     y log material la usa en `kin_to_cauchy`.
     - Hay una función de familia nueva, `kin_update_norot`. TL/UL/SL la implementan con
       `kin_finite_update_norot`, que deja R = I; log material y pequeñas deformaciones la reenvían a
       `kin_update`.
     - Los autovalores salen de `sym_eigenvalues`, que va por el mismo camino que `sym_eigen`, así que
       `status` y el resto del estado dan los mismos bits.
     - El kernel usa `kin_update_norot` salvo que `material_kin_rotation(mat)` sea verdadero.
   - (g) Aclaración de D-018.16b (decisión del usuario): los argumentos `optional` están prohibidos sólo
     en las rutinas host que lanzan regiones target (API por lotes), donde se resuelven en el host. Los
     procedimientos `@:DEVICE()` llamados desde dispositivo pueden tenerlos: hay 20 hoy y pasan en los
     tres targets. No se agrega lint.
   - (h) Etapa D: `tools/kernel_resources.sh` mide los recursos por kernel de material, en amdgcn a
     partir del asm GCN y en spir64 a partir de los dumps de IGC. Resultados del 2026-11-03:
     - amdgcn: VGPR 48 / AGPR 24 / SGPR 62 en todos los kernels de lote. Las funciones de punto no se
       inlinean y llevan pila de 9.5–27 KB por work item.
     - spir64: spill sp 0–4.2 KB y dp 9.9–37.7 KB, con dp emulado.
     - nvptx no se mide porque falta el `ptxas` de CUDA. La lista de exclusiones de GPU pasa a
       `tools/gpu_skip_tests.sh`.
8. **Paso 6, benchmark (2026-11-03):** `benchmarks/gpu/bench_batch_gpu.fypp`, en sp y con
   `muscle_update_<m>_dev` sobre datos residentes. La tabla está en `benchmarks/README.md`.
   - Sólo cut_vm_bilinear y cut_yld2004_ludwik corren en Release en la UHD 770:
     - cut_vm: 2–3e8 puntos/s, ~33–41× el CPU de 32 hilos;
     - cut_yld2004: ~2.2e6 puntos/s, 0.5–0.8× el CPU de 32 hilos.
   - Las transferencias de 1e6 puntos tardan ~12 veces el kernel de cut_vm: con datos residentes no se
     transfieren por paso (D-008).
   - El chequeo de GPU contra host se hace dentro de 43 ε, `TOL_DEVICE_sp` = 64. Los puntos que en GPU
     hacen una iteración de Newton más o menos (≤ 0.4 %) se cuentan aparte, sin fallar: el benchmark
     no reemplaza a los tests target.
9. **Cierre de la Fase 8 (2026-11-04):** los tests target de materiales y familias comparan cada bloque
   (tensión, tangente, strain_p, ep, strain_e, h, f, gradiente, hessiana) con su propia escala.
   Antes, la escala común era la del tangente (~2e5), que en sp dejaba sin comprobar la tensión y la
   deformación plástica.
   - Para strain_p y ep la escala es max(|εp|, |ε total|) del punto (decisión del usuario), porque
     εp = ε − C⁻¹:σ cancela. En la GPU Intel daba 84–113 ε contra 64 ε en un punto apenas plástico.

## D-020 — Fase 9: legibilidad y mantenibilidad («FYPP mínimo»)
**Estado:** Aceptada (2026-11-04, decisión del usuario)
**Contexto:** al cerrar la Fase 8, el usuario observó que se abusó de FYPP y se perdió la esencia de
Fortran: macros como `MODEL_FUNCTION`, `MODEL_PARAMS_TYPE` o `MODEL_INIT` generan las cabeceras de las
funciones, los tipos y sus cierres a partir de metadatos de Python (`muscle_meta.py`), de modo que una
plantilla de 48 líneas no muestra el tipo, las firmas, los `intent` ni los `end ...`. Para entender un
modelo hay que leer tres archivos, y los errores de compilación apuntan a líneas que el autor nunca
escribió. Es sobreingeniería respecto de lo que se necesitaba (sp/dp, 2D/3D y despacho estático).
**Decisión:**
1. FYPP se usa sólo para lo que Fortran no puede hacer: instanciar por kind y dimensión,
   interruptores de build, código idéntico sin lógica (`@:DEVICE()`) y la generación de un módulo por
   material del registro. Tipos, firmas, declaraciones, cuerpos y cierres se escriben en Fortran.
2. El criterio está en `.claude/rules/fypp-templates.md` («FYPP mínimo»). El código nuevo lo cumple desde
   ya. El existente se reescribe en la Fase 9.
3. Los metadatos van en la dirección contraria: el Fortran es la fuente y el registro de Python se
   comprueba contra lo que el Fortran declara (test), y se usa sólo donde hace falta generar algo.
4. **La Fase 9 pasa a ser sólo de legibilidad y mantenibilidad**, sin funcionalidad nueva. Los
   pendientes de las Fases 7–8 (FFT, ejecución en NVIDIA/AMD, spill, 2D, etc.) quedan como «Pendientes
   fuera del plan» en PLAN.md, a pedido del usuario.
5. Se itera con el usuario: se parte por un piloto (`muscle_hard_swift`), se revisa cómo quedó y se
   ajusta el criterio antes de seguir con las demás familias.
6. Red de seguridad: los goldens (mismos bits en dp) y los tests target. Un cambio de bits es una
   decisión aparte.
**Consecuencias:** se pierde parte de la unificación automática (por ejemplo, el orden de props de
`mat_init` deja de generarse de una tabla única y pasa a comprobarse con un test). A cambio, un modelo
se escribe y se lee en un solo archivo.

**Enmienda (2026-11-04, tras revisar los puntos estrictos con el usuario):**
- La regla se afloja en cuatro puntos, que están en `.claude/rules/fypp-templates.md`:
  1. se permite el texto fijo con sustitución de nombres, aunque contenga funciones completas
     (`TENSOR_COMMON_PROCEDURES`), y se prohíbe armar declaraciones a partir de metadatos o lógica
     (`MODEL_FUNCTION`): el criterio es dónde vive la información;
  2. se permiten expresiones como argumento de macro (`FD_*`), porque reemplazan a los closures;
     las sentencias y bloques siguen prohibidos;
  3. se permite un nivel de `#:include` para el cuerpo de un solver compartido por materiales, porque
     Fortran no instancia un cuerpo para tipos distintos sin plantillas y el despacho dinámico está
     prohibido en GPU;
  4. no se exige un test que lea Fortran desde Python: el compilador verifica los nombres que asigna
     `mat_init`, y «el Fortran es la fuente» pasa a «sin información duplicada que el compilador no
     pueda verificar».
- El alcance de la Fase 9 incluye tensores, operadores y math (decisión del usuario), con el criterio
  laxo: lo que sea texto fijo con sustitución de nombres se conserva.

**Registro de la ejecución (2026-11-04, noche; el usuario delegó las decisiones menores):**
- **Piloto aprobado por el usuario** (hard swift/bilinear + `muscle_hard_laws`): tipo, firmas, cuerpos y
  cierres en Fortran, `#:for k in REAL_KINDS` por bloque (tipos, `public`, funciones), `@:DEVICE()` justo
  después de `implicit none`, `use` e `interface` de la familia escritos a mano con un `#:for k`. Agregar
  una ley = añadir su nombre a los `use` y a las interfaces de la familia (el compilador lo verifica).
- **Red de seguridad** (`tools/gen_snapshot.sh`): el `.f90` generado, normalizado (sin `#line`,
  comentarios, espacios, continuaciones; listas de nombres una por línea), debe ser igual antes y después
  salvo el orden de las líneas (REORDERED, aceptado: el `public` queda agrupado por kind); la tabla de
  símbolos `nm` idéntica. Más los goldens dp bit a bit y el `check.sh` completo.
- **Hecho (un commit por unidad):** hard, elas, vp/os, visco, yield (incluida la evaluación conjunta a
  mano), kin, y los cinco cuerpos de solver sin `SOLVE_DUMMIES` (firma de `SOLVE` escrita en cada cuerpo).
  Poda: borradas `muscle_models_common.fypp`, todas las macros de contrato y 13 funciones de
  `muscle_meta.py`; los modelos de prueba de los tests se reescribieron a mano. Cifras
  (`tools/fypp_metrics.sh`): llamadas a macros prohibidas 207 → 0; `include/` 1190 → 769 líneas;
  `muscle_meta.py` 2467 → 2287 líneas (111 → 98 funciones).
- **Se conserva, por la regla laxa:** `muscle_meta.py` como datos (`MODEL_PARAMS`, `MODEL_INITS`,
  registro de materiales, `MATERIAL_SOLVERS`, `material_*`, `model_*` de nombres) para generar el módulo y
  la tabla de props de cada material; `ModelFamily`/`ModelFunc` los usan sólo `tests/gpu/target_cases.py`
  y `test_family_target`. Los nombres de `MODEL_PARAMS` se comprueban contra el tipo escrito a mano
  (`tests/tools/test_muscle_model_params.py`).
- **Deuda anotada:** `MAX_VISCO_BRANCHES = 10` está como literal en `muscle_visco_types.fypp` y en
  `muscle_meta` (lo usa `muscle_materials`); campos sin consumidor de `ModelFunc` (`deriv`, `default`) y
  `ModelFamily` (`fd_scale`, `fd_lower`); `model_params` sólo la usan los tests.
- **Corregido en el camino:** `muscle_yield_druckerprager` declaraba las tolerancias `sp` sin guarda de
  kind (rompía el build `-DMUSCLE_REAL_KINDS=dp` en Release); lo detectó `check.sh`.
- **Intermitente de IGC:** SIGSEGV del JIT de la Intel en un kernel de yld2004 una vez durante `check.sh`
  (pasa 3/3 al repetirlo); en PROBLEMAS.md.
- **Revisión final de la Fase 9 (migration-reviewer, 2026-11-04):** aprobada con cambios menores, sin
  regresión numérica ni de build (181/181 en `sp;dp` y 149/149 con `-DMUSCLE_REAL_KINDS=dp`).
  - **Materiales (E5): se quedan como están.** `src/materials/muscle_material*.fypp` y los helpers
    `material_*` de `muscle_meta.py` generan un módulo por entrada del registro (tipo del material,
    tabla de props, APIs de punto y lote): es la generación que la regla permite porque Fortran no puede
    declarar un módulo por entrada. Los candidatos a mejorar (queda como deuda, a pedido del usuario):
    `material_init_call`, `material_prop_check` y `material_domain_checks` arman sentencias Fortran en
    Python.
  - **Criterio laxo, deuda explícita (E6):** `TENSOR_COMMON_TBP` (arma resultado/forma de `vals` desde
    listas y `layout()`), `OP_USES`/`OP_INTERFACES` y los `init2_*_tbp` (lista de argumentos desde
    `layout()`). Las firmas de los seis módulos de operadores sí se escribieron a mano.
  - Corregidos tras la revisión: include del tangente de closest point a un solo nivel, `muscle_fd` movido
    a `tests/common` (sólo lo usan tests), y el test de `MODEL_PARAMS` ahora compara también tipo y
    valor por defecto, las listas de las familias y `MAX_VISCO_BRANCHES`.
- **Investigación de la alternativa D (tipos vacíos en lugar de condicionales de variante), 2026-11-04:**
  aplicada de verdad a VISCO y a KIN en worktrees; goldens dp idénticos, 181/181, sin costo de compilación,
  +8.6 % de librería; pero −34 % de `#:if` a cambio de +865 líneas, 64 específicos por pareja generados
  desde metadatos (contrario a D-020) y funcionalidad nueva (kin + visco). **No se adopta por ahora**; el
  parche queda en `build-snapshots/d_visco_kin.patch`. Análisis completo y plan de mejora de la
  comprensión: `docs/fypp/ANALISIS_COMPRENSION.md`.

## D-021 — Goldens con tolerancia por categoría (2026-11-08)
**Contexto:** los goldens (`tests/reference/`) se capturaron para migrar el código OO a FYPP bit a bit. Con
las Fases 0–9 cerradas ese objetivo se cumplió; la comparación bit a bit fallaba con otros flags
(`Profile_unittests`: `kin_golden` 4e-16, `math_golden_derivatives` 2.4e-13) y bloquearía optimizaciones
legítimas. Decisión del usuario: los goldens deben comprobar la funcionalidad, no los bits.
**Decisión:**
1. `golden_check` compara con tolerancia relativa por categoría, respecto del mayor valor real del
   registro: álgebra directa 1e-13 (por defecto), diferencias finitas 1e-9, solvers iterativos 1e-9
   (10 × `NEWTON_TOL_dp`); en los registros de solver, `status` exacto e iteraciones ±1. Un NaN coincide con
   cualquier NaN; un NaN de un solo lado nunca está dentro de la tolerancia.
2. Modo estricto opcional `MUSCLE_GOLDEN_STRICT=1` (bit a bit, sólo con sentido en gfortran Debug) y modo
   informe `MUSCLE_GOLDEN_REPORT=1`.
3. El nombre «golden» se mantiene por ahora; el renombre a «reference» y el reordenamiento de los tests
   quedan para más adelante (ESTADO_OPTIMO, 3.8).
4. Tests de mutación: `*_golden_mutation` perturba un parámetro en 1e-6 relativo y debe fallar;
   `*_golden_mutation_strict` perturba 1 ulp con el modo estricto. (Antes: 1 ulp en modo normal.)
5. El driver de los goldens de operadores (`ops_golden_cases.fypp`) usa la rutina central en lugar de su
   copia propia de la lectura y la comparación.
**Calibración (VERIFICADO, máximos medidos):** álgebra 2.0e-14 (tangente corotacional, cpb06ex2),
diferencias finitas 3.1e-11, solver 4.7e-13, tangente numérica de solver 1.2e-12, status/iteraciones sin
diferencias, en `Profile_unittests` y Release (-O3). Debug, Release y `Profile_unittests`: 186/186 tests;
Debug con `MUSCLE_GOLDEN_STRICT=1`: los 23 tests de goldens pasan bit a bit.
**Fuera de la configuración acordada:** con FMA (`-march=native -O3`) fallan tests por razones que la
tolerancia no debe tapar (PROBLEMAS, «Tests»): registros cuyo valor exacto es 0 (sin escala), cantidades no
únicas (autovectores con autovalores repetidos), casos mal condicionados y la tangente numérica en un punto
no suave. Quedan documentados para rediseñar esos tests.
- **Repetición aceptada y explicitud (decisión del usuario, 2026-11-08):** se mantienen escritas a mano, aunque se
  repitan, la evaluación conjunta por defecto de los criterios de fluencia (5 copias de 36 líneas que sólo
  cambian el nombre del modelo) y la firma de `SOLVE` con su preámbulo en cada cuerpo de solver: cada archivo
  se lee completo. Los módulos de familia siguen listando los modelos de forma explícita (`use`, `public`,
  `interface`), sin generarlos con `#:for`; la prueba del recién llegado mostró que agregar un modelo así es fácil.
- **Hallazgo tardío de la Fase 9 (2026-11-08):** `muscle_kin_finite.fypp` seguía generando 24 funciones (medidas
  de deformación, conversiones, tangentes, rotaciones) a partir de tablas de Python con expresiones Fortran como
  datos, el patrón que D-020 prohíbe; no se detectó porque no usaba las macros `MODEL_*` que medía
  `tools/fypp_metrics.sh`. Se reescribió a mano (929 → 1288 líneas; sólo quedan como `#:set` los umbrales por
  kind `EIG_MIN` y `LOG_SERIES`, que son datos). Lección: la métrica de «macros prohibidas» no cubre la lógica
  inline (`#:for` sobre tablas con código); revisar a ojo los `#:set` de las plantillas.

## D-022 — Todos los inicializadores de modelos devuelven `status` (2026-11-08)
**Contexto:** el inventario de convenciones (ESTADO_OPTIMO, 3.6) encontró que los inicializadores de visco
devolvían `status` y los de elasticidad y fluencia no. **Decisión del usuario:** agregar `status` a todos.
**Hecho:** `elas_linear_init`, `yield_cpb06_init`, `yield_cpb06ex2_init`, `yield_druckerprager_init` y
`yield_yld2004_init` reciben `integer, intent(out) :: status` como último argumento (hoy siempre `STATUS_OK`:
no validan nada; las validaciones quedan para la etapa de robustez). `MODEL_INITS` los marca con status, así que
`mat_init_<name>` lo pasa y lo comprueba como con visco. Cambia la API pública: quien llame a estos cinco
inicializadores debe pasar `status`. Sin cambio numérico (goldens iguales); los módulos generados cambian sólo en
la firma, la asignación, el `use`, la llamada con `status=status` y la comprobación en `mat_init`.

## D-023 — Tests unitarios en vez de referencias; criterio de mutaciones (plan R, 2026-11-09)
**Contexto:** los 22 tests «golden» (~8 600 registros) se capturaron para migrar bit a bit; sólo detectan
cambios, no dicen si el resultado es correcto. **Decisiones del usuario:** (1) se acepta un mínimo de
referencias donde no hay verdad independiente; (2) **criterio de mutaciones**: una referencia sólo se borra si
los tests unitarios detectan todas las mutaciones del código que ella detecta; (3) se hace todo el plan R
(ESTADO_OPTIMO, Bloque R); (4) se renombra «golden» → «reference».
**Herramienta:** `tools/mutation_check.py <area>` (áreas en `tests/mutations/<area>.toml`): aplica una a una
mutaciones de las fuentes (operador binario, literal real, subíndice, relacional; muestreadas con semilla o
todas, más las escritas a mano), recompila un build dp Debug y reporta qué grupo las detecta. Se aplicó con
**todas** las mutaciones de cada área (≈ 3 600), no con una muestra.
**Resultado:** quedan 3 tests de referencia con 83 registros: `math_reference` (6: autovectores elegidos en
empates exactos, elección no especificada; la rama `|1 - d| <= EPS` de `eigenvals`, PROBLEMAS) y las secuencias
reducidas de `closest_point_time` (50) y `closest_point_perzyna` (27: amortiguamiento, búsqueda lineal y
función de mérito del Newton, sin verdad independiente). Todo lo demás lo cubren tests unitarios nuevos o
ampliados (ESTADO_OPTIMO R1–R6). API renombrada: `muscle_test_reference` (`reference_init/check/finish`,
`REFERENCE_TOL_*`), `MUSCLE_REFERENCE_STRICT/REPORT`, `*_reference_mutation`, `muscle_reference_test` en CMake;
`mat_sequence_common.fypp`, `case_params.py` (`CASE_PARAMS`), `solver_case`, `muscle_test_solver_sequence`.
**Consecuencias:** un modelo o solver nuevo lleva tests unitarios y el chequeo de mutaciones de su área, no un
golden; un material nuevo del registro de prueba entra al test de secuencia de su solver con `solver_case(0,
...)`. Etiquetas de ctest: `kinds`, `legacy`, `reference`, `target`. Hallazgos registrados en PROBLEMAS:
`check_approx` y `.approx.` aceptaban Inf (corregidos), y quedan marcados sin corregir la rama de `eigenvals`, el jacobiano de
`closest_point_2` sin factor 2 de corte y la no convergencia de `closest_point_time` + JC con poca sobretensión.

## D-024 — Se retira el pack/unpack genérico de `muscle_plastic_state` (2026-11-09)

**Contexto:** `plastic_state_pack/unpack` (`PLASTIC_STATE_SIZE = 13`: tensión, deformación plástica, εp) era la
traducción de `Plastic_material_history%pack_to_fea/unpack_from_fea` (Fase 5). Desde D-018.7 cada material
tiene su propio pack (`mat_pack_<name>`, `MAT_<NAME>_NSTATE`: deformación plástica, εp y el estado viscoelástico,
sin la tensión, que va en `sig`), que es el contrato documentado en `pages/integration.md` y `pages/migration.md`.
El genérico quedó sin uso (sólo su test) y su documentación lo presentaba como «el» layout del array de
historia FE (re-análisis de comprensión, `REANALISIS_COMPRENSION.md`).
**Decisión del usuario:** opción (a): se borran `plastic_state_pack/unpack` y `PLASTIC_STATE_SIZE`; se conserva
el tipo `plastic_state_<k>` (estado en memoria de todos los kernels). No cambia ningún material ni el layout de
`hsv`. El test `plasticity-plastic_state_kinds` comprueba ahora los valores por defecto y la copia en un
procedimiento `@:DEVICE()`; el layout de los materiales lo comprueban los tests de `tests/materials/`.

## D-025 — Se retiran los 38 tests heredados `.f90`; nombres de tests espejo de `src` (2026-11-10)

**Contexto:** re-análisis de la calidad de los tests (`REANALISIS_TESTS.md`): con `mutation_check --samples all`
en todas las áreas (≈ 4 300 mutaciones, los heredados como grupo «referencia»), los heredados sólo detectaban 14
mutaciones que los demás tests no detectaban (normalización de CPB06, C′ de CPB06ex2, `kin_pull_back_tangent`);
los tests unitarios agregados en el commit `6723289` las detectan todas. Los heredados eran 9 650 líneas con 49 %
en bloques repetidos y la regla «no se mueven ni se reescriben» de la migración.
**Decisión del usuario:** (1) se retiran los 38 tests de la etiqueta `legacy` (tipos de tensores 17, incluidos los
2 de OpenMP; operadores 5; fluencia 6 y su biblioteca `yield_test_numeric`; endurecimiento 4; elasticidad 2,
incluido el de OpenMP; math 2; viscoplásticas 1; cinemáticas 1), criterio de D-023. El uso en hilos lo cubren
`materials-mat_batch_openmp`, los tests target y el lint/auditoría de kernels (sin `save`). Desaparece la etiqueta
`legacy` (quedan `kinds`, `reference`, `target`). (2) Nombres: los tests de `src/<ruta>/muscle_<x>.fypp` son
`tests/<ruta>/test_muscle_<x>[_<tema>].fypp` (siempre con `muscle_`); un test de varios módulos lleva el de su
tema (`test_muscle_material_<solver>`, `test_muscle_fypp_smoke`). Se renombraron 44 archivos con `git mv`; los
nombres de ctest (`<carpeta>-<nombre>`) no cambian, así que los filtros `-R` y las áreas de `mutation_check`
siguen valiendo.

## D-026 — Directivas FYPP indentadas al nivel del código que controlan (2026-11-11)

**Contexto:** las ~3 000 directivas de control de las plantillas estaban en la columna 0, lo que corta la lectura
del Fortran (re-análisis de convenciones, tarea V9). FYPP 3.1 acepta espacios antes de una directiva y borra la
línea entera al preprocesar, así que su sangría no llega al código generado.
**Decisión del usuario:** una directiva de control (`#:if`, `#:for`, `#:set`, `#:assert`, `#:call`, `#:block`…)
toma la sangría de la línea de código que le sigue; `#:elif`, `#:else` y el cierre, la de su apertura. El
contenido **no** gana un nivel extra (opción A): la sangría sigue mostrando el anidamiento del Fortran, el
generado no cambia y las líneas no se alargan. Quedan en columna 0 `#:include`, `#:def`/`#:enddef`, la guarda más
externa de un archivo entero (p. ej. `#:if defined('MAT')` de un cuerpo incluido) y los comentarios `#!`.
Aplicado a las 119 plantillas de `src/`, `include/`, `tests/` y `benchmarks/` con un script de una sola vez
(3 164 líneas, sólo espacios iniciales); `gen_snapshot.sh` idéntico y el Fortran generado de los tests idéntico.
Regla en `pages/conventions.md` y `.claude/rules/fypp-templates.md`.

## D-027 — Módulos de las familias `os`, `elas` y `yield` con el prefijo de sus funciones (2026-11-11)

**Contexto:** re-análisis de convenciones (V6, I6): tres familias no seguían la regla «la familia `x` vive en
`muscle_x_<modelo>` y `muscle_x_laws`», y eso confundió a la prueba del recién llegado (ley Perić).
**Decisión del usuario:** renombrar archivos y módulos, **sin alias** de compatibilidad (la API FYPP no se publicó:
vive en esta rama sin integrar; `main` es el código OO anterior):

| Antes | Ahora |
|---|---|
| `muscle_vp_overstress_laws` | `muscle_os_laws` |
| `muscle_vp_perzyna`, `muscle_vp_chaboche` | `muscle_os_perzyna`, `muscle_os_chaboche` |
| `muscle_elasticity_laws`, `muscle_elasticity_linear` | `muscle_elas_laws`, `muscle_elas_linear` |
| `muscle_visco_frozen` | `muscle_elas_frozen` (sigue en `src/viscoelasticity/`) |
| `muscle_yield_criteria` | `muscle_yield_laws` |
| `muscle_anisotropic_yield` | `muscle_yield_anisotropic` |

Las carpetas no cambian. Desaparece el caso especial `MODEL_MODULES` de `muscle_meta.py` (el módulo de un modelo es
siempre `ModelFamily.module`). Los tests siguen D-025 (`test_muscle_elas_laws_kinds`, `test_muscle_yield_laws_kinds`,
`test_muscle_os_chaboche`); los nombres de ctest no cambian. Verificado: el Fortran generado y la tabla de símbolos son
idénticos a los de antes con los nombres cambiados. Los documentos anteriores de `docs/fypp/` conservan los nombres de
su fecha.

## D-028 — Autoproyecciones robustas en math; rotación de Jacobi sin desborde (2026-11-11)

> **Reemplazada en parte por D-032 (2026-11-12):** `eigenprojections`/`eigenvals_refined` siguen siendo los genéricos
> públicos, pero sobre `spectral_decomposition`; `jacobi_band`, `jacobi_eigen` y `JACOBI_GAP_LO/HI_<k>` se retiraron. La
> tangente sin desborde de (2) sigue viva en la rotación plana de `spectral_decomposition`.

**Contexto:** revisión de las auxiliares privadas de `muscle_kin_finite` (la de `is_identity`, D-011): las rutinas
de autovalores y autoproyecciones de un tensor simétrico (`sym_eigen`, `sym_eigenvalues`, `jacobi_band`,
`jacobi_eigen`, D-018.19a) son álgebra general, no cinemática.
**Decisión del usuario:** (1) se mueven a `muscle_math_spectral_derivs` como genéricos públicos `eigenprojections(s,
eig, proj)` y `eigenvals_refined(s, eig)` (las de Jacobi quedan privadas); las constantes `KIN_JACOBI_GAP_LO/HI_<k>` de
`muscle_kinds` pasan a `JACOBI_GAP_LO/HI_<k>`, sin alias; la copia de Jacobi arma la matriz con la asignación
`ten_3D2O = ten_3D2Osym` y `VEC9_TO_MAT33` (mismos valores que el antiguo `mat33`). Test nuevo
`math-math_spectral_projections` con verdad en precisión cuádruple (Q diag(λ) Qᵀ conocidos): autovalores, cada
proyección dentro de eps/gap (la forma cerrada, eps/gap², no lo cumpliría), sumas, ortogonalidad, pares cercanos y
repetidos arriba, abajo y al medio, dos rotaciones y escalas 1e3 y 1e-3. (2) **Corrección** encontrada por ese test:
la rotación de Jacobi calculaba `theta = (m(q,q) - m(p,p))/(2 m(p,q))` y sólo saltaba `|m(p,q)| <= tiny`; en sp, con
autovalores ~1e3, `m(p,q)` cerca de `tiny` hace desbordar theta (excepción de punto flotante; Inf sin trampa). Ahora,
para theta > 1e15, `t = m(p,q)/(m(q,q) - m(p,p))` sin formar theta (el mismo valor que el `0.5/theta` anterior, sin
desborde). Cambia sólo el último bit de t en esas rotaciones (ángulos ~1e-15); los tests de cinemática (verdad en
precisión cuádruple) siguen en verde. En la cinemática en sp el desborde podía darse con estiramientos grandes.
**Ampliación (2026-11-11, misma revisión):** también salen de `muscle_kin_finite` a `muscle_math_operations`
`polar_rotation(F)` (rotación de la descomposición polar F = R U por Newton escalado; antes sólo existía en sp, ahora
en todos los kinds) y `log_divided_difference(a, b)` = (ln a − ln b)/(a − b) con su serie cerca de a = b (antes
`log_ratio`; el umbral `LOG_SERIES` por kind se mueve con ella). Test `math-math_polar_log` con verdad en precisión
cuádruple. Resultados de la cinemática sin cambios (mismo código).
**Ampliación 2 (2026-11-11, opción B del usuario):** las conversiones entre tensor y matriz 3x3 de `muscle_kin_finite`
(`mat33`, `voigt6`, los índices escritos a mano del triángulo superior en `kin_pk1_to_pk2`/`kin_pk1_to_cauchy`) y la
macro `VEC9_TO_MAT33` pasan a funciones libres públicas de los tipos de tensores (`muscle_tensor_functions`):
`to_matrix(a)` (`2Osym` y `2O`, 2D y 3D; en 2D los acoplamientos fuera del plano valen 0) y `from_upper(m)` (sólo
`ten_3D2Osym`, genérico en kind: en 2D no se distinguiría por el argumento y no hace falta hoy). `from_upper` devuelve un
tensor, no un arreglo de 6 valores, para no exponer el orden Voigt a quien lo usa; lee sólo el triángulo superior (la
parte simétrica (m + mᵀ)/2 es otra operación). Son copias: código generado distinto sólo en `kin_finite`, las dos de
math y los tipos 2O/2Osym, y resultados iguales bit a bit (tests y referencias estrictas). `VEC9_TO_MAT33` se borró
(sin usos); `MAT33_TO_VEC9` queda (un 2O no tiene `from_matrix`). Tests en `tensors-tensor_functions` (accesores con
nombre en su entrada de la matriz, ceros en 2D, `square`/`transpose` contra `matmul`/`transpose` de la matriz,
`from_upper(to_matrix(a)) = a`, triángulo inferior ignorado) y en los tests `target` de tensores; mutaciones: 2/2 de la
herramienta y 5/5 de índices hechas a mano (los índices los genera FYPP y la herramienta no los muta).

## D-029 — `src/math/` ordenada por tema y por capa (2026-11-11)

**Contexto:** `REANALISIS_MATH.md`: los nombres de los módulos de math no decían lo que contenían (`operations` era
un cajón; `spectral_derivs` tenía también la descomposición, Jacobi y una utilidad de bits; `derivatives` sólo
`dsquare`), lo espectral estaba partido en dos módulos y las utilidades de dispositivo vivían en módulos ajenos.
El usuario pidió hacer la opción C (unificar los algoritmos espectrales, `PLAN_ESPECTRAL.md`) con el criterio de que
el kernel no crezca mucho y de que, **ante todo, la carpeta quede ordenada y se entienda a primera vista**; lo que
no se pueda unificar barato queda separado pero ordenado. Por eso la fase B (orden) va antes que la unificación.
**Decisión del usuario (fase B, sin alias, como D-027; resultados idénticos bit a bit):**
- `muscle_math_fp` (nuevo): `QUIET_NAN` (de `inverses`) y `fp_exponent` (de `spectral_derivs`), reemplazos de
  IEEE/intrínsecos para el dispositivo; bajo los tensores.
- `muscle_math_eigen` (nuevo): `eigenvals` (de `operations`), `spectral_decomposition`, `eigenprojections`,
  `eigenvals_refined` y el Jacobi privado (de `spectral_derivs`).
- `muscle_math_spectral_derivs`: derivadas (`dEigenvalues_dTensor`, `d2Eigenvalues_dTensor2`, `spectral_gradient`,
  `spectral_hessian`) y `dsquare` (se retira `muscle_math_derivatives`).
- `muscle_math_operations` → `muscle_math_matrix_functions` (`polar_rotation`, `log_divided_difference`); el nombre
  `muscle_math_tensor_functions` se descartó por parecerse a `muscle_tensor_functions` (pregunta al usuario).
- Include `muscle_math_lu_inverse.fypp` → `muscle_math_lu_bodies.fypp` (tiene también el cuerpo del solver).
- Estilo: palabras clave en minúscula en `inverses` (los nombres `M33INV`, `AINV`, `OK_FLAG`, `KAPPA` se conservan),
  sangría y `<=` en `eigenvals`.
- Tests espejo: `math-math_inverses`, `math-math_fd_known_functions` (las diferencias finitas sobre funciones
  conocidas; el test de `dsquare` pasa a `math-math_spectral_cases`, ahora en sp y dp), `math-math_eigen_eigenvals`,
  `math-math_eigen_projections`, `math-math_matrix_functions`.
- `pages/math.md`: las dos capas (bajo y sobre los tensores) y la tabla «necesito X → módulo Y».
`muscle_math_eigen` usa `dEigenvalues_dTensor` (las autoproyecciones de Cayley-Hamilton de `eigenprojections`):
la dependencia eigen → spectral_derivs desaparece si la fase F2 de `PLAN_ESPECTRAL.md` unifica la cinemática.
Verificación por paso: código generado normalizado igual línea por línea salvo el andamiaje de los módulos
(`gen_snapshot`), ctest 133/133, referencias estrictas, mutaciones de `dsquare` 45/45 desde su nuevo test.

## D-030 — Sin temporales de tamaño variable en código de dispositivo; `NORM_INF` (2026-11-12)

**Contexto:** al reemplazar `scale()` por `fp_scale` (`PLAN_ESPECTRAL.md`, problema de IGC con la emulación FP64,
D-031) tres kernels sp de material (`tm_jc_cpb06_ludwik`, `pz_cha_yld2004_ludwik`,
`pz_cpb06ex2_bilinear`) perdían el dispositivo, como en el intento de D-019.6. Diagnóstico (build sólo sp, sin dp
ni emulación en la imagen): el fallo seguía con 1 work item, `-ze-opt-disable` y la fusión de EUs apagada;
`IGC_StackOverflowDetection=1` informa «Stack overflow detected!» y con más memoria privada por work item pasa.
Causa: `COND_INF` calculaba `maxval(sum(abs(A), dim=2))` sobre un `A(:,:)` de forma asumida; en el dispositivo el
temporal de `sum(..., dim=2)` es un `alloca` de tamaño dinámico (el único de la imagen) que IGC no puede acotar. Con
stack calls en el kernel (las de la emulación FP64 de `scale()`) IGC reserva 1,4 MB por hilo y alcanza; sin ellas
reserva 4 KB por work item y el kernel la desborda. La misma expresión estaba escrita 6 veces (`COND_INF`,
`M33INV`, `M44INV`, el cuerpo del solver LU con `KAPPA`).
**Decisión del usuario:**
- `NORM_INF_<k>(A)` en `muscle_math_inverses` (genérico `NORM_INF`): \( \|A\|_\infty \) con bucles explícitos,
  cada fila sumada de izquierda a derecha (el orden de `sum(abs(A), dim=2)`: mismos bits en CPU sin NaN), sin
  temporales; NaN en cuanto una suma de fila es NaN. Esto corrige un error latente: `maxval` de gfortran salta las
  filas NaN, así que `M33INV`/`M44INV` aceptaban la inversa adjunta con un NaN en una fila (κ finito) en lugar de
  recurrir a LU (`PROBLEMAS.md`, test `test_adjugate_nan`); `COND_INF` da NaN si una fila de `A` o de su inversa lo
  es. Reemplaza las 6 apariciones; `muscle_math_lu_solve` la importa (dependencia nueva dentro de la capa baja de
  math, sin ciclo).
- Regla (`gpu-static-dispatch.md`): en código de dispositivo no hay arrays temporales de tamaño conocido sólo en
  ejecución (reducciones con `dim=` o expresiones de arrays sobre argumentos de forma asumida): bucles explícitos.
  Una comprobación automática (buscar `alloca` dinámicos en las imágenes SPIR64 de la auditoría) queda propuesta.
**Verificación:** ctest 133/133 (referencias incluidas: mismos bits); tests nuevos `test_norm_inf` y
`test_adjugate_nan` (falla con la norma anterior) en
`math-math_inverses` (valores exactos a mano, fila mayor primera/última/al medio, empates, rectangular, sin filas, NaN
e Inf en cada fila, y bit a bit contra `maxval(sum(abs(A), dim=2))` en matrices aleatorias); mutaciones de
`NORM_INF` 5/6 (la sobreviviente, `<=` → `<`, es equivalente: en un empate asigna el mismo valor). En la GPU Intel,
con `fp_scale` (D-031) los tres kernels pasan, y con las dos correcciones juntas vuelven a pasar 8 kernels sp
excluidos (closest point, closest point 2, closest point time y Perzyna; no se midió cuántos por cada una).

## D-031 — `fp_scale` en lugar de `scale()` en el código de dispositivo (2026-11-12)

**Contexto:** los candidatos de Jacobi de `PLAN_ESPECTRAL.md` (fase F0) daban resultados erróneos en la GPU Intel.
Investigación (reproductores e informe en `docs/fypp/igc_mre/`): ifx compila `scale()` sobre `real32` con
operaciones dp (`for_scale4_v`), así que un kernel sp llama a la emulación FP64 de IGC; con esas llamadas dentro de
una rama divergente y arreglos locales indexados con tablas `PARAMETER`, cuando sólo uno de los dos hilos SIMD8 de
una pareja de EU fusionadas ejecuta la rama, el otro devuelve valores erróneos. El mismo código máquina es correcto
con la fusión de EUs desactivada; se reproduce con lo último publicado (IGC 2.42.6, NEO 26.40). D-019.6 había
mantenido `scale()` porque su reemplazo hacía perder el dispositivo; la causa de eso era el temporal de `COND_INF`
(D-030), no el reemplazo.
**Decisión del usuario:**
- `fp_scale_<k>(x, n)` en `muscle_math_fp` (genérico `fp_scale`, `elemental`): \( x\,2^n \) con operaciones enteras
  sobre los bits IEEE, redondeo al par más cercano; bit a bit el resultado IEEE (subnormales, desborde a Inf,
  subdesbordamiento total, ceros, Inf y NaN sin cambio). No usa dp ni depende de `fp-model` (una versión con
  productos por potencias de 2, al estilo de `scalbnf` de musl, cambia bits con la reasociación de `-fp-model fast`
  de ifx en los casos extremos).
- Reemplaza las 20 llamadas a `scale()` de `muscle_math_eigen` y de los criterios cpb06, cpb06ex2 y yld2004.
- Lint `scale` en `tools/check.sh` (como `exponent`): `scale()` prohibido en el código de kernels.
- La auditoría SPIR64 (`tools/audit_kernels.sh`) ya no permite `for_scale*_v`.
- Los reproductores y el informe del error de IGC quedan en `docs/fypp/igc_mre/` (para reportarlo a Intel y volver a
  probar con drivers nuevos).
**Verificación:** `fp_scale` = `scale` bit a bit en 176 millones de casos en el host (sp y dp; gfortran, e ifx con
`-no-ftz`: ifx a -O2 lleva los subnormales a cero); test `fp_scale = scale` en `math-math_kernels_target` (valores
especiales, subnormales con todos los casos de redondeo, potencias que llevan al rango subnormal bit a bit, desborde
y subdesbordamiento total) y mutaciones de `fp_scale` 14/14. ctest 133/133 (mismos bits). En la GPU Intel, los 40
tests target pasan y los 8 kernels sp excluidos de D-030 vuelven a pasar (medido con las dos juntas); los
reproductores dan 0 errores con
`fp_scale` y la imagen ya no necesita FP64.

## D-032 — Una sola descomposición espectral: K (Kopp mejorado) en `spectral_decomposition` (2026-11-12)

**Contexto:** `PLAN_ESPECTRAL.md` (opción C de `REANALISIS_MATH.md`): hay dos estrategias para autovalores
cercanos, la rotación plana de `spectral_decomposition` (criterios anisótropos) y la banda de Jacobi de
`eigenprojections` (cinemática finita). Las dos pierden precisión en casos medidos (C = I + 1e-11 en dp, escalas
extremas, excepciones en sp; `ESPECTRAL_F0.md` §1). Mediciones de F0 y su repetición tras D-030/D-031
(`ESPECTRAL_F0.md` §8): todos los candidatos son correctos en la GPU Intel; K es el único dentro de 1,2x en todo
(GPU 1,14x, CPU 0,87–0,92x en los criterios y 0,31–0,50x en la cinemática) y pasa `math_reference`.
**Decisión del usuario:**
- **Algoritmo K:** desviador escalado por potencia de 2 (exacto, `fp_scale`), autovalores en forma cerrada (Kopp
  2008), el vector del autovalor más separado por producto vectorial, el par restante por una rotación plana sin
  desborde (D-028) y los autovalores finales como cocientes de Rayleigh. Errores de 2 a 4 ε en todas las brechas.
- **Nombre:** K es la implementación de `spectral_decomposition` (misma firma; sin genérico nuevo ni envoltorio).
- **Alcance:** los dos consumidores. `eigenprojections` y `eigenvals_refined` se escriben sobre
  `spectral_decomposition` (\( P_a = v_a \otimes v_a \)); se retiran `jacobi_band`, `jacobi_eigen` y
  `JACOBI_GAP_LO/HI_<k>` de `muscle_kinds` (D-018.19a queda reemplazada).
- **F4 (b):** la f de los criterios anisótropos usa los λ de `spectral_decomposition`, la misma fuente que su
  gradiente y hessiana. Desaparece la doble fuente de autovalores, y con ella el problema «marcado, sin corregir» de
  `eigenvals` para los criterios (autorizado por el usuario). `eigenvals` sigue disponible para quien sólo quiera
  valores.
- **F5:** `dEigenvalues_dTensor` y `d2Eigenvalues_dTensor2` se mantienen como API pública (derivadas por
  Cayley-Hamilton, sin vectores), aunque la librería deje de usarlas.
**Fases (en `PLAN_ESPECTRAL.md`):** F1+F3+F4 en un paso (criterios), F2 (cinemática), F6/F7 (limpieza y cierre).
**F2 hecho (2026-11-12):** `eigenprojections` y `eigenvals_refined` sobre `spectral_decomposition`
(\( P_a = v_a \otimes v_a \)); se retiraron `jacobi_band`, `jacobi_eigen`, `JACOBI_GAP_LO/HI_<k>` y la dependencia
de `muscle_math_eigen` sobre `muscle_math_spectral_derivs`. Cambio numérico de la deformación de Hencky
\( \tfrac12\ln\mathbf{C} \) frente a una verdad de 50 dígitos (400 F aleatorias, |F − I| de 1e-1 a 1e-13; script en el
scratchpad): en dp el error relativo baja de 3 a 100 veces (1e-1: 5.3e-13 → 4.5e-15; 1e-6: 1.3e-9 → 3.9e-10) y cerca
de la identidad deja de perderse el desviador (|F − I| ≤ 1e-10: 120 % → 2e-6 … 3e-3, el límite ~ε/|F − I| de formar
C); en sp baja de 2 a 8 veces hasta 1e-6, y por debajo de ~1e-7 (C = I en sp) da H = 0 en lugar de valores
arbitrarios. R cambia ≤ 5e-11 (dp). La tensión y la tangente de las cinemáticas finitas dependen de R y H y no se
midieron aparte; los tests de cinemática y de materiales pasan con sus tolerancias. El salto en los bordes de la
antigua franja desaparece: el test de autoproyecciones (verdad en cuádruple) cubre brechas a ambos lados (2e-2 y
1e-11, 1e-13 en dp) y casi identidad, y falla con la franja; test nuevo `test_near_identity` de la deformación de
Hencky en `kinematics-kin_finite` (falla con la franja en dp).

## D-033 — Skills de Claude Code versionadas en el repo; `check.sh --doc` (2026-11-13)

**Contexto:** la skill `muscle-nuevo-modulo` vivía en un repo aparte y quedó desactualizada con la migración y con
D-030 a D-032; el usuario pidió además skills para revisar algoritmos y documentación. `tools/check.sh` sólo buscaba
FORD en el `PATH` y nunca generaba la documentación donde FORD está en el entorno conda de `environment.yml`.
**Decisión del usuario:**
- Las skills del proyecto se versionan en `.claude/skills/` (`.gitignore`: `.claude/*` y `!.claude/skills/`); el
  resto de `.claude/` (settings, reglas, agentes) y `CLAUDE.md` siguen locales. Un cambio del proyecto que altere una
  regla, una ruta o una decisión que una skill cita actualiza la skill en el mismo PR (`pages/contributing.md`).
- Skills: `muscle-nuevo-modulo` (crear un módulo), `muscle-revisar-algoritmo` (revisar precisión, tests, costo CPU/GPU
  y kernel de un algoritmo) y `muscle-revisar-documentacion` (documentación al día con el código).
- `tools/check.sh` busca FORD en el `PATH`, en el entorno conda `ford` y en `build-check/.ford-venv`; modo nuevo
  `--doc` (lint + FORD + `doc_drift.py`, ~2 min) antes de un PR que toque documentación o API pública.
