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
    comportamiento dp sólo en ese borde).

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
