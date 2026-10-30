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
