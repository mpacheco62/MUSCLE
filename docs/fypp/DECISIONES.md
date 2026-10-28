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
**Estado:** Propuesta
**Contexto:** Fortran no tiene genéricos de tipo; los PDT (`type :: t(k)`) tienen bugs en
gfortran y soporte incierto en compiladores GPU.
**Decisión:** Tipos concretos `ten_3D2Osym_sp` / `ten_3D2Osym_dp`, procedimientos con sufijo,
genéricos que agrupan ambos.
**Consecuencias:** Más símbolos; resolución de genéricos por kind; sin coste en tiempo de ejecución.

## D-003 — Nombres históricos = doble precisión
**Estado:** Propuesta
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
**Pendiente:** verificar `#line` con ifx/IntelLLVM (`-fpp`) cuando esté instalado.
