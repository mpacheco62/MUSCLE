# Plan de migración de MUSCLE a FYPP

Estado: **borrador v1 (2026-10-25)** · Rama base: `dev` · Línea base: 60/60 tests (2026-10-28, `dev` en `db0bc91`, gfortran 15.2, Debug)

## 1. Objetivos

| # | Objetivo | Criterio de "hecho" |
|---|----------|---------------------|
| O1 | Tensores en precisión simple y doble | Cada tipo `ten_*`/`iden_*` existe en `_sp` y `_dp` desde una sola plantilla; tests para ambos kinds |
| O2 | Mismo código 2D/3D | Ninguna operación tiene implementación 2D y 3D escritas a mano por separado; tests de embedding 2D↔3D |
| O3 | Despacho estático, compilable para GPU | Ruta de kernel (tensores → modelos → solver) sin `class`, `deferred`, `select type`, closures, allocatables, I/O ni LAPACK; compila con gfortran OpenMP offload (`-fopenmp -foffload=nvptx-none` / `amdgcn-amdhsa`) y luego ifx (`-fopenmp-targets=spir64`), y un test corre en GPU NVIDIA, AMD e Intel |
| O4 | Sin regresiones | API pública de doble precisión intacta (alias `ten_3D2Osym` ⇒ `ten_3D2Osym_dp`), todos los tests en verde en cada commit |

## 2. Diagnóstico del código actual (2026-10-25)

**Tensores** (`src/tensors`, ~10.3k líneas)
- 16 módulos de tipo (10 tensores + 6 identidades) + 6 de operadores. Sólo `real(real64)`.
- 2D y 3D duplicados a mano: `ten_2D2Osym`(4) vs `ten_3D2Osym`(6), `ten_2D2O`(5) vs `ten_3D2O`(9),
  `ten_2D4O3sym`(10) vs `ten_3D4O3sym`(21), `ten_2D4O2sym`(4x4) vs `ten_3D4O2sym`(6x6), `ten_2D4O`(5x5) vs `ten_3D4O`(9x9).
- Operadores por pares de tipos: 222 `module procedure` (sólo `addition_subtraction` tiene 1886 líneas).
  Es el caso ideal para generación por tablas.
- Ya son `pure` y sin polimorfismo → buen punto de partida para GPU. Problemas: UDTIO
  `write(formatted)`, `print`/`write` en varios tipos, `inverse` usa `M33INV` (OK, propio).

**Math** (`src/math`)
- `muscle_math_derivatives`: derivadas numéricas con `procedure(f_...) :: func` → los
  llamadores pasan **procedimientos internos** (closures). El linker ya avisa de stack
  ejecutable en `muscle_yield_base.f90.o` por los trampolines. Incompatible con GPU.
- `muscle_math_linsolvers`: LAPACK/MKL con arrays `allocatable` (`real(8)` → `real(dp)` en la Fase 1).
- `muscle_math_lu_inverse.inc`: LU propia vía `include` → reutilizable para GPU (Fase 1: macro FYPP `muscle_math_lu_inverse.fypp`).

**Modelos** (elasticidad, plasticidad, viscoelasticidad, cinemática)
- Todo cuelga de tipos abstractos con bindings `deferred`: `Base_elasticity`,
  `Base_hardening_laws`, `Base_yield_critera`, `Base_viscoplastic_law`,
  `Base_overstress_law`, `Base_viscoelasticity`, `Base_kinematics`.
- `Base_elasticity` duplica cada operación en `_3D`/`_2D` (+ variantes `_transient`).
- Criterios de fluencia, endurecimiento y solvers son sólo 3D.
- Componentes `class(...), allocatable` en solvers (`Closest_point*`) y en leyes
  viscoplásticas (`hard_law`).

**Solvers**: `closest_point_iter` monta un sistema 7x7 y lo resuelve con MKL/LAPACK.
Versiones `closest_point`, `_2`, `_time`, `_perzyna`, `cutting_plane` + un `*.f90_bkup`.

**Build**: CMake + Ninja, `-cpp` ya activo (gfortran), sin FYPP. FYPP no instalado.

## 3. Arquitectura objetivo (resumen)

```
include/muscle_common.fypp         REAL_KINDS, DIMS, macros de nombres, @:DEVICE(), registro de modelos
tools/fypp/muscle_meta.py          Layouts de tensores (índices, pesos, etiquetas) + helpers de expresiones
cmake/MuscleFypp.cmake             muscle_add_fypp_sources(), add_muscle_fypp_test()
src/common/muscle_kinds.fypp       sp, dp, tolerancias y EPS_ABS por kind
src/tensors/types/*.fypp           Un template por familia (2Osym, 2O, 4O3sym, 4O2sym, 4O, iden_*) × DIMS × KINDS
src/tensors/operators/*.fypp       Operadores generados desde tablas de pares de tipos
src/tensors/muscle_tensors.f90     Fachada: re-exporta *_sp, *_dp y alias históricos (= dp)
src/<modelos>/*.fypp               Modelos concretos, genéricos en kind y dimensión
src/dispatch/*.fypp                Unión etiquetada + select case generado desde el registro
src/solvers/*.fypp                 Solvers monomorfizados / sobre la unión etiquetada
src/oo_facade/ (CPU)               API OO actual como envoltorio fino (transición)
```

Decisiones de diseño clave (detalle y alternativas en `DECISIONES.md`):
- D-001 FYPP como único generador; `.f90` generados no se versionan.
- D-002 Tipos distintos por kind con sufijo `_sp`/`_dp`, no PDT.
- D-003 Alias de compatibilidad: nombres históricos = doble precisión.
- D-004 Layouts actuales se conservan; 2D/3D se unifican con metadatos de layout.
- D-005 GPU: OpenMP offload portable (gfortran primero, luego ifx; NVIDIA, AMD e Intel). **Aceptada.**
- D-006 Despacho estático; límites a la monomorfización; `select case` fuera del kernel en GPU.
- D-007 Modelos sin OO dinámico: contrato de modelo; sin fachada OO. **Aceptada.**
- D-008 Materiales declarados (un kernel por material) + API por lotes SoA residente en GPU. **Aceptada.**

## 4. Fases

Cada fase se trabaja en ramas `fypp/f<N>-<tema>` desde `dev`, con commits pequeños y tests
en verde. **No se hace push** (lo hace el usuario). Marcar `[x]` al cerrar cada tarea y
anotar la fecha/commit.

### Fase 0 — Infraestructura de build
- [x] Añadir `fypp>=3.1` a `environment.yml` (pip). En esta máquina: fypp 3.1 vía apt. (2026-10-28)
- [x] `cmake/MuscleFypp.cmake`: `find_program` + versión ≥ 3.1, `muscle_fypp_preprocess()`,
      `muscle_add_fypp_sources()` (funciona aunque el target esté en otro directorio, varias
      llamadas por directorio), `DEPENDS` opcional, `MUSCLE_FYPP_DEFINES` (`ENABLE_UDTIO`,
      `USE_MKL`), `MUSCLE_REAL_KINDS` (por defecto `sp;dp`, exige `dp`), marcadores `#line`
      (D-009). `add_muscle_fypp_test()` en `tests/CMakeLists.txt`. (2026-10-28)
- [x] `include/muscle_common.fypp` (`ALL_KINDS`, `REAL_KINDS` validado, `REAL_TYPES`, `DIMS`, `tname`, `@:DEVICE()`) y
      `tools/fypp/muscle_meta.py` (`tname`). (2026-10-28)
- [x] `tools/check.sh [--lint-only] [BUILD_DIR]`: lint (sin comentarios ni cadenas; reglas GPU
      sólo en `src/`/`include/`, excepciones `! lint: allow(...)` / `! lint: host-only`) +
      build propio `build-check` (Debug, `sp;dp`) + ctest. (2026-10-28)
- [x] Prueba de humo sólo en tests: `tests/common/test_fypp_smoke.fypp` (`common-fypp_smoke`).
      Verificado: regenera al cambiar plantilla o `muscle_common.fypp`, el test detecta fallos,
      Release con `MUSCLE_REAL_KINDS=dp`, Ninja Multi-Config comparte el `.f90`, `@:DEVICE()`
      emite `!$omp declare target` con `-DMUSCLE_GPU`. El camino "plantilla dentro de
      `muscle`" se validó y lo ejercitará `muscle_kinds` en la Fase 1. (2026-10-28)
- [x] Revisión independiente (`migration-reviewer`): aprobada con cambios menores, aplicados
      (multi-config, defines, lint, validación de kinds, humo fuera de la librería). (2026-10-28)
- [x] Limpieza previa: borrar `closest_point_2.f90_bkup`; tabs → 4 espacios y espacios finales en
      `closest_point_2`, `math_linsolvers`, `vp_base` (`git diff -w` vacío). (2026-10-28)
- [x] Agentes `.claude/agents/`: `fypp-migrator`, `test-builder`, `migration-reviewer`. (2026-10-28)
- **Salida**: build idéntico + infraestructura FYPP funcionando; todos los tests en verde.

### Fase 1 — Kinds y utilidades comunes
- [x] `src/common/muscle_kinds.fypp` (`sp`, `dp`, `TOL_APPROX_sp = 1e-5`, `TOL_APPROX_dp = 1e-12`);
      `real(8)` → `real(dp)` en `muscle_math_linsolvers` y tests de OpenMP (en
      `test_muscle_openmp.f90`, `real64`: el test tiene una variable `dp`). (2026-10-28)
- [x] `tools/fypp/muscle_meta.py`: `layout(family, dim)` para los 10 tipos (los 4º orden se
      derivan de los de 2º orden con la regla de cada familia). Tests: Python
      (`tools-muscle_meta`, órdenes copiados de los `init2`) y Fortran generado
      (`tensors-types-layouts`, contra `init`/`get` de los tipos actuales, con prueba de
      mutación). Sin discrepancias. Los pesos de corte se añaden en la Fase 2. (2026-10-28)
- [x] `tests/common/muscle_test_utils.fypp`: `check_approx` real en todos los kinds; checks de
      tensores generados desde el layout (etiquetas `xxyz` en vez de `xx:yz` para 4º orden). (2026-10-28)
- [x] Revisión independiente: aprobada con cambios menores, aplicados — consultas de layout
      para la Fase 2 (`Layout.index`, `Comp.mult`, `Comp.voigt`, `Layout.base`, con tests),
      `EPS_ABS_<k>` en `muscle_kinds`, `DIMS`/`FAMILIES` desde `muscle_meta`, tests de contrato
      de inversas también en dp, fypp sin bytecode en el árbol fuente, D-010. (2026-10-28)
- [x] `src/math/muscle_math_inverses.fypp` + `muscle_math_lu_inverse.fypp` (macro): rutinas
      `_sp`/`_dp` con genéricos de nombre histórico; `QUIET_NAN` genérico = `dp`. (2026-10-28)
- [x] Tests `sp` de las inversas: `test_muscle_math_inverses.fypp` (instancia dp idéntica a la
      original; sp con tolerancias κ·ε, Hilbert sólo n=3,4, escalado 1e-12/1e+11, grupos nuevos
      `KAPPA` opcional y respaldo LU con matrices u·vᵀ+2⁻ᵖI). Umbral `ADJ_KAPPA_MAX_sp = 1e3`
      validado con 30k matrices aleatorias por kind: el adjunto en sp se comporta como en dp
      (error/(κε) ≤ 0.13 medio, ≈4–10 máx. hasta κ≈1e3). (2026-10-28)
- [x] Regla `legacy-kind` de `tools/check.sh` en todos los `.f90` (sin `real(8)`, `real*8`,
      `kind=8`). (2026-10-28)
- Hallazgos para la Fase 2–3: tolerancias por defecto de `is_approx` distintas por tipo
  (`1e-12`, pero `1e-8` en `4O2sym`) y escritas como literal real por defecto (`1E-12` es
  `real32` ≈ 9.9999998e-13): al templar, usar `1.0e-12_dp` cambia el umbral en ~2e-8 relativo
  (registrar en DECISIONES al migrar). `EPS_ABS=1e-30` también es literal por defecto.
- **Salida**: sin cambios de comportamiento; todos los tests en verde.

### Fase 2 — Piloto: tensores simétricos de 2º orden
- [x] Plantilla única `muscle_tensor_2osym.fypp` → `ten_3D2Osym_{sp,dp}`, `ten_2D2Osym_{sp,dp}`
      (patrón D-011). (2026-10-29)
- [x] Alias históricos (módulos de compatibilidad + `muscle_tensors.fypp`); tests existentes
      sin cambios. (2026-10-29)
- [x] Tests `sp` + embebido 2D↔3D (`tensors-types-2osym_kinds`); verificación bit a bit contra
      el original. (2026-10-29)
- [x] Benchmark `closest_point` Release: 0.0958 s → 0.0965 s (mín.), dentro del ruido. (2026-10-29)
- [x] Revisión independiente: "cambios necesarios", aplicados (procedimientos libres `type` +
      envoltorios `_tbp`, `@:DEVICE()`, nombres con kind al final, `ref`/`join_terms`/`specifics`,
      UDTIO aparte, reglas nuevas de `check.sh`). (2026-10-29)
- Notas UDTIO: con `-DMUSCLE_UDTIO=ON` y gfortran 15, `elasticity-linear_openmp` y
  `solvers-return_mapping-openmp_ts` fallan (el bug de UDTIO+OpenMP que motiva el AUTO); por eso
  `check.sh` sólo compila ese camino. gfortran 15 ignora los `/` del modo matriz en unidades internas.
- [x] Revisión del enfoque con el usuario antes de escalar: patrón D-011 aprobado. (2026-10-29)
- **Salida**: patrón validado y documentado en `DECISIONES.md`.

### Fase 3 — Resto de tipos
- [x] `2O` general (3D 9 comps / 2D 5 comps): `muscle_tensor_2o.fypp`, test `2o_kinds`. (2026-10-29)
- [x] `4O3sym` (21 / 10): `muscle_tensor_4o3sym.fypp`, test `4o3sym_kinds` (con UDTIO nuevo). (2026-10-29)
- [x] **PRIORIDAD antes de la segunda tanda (pedido del usuario, 2026-10-29): commit `refactor`,
      verificado bit a bit, registrado como D-012 (enmienda a D-011):**
      1. Macros comunes en `include/` para lo idéntico en las plantillas de tipo (~145 líneas:
         assign, sum, sub, subU, mul ×2, div, is_approx, approx, inv, envoltorios `_tbp` e
         interfaces de operadores) y aplicarlas a 2Osym, 2O y 4O3sym.
      2. `KAPPA` opcional en `M55INV`/`M66INV`/`M99INV` (con `COND_INF` tras el LU, mismos bits) para
         quitar los `#:if d` de las inversas de 4º orden; confirmar que no cambia el camino de `M44INV`.
      3. Helpers de tests a `tests/common/` (Park-Miller `rnd`, `inv_ref`/`eye`/`norm_inf`/`det_ref`
         n×n, `check_scaled`, `record`) y tests de `approx`/`free_generics` como macros.
      4. Consistencia: `Layout.matrix()` en las inversas en vez de Python inline (ampliar su test a
         4O/4O2sym), llamar siempre a la rutina específica `_<k>` de `muscle_math_inverses`, y
         comprobación de `pure` en el test de 4O3sym.
      Detalle completo: revisión de la tanda 1 (hallazgos 1–5).
      Hecho 2026-10-30 (D-012), con revisión aplicada: `init_tbp` con la forma del layout,
      `tol_default` por familia, test directo del `KAPPA` de `M55/M66/M99INV`.
- Instrucciones del usuario para la segunda tanda (2026-10-30, ejecución programada):
  sólo 4O2sym y 4O (identidades después); commits locales si `check.sh` pasa y el revisor
  aprueba (sin autoría de Claude, sin push); 4O2sym gana salida formateada igual que 4O3sym
  (lista + matriz 6x6/4x4, sólo con `ENABLE_UDTIO`); ante una decisión no cubierta o un cambio de
  comportamiento inevitable, detener esa parte, no commitearla y dejar la pregunta al usuario.
- [x] `4O2sym` (6x6 / 4x4) y `4O` general (9x9 / 5x5) — segunda tanda (ejecución programada,
      2026-10-30): plantillas `muscle_tensor_4o2sym.fypp`, `muscle_tensor_4o.fypp`; tests
      `4o2sym_kinds`, `4o_kinds`; macros de test con soporte `vals(n,n)`; 68/68; revisión aprobada.
- [x] **Pregunta abierta para el usuario (tanda 2) — resuelta 2026-10-30 (D-013: `ADJ_KAPPA_MAX_sp = 1e2`,
      tolerancia sp de 4O2sym se mantiene y se documenta):** el `test-builder` de 4O2sym midió que con
      `ADJ_KAPPA_MAX_sp = 1e3`, `M44INV_sp` por adjunta es ~8× menos preciso que LU para κ∞ en
      [3e2, 1e3] (p99 3.5e-5 vs 4.6e-6; máx 9.4e-5 vs 1.1e-5), al límite de la tolerancia sp 1e-4 de
      4O2sym; sugiere ~1e2 para sp (contradice la medición de la Fase 1, que lo daba por bueno).
      No se cambió nada. Además: la tolerancia sp 1e-4 de 4O2sym es adecuada hasta κ∞≈1e3 y
      estricta para κ≥1e4 (documentar: pasar `tol ≈ 2·κ·ε_sp` al comparar inversas mal
      condicionadas). Tolerancia por defecto de `4O2sym`
      decidida (D-010): 1e-8 dp / 1e-4 sp, con `tol_default`. Las macros de test
      (`muscle_test_tensor_cases.fypp`) indexan `vals` como vector: adaptarlas a `vals(n,n)` (p.ej.
      un `perturb_<T>`/`set_<T>` por índice lineal aportado por cada test).
- [x] Identidades `iden_2O`, `iden_2OS`, `iden_4O3T(S)`, `iden_4O4T(S)` (las escaladas llevan kind).
      Decisiones D-014 (2026-10-30). Plantilla `muscle_tensor_iden.fypp` + 6 módulos de
      compatibilidad; test `iden_kinds`; 69/69; bit a bit = original en dp; revisión aprobada.
      Álgebra de las de 4º orden (`+`, `-`, `.approx.`, `%init`): el usuario decidió no añadirla por
      ahora (2026-10-30); se agrega si un modelo la necesita (lista `ALGEBRA`, una línea + test).
- [x] UDTIO fuera de los módulos de tipo: `muscle_tensor_io` (todas las familias e identidades).

### Fase 4 — Operadores entre tipos
- [ ] Operadores mixtos tensor–identidad (hoy sólo dp con nombres históricos): al generarlos, las
      identidades sin escalar (`iden_2O`, `iden_4O3T`, `iden_4O4T`) se combinan con tensores de
      ambos kinds, y cada `iden_*S_<k>` sólo con tensores de su mismo kind (D-014).
- [x] `muscle_tensor_functions` (D-011): re-exporta las funciones libres de todas las familias
      (`norm`, `det`, `square`, `transpose`, `is_approx`, `inverse`, `convert_3sym`, `convert_2sym`)
      en un solo `use`, sin exportarlas desde `muscle_tensors`; test `tensors-functions`. Hecho al
      cerrar la Fase 3 (2026-10-30).
- La Fase 4 se trabaja en la rama nueva `fypp/f4-operators` desde `fypp/f3-types` (decisión del
      usuario, 2026-10-30).
- [ ] Tabla de pares por operador (`ddot`, `tdot`, `tdotsym`, `contraction_single`, `transform`,
      `assignment`, `addition_subtraction`) en `muscle_meta.py`.
- [ ] Generación de cada módulo de operadores desde su tabla × DIMS × KINDS.
- [ ] Verificar que el conjunto de interfaces `dp` generado ⊇ el actual (script que compara
      la lista de `module procedure`).

### Fase 5 — Math sin closures
- [ ] `muscle_math_derivatives`: reemplazar `procedure(...) :: func` por macros FYPP que
      generan la diferencia finita llamando directamente a la función concreta
      (`@:FD_GRADIENT(stress_eq_vonmises_dp, ...)`). Eliminar los procedimientos internos.
- [ ] `muscle_math_spectral_derivs`, `muscle_math_operations` genéricos en kind.
- [ ] Solver lineal pequeño propio (LU desenrollado ≤ 9x9) para uso en kernel; LAPACK queda
      sólo como opción CPU.

### Fase 6 — Modelos genéricos en kind y dimensión
Orden sugerido (de menos a más dependencias): hardening → viscoplastic → yield → elasticity
→ viscoelasticity → kinematics.
- [ ] Cada modelo concreto pasa a `type` concreto + procedimientos de módulo `pure`
      (`stress_eq_hill48_dp(params, stress)`), genérico en `k` y, cuando aplique, en `d`.
- [ ] `Base_elasticity`: `stress_3D`/`stress_2D` (+ `_transient`) desde una plantilla.
- [ ] Unificar firmas por familia (contrato D-007): mismos nombres y argumentos en todos los
      modelos de una familia. Condición previa de la Fase 7 (D-008).
- [ ] Retirar los tipos abstractos (`Base_*`): sin fachada OO (D-007/D-008).
- [ ] Maxwell: ramas en array fijo (`MAXWELL_MAX_BRANCHES`), sin allocatables.

### Fase 7 — Materiales declarados y API por lotes (D-008)
- [ ] Macro de contrato `#:call yield_model(...)` (y análogas por familia) — ver D-007.
- [ ] Cuerpos de solver (`closest_point`, `perzyna`, `time`) escritos una vez contra nombres
      genéricos (`yield_f`, `sigy`, ...) e instanciados por material vía `use` con renombre.
- [ ] Generador de materiales: registro `MATERIALS` (`-DMUSCLE_MATERIALS=...`), un módulo por
      material, `NSTATE` en compilación, registro de prueba que cubre cada modelo.
- [ ] API por lotes `muscle_update_<mat>` con arrays SoA y `map(present: ...)`; multifase
      por listas de índices.
- [ ] Estados (`Closest_point_2_data`, history) como arrays SoA de tamaño fijo por material.
- [ ] Adaptadores UMAT generados por material.
- [ ] (Opcional, posterior) uniones etiquetadas con `select case` generado si se necesita
      elegir el modelo en tiempo de ejecución.

### Fase 8 — Build GPU
- [ ] Opción CMake `MUSCLE_GPU={OFF,NVPTX,AMDGCN,SPIR64}` → flags de offload de gfortran
      (`-foffload=nvptx-none|amdgcn-amdhsa`) e ifx (`-fopenmp-targets=spir64`).
- [ ] Instalar los compiladores de offload de gfortran (paquetes `gcc-15-offload-nvptx` /
      `gcc-15-offload-amdgcn`): gfortran 15.2 está configurado con offload pero faltan los backends.
- [ ] `@:DEVICE()` expande a la directiva correspondiente.
- [x] Verificar `ieee_value` (usado por `QUIET_NAN_<k>`) en offload: funciona con ifx/spir64,
      pero con gfortran nvptx/amdgcn el enlace falla (`_gfortran_ieee_procedure_entry/exit` no
      existen en el dispositivo). Con `transfer(int(z'7FC00000', int32), 1.0_sp)` /
      `z'7FF8000000000000'` (bits idénticos a `ieee_value`) enlaza en ambos. (2026-10-30)
- [x] `QUIET_NAN_<k>` con `transfer` en vez de `ieee_value` en `src/math/muscle_math_inverses.fypp`
      (bits idénticos; sin el coste de ~150 ns del módulo IEEE). (2026-10-30)
- Prototipo GPU (2026-10-30, scratchpad, sin cambios en el repo): kernel con las 5 familias,
  identidades e inversas.
  - **Intel UHD 770 (ifx, spir64):** corre en GPU en sp; error máx. vs CPU 2.9 ε (3D) / 1.9 ε (2D)
    relativo a la escala del punto (FMA/orden). 200k puntos: kernel 0.21 s vs CPU 0.30 s (1 hilo) /
    0.039 s (8 hilos); 1ª llamada 4.2 s (JIT). dp sólo con emulación (`OverrideDefaultFP64Settings=1
    IGC_EnableDPEmulation=1`): kernels chicos funcionan, el completo no carga (límite de la
    emulación). Bugs de ifx/IGC encontrados: ver `.claude/rules/gpu-static-dispatch.md`.
  - **NVIDIA (nvptx-none) y AMD (amdgcn-amdhsa), gfortran 15.2:** todos los módulos compilan; el
    programa enlaza en sp y dp tras el cambio de `QUIET_NAN`. Sin hardware para ejecutar.
  - Pendiente: uso de registros por kernel (PTX/GCN) y ejecución en GPU con FP64 nativo.
- [ ] Test/benchmark: N puntos de integración en un `parallel loop` en GPU vs CPU.
- [ ] Prototipo FFT → MUSCLE → iFFT (cuFFT primero) con datos residentes: medir tráfico
      host↔device (sólo escalares) y rendimiento por vóxel; registros/spills por kernel.
- [ ] Auditoría automática (`tools/audit_kernels.sh`): grep de patrones prohibidos en la ruta de kernel.

### Fase 9 — Documentación y cierre
- [ ] FORD sobre las fuentes generadas (o `preprocessor: fypp` en `my_project.md`). **Urgente
      desde la Fase 1–2**: `my_project.md` usa `src_dir: ./src` y no ve los `.fypp`, así que
      los módulos migrados desaparecen de la doc (añadir `build/fypp/src` o un paso de generación).
- [ ] Actualizar `README.md`, `pages/tensors.md`, `pages/contributing.md` (cómo añadir un modelo al registro).
- [ ] Reescribir los usuarios externos de la API OO contra la nueva API.

## 5. Riesgos

| Riesgo | Mitigación |
|--------|-----------|
| Cambios de bits en `dp` por reordenamiento de operaciones | Plantillas emiten las mismas expresiones; tolerancia ≤1e-14 documentada |
| Explosión de tiempo de compilación (2 kinds × 2 dims × N pares) | `MUSCLE_REAL_KINDS=dp` para desarrollo; módulos de operadores divididos |
| Ambigüedad de genéricos al añadir `sp` (p.ej. `a * T` con literal sin kind) | Regla de literales con kind; compilar tests con `-Wconversion` |
| Soporte de OpenMP offload (gfortran/ifx) para operadores sobrecargados, funciones que devuelven tipos derivados y type-bound procedures | Prototipo temprano (al final de Fase 2) sobre `ten_3D2Osym_{sp,dp}` en nvptx y amdgcn |
| Precisión simple insuficiente en Newton de los solvers | Tolerancias por kind; evaluar si solvers `sp` tienen sentido o sólo tensores (decisión pendiente) |
| Desfase entre doc FORD y código generado | Fase 9; mantener `!!` dentro de plantillas |

## 6. Preguntas abiertas para el usuario

1. ~~GPU~~ → OpenMP offload, gfortran luego ifx, NVIDIA/AMD/Intel (D-005, 2026-10-25).
2. ~~Precisión~~ → **todo** (tensores, modelos, solvers) en `sp` y `dp` (2026-10-25).
3. ~~API OO~~ → sin fachada OO; materiales declarados + API por lotes en GPU (D-007, D-008, 2026-10-28).
4. ~~Git~~ → `.claude/` se mantiene en `.gitignore` (2026-10-25).

## 7. Registro de avance

| Fecha | Fase/tarea | Rama/commit | Notas |
|-------|-----------|-------------|-------|
| 2026-10-25 | Plan v1, CLAUDE.md, reglas | — (sin commit) | Línea base 53/53 |
| 2026-10-25 | D-005 aceptada, sp/dp en todo, D-007 propuesta | — (sin commit) | |
| 2026-10-28 | D-007 aceptada, D-008 (materiales declarados, API por lotes GPU) | — (sin commit) | CLAUDE.md en .gitignore |
| 2026-10-28 | Actualización desde `origin/dev` (`db0bc91`, +21 commits) | `fypp/f0-infra` | Nuevos: vp `chaboche`, viscoelástico `kelvin`; línea base 60/60 |
| 2026-10-28 | Fase 0 completa (infra FYPP, check.sh, agentes, limpieza) | `fypp/f0-infra` (sin commit) | 61/61; revisión aplicada; D-001 y D-009 aceptadas |
| 2026-10-28 | Fase 1 completa (kinds, layouts, test utils, inversas sp/dp) | `fypp/f1-kinds` (sin commit) | 63/63 (sp;dp y sólo dp); revisión aplicada; D-010 |
| 2026-10-29 | Fase 2 completa: piloto 2Osym (patrón D-011 aprobado) | `fypp/f2-ten2osym-pilot` | 64/64 (sp;dp y sólo dp), bit a bit = original salvo cambios D-011 |
| 2026-10-29 | Fase 3, tanda 1: 2O y 4O3sym (2 migradores + 2 test-builder en worktrees) | `fypp/f3-types` | 66/66; bit a bit = original; revisión: aprobar con cambios menores |
| 2026-10-30 | Refactor D-012 (macros comunes, KAPPA en inversas LU, helpers de tests) | `fypp/f3-types` | 66/66; bit a bit = eaee906; revisión aplicada |
| 2026-10-30 | Fase 3, tanda 2: 4O2sym y 4O (ejecución programada) | `fypp/f3-types` | 68/68 (sp;dp, sólo dp, UDTIO); bit a bit = original; revisión: aprobar con cambios menores |
| 2026-10-30 | D-013: `ADJ_KAPPA_MAX_sp = 1e2` (medición de confirmación), doc. tolerancia sp 4O2sym | `fypp/f3-types` | 68/68; dp sin cambios |
| 2026-10-30 | Fase 3: identidades (D-014) | `fypp/f3-types` | 69/69 (sp;dp, sólo dp, UDTIO); bit a bit = original salvo D-014.5; revisión aplicada |
| 2026-10-30 | `muscle_tensor_functions` (cierre de la Fase 3) | `fypp/f3-types` | 70/70 |
| 2026-10-30 | Prototipo GPU: Intel (corre, sp), NVIDIA/AMD (compila y enlaza) | scratchpad | `QUIET_NAN` con `transfer` aplicado; 2 bugs de ifx/IGC documentados |
