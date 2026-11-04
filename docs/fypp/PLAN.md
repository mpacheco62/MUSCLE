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

### Fase 4 — Operadores entre tipos (cerrada 2026-10-30)
- [x] Operadores mixtos tensor–identidad (hoy sólo dp con nombres históricos): al generarlos, las
      identidades sin escalar (`iden_2O`, `iden_4O3T`, `iden_4O4T`) se combinan con tensores de
      ambos kinds, y cada `iden_*S_<k>` sólo con tensores de su mismo kind (D-014).
- [x] `muscle_tensor_functions` (D-011): re-exporta las funciones libres de todas las familias
      (`norm`, `det`, `square`, `transpose`, `is_approx`, `inverse`, `convert_3sym`, `convert_2sym`)
      en un solo `use`, sin exportarlas desde `muscle_tensors`; test `tensors-functions`. Hecho al
      cerrar la Fase 3 (2026-10-30).
- La Fase 4 se trabaja en la rama nueva `fypp/f4-operators` desde `fypp/f3-types` (decisión del
      usuario, 2026-10-30).
- Plan de trabajo (2026-10-30; organización y kinds en D-015): 0) inventario de firmas dp actuales
      (script en `tools/`); 1) `test-builder`: golden dp de los 222 procedimientos contra el original
      (`+`/`-` y `=` hoy sin tests), commit aparte; 2) tabla de pares en `muscle_meta.py` + test;
      3) migración en tandas — A (paralelo): `addition_subtraction`+`assignment` |
      `contraction_single`+`dyadic`+`transform`; B: `contraction_double` — dp bit a bit, `check.sh`,
      revisión, commit; 4) tests `*_kinds` por módulo y embebido 2D/3D de cada operador; 5) borrar los
      12 módulos de compatibilidad (la fachada mantiene los alias; tests sólo cambian el `use`);
      6) humo GPU: prototipo con `.ddot.` en nvptx/amdgcn/ifx; 7) cierre de docs.
- Para el paso 3: los helpers privados 2D (`add_i3_2d4o2sym`, `add_i3_2d4o3sym`, `add_i4_2d4o2sym`,
      `add_i4_2d4o3sym`, `add_iden_2d2o` en `addition_subtraction`; `unpack_2d4o3sym` en
      `contraction_double`) no están en la tabla. En la tanda A los de `addition_subtraction` se
      eliminaron: su lógica sale de `iden_combination` (mismos bits). `transform_3D2O_3D4O3sym` ya no
      expande a 6x6 + kernel 4O2sym + `convert_3sym`: calcula el triángulo superior directo (mismos
      bits, ~20% más rápido) y el módulo deja de depender de `assignment`. Los cuerpos se
      eligen con `OpEntry.tag_a/tag_b/tag_res` + `op_iden`; los espejos no conmutativos (`-`,
      `2O*2Osym`, `4O2sym:2Osym`, `2Osym⊗I`) no reutilizan el cuerpo de la fila directa.
- [x] Tabla de pares por operador (`ddot`, `tdot`, `tdotsym`, `contraction_single`, `transform`,
      `assignment`, `addition_subtraction`) en `muscle_meta.py` (`OP_TABLE`, `op_entries`): expansión
      dp = las 222 firmas del inventario; sp 218 (D-015.4); test en `tests/tools/test_muscle_meta.py`;
      revisión aplicada. (2026-10-30)
- [x] Generación de cada módulo de operadores desde su tabla × DIMS × KINDS.
  - [x] Tanda A (2026-10-30): `addition_subtraction`, `assignment`, `contraction_single`, `dyadic`,
        `transform` → `.fypp` con `include/muscle_tensor_ops_common.fypp` (`OP_USES`,
        `OP_INTERFACES`, `OP_FUNCTION`, `OP_SUBROUTINE`) y helpers en `muscle_meta.py`; dp bit a
        bit = original (golden), salvo D-015.6; `check.sh` comprueba firmas y literales generados;
        revisión aplicada; 77/77.
  - [x] Tanda B (2026-10-30): `contraction_double` (`.ddot.`) → `.fypp` (1025 → 150 líneas);
        dp bit a bit = original a -O0/-O2/-O3 (golden + 42 procs × escalas 1e-320..1e308);
        unificado en 2D/3D con una sola regla (D-015.8, aprobado: sólo `4O3sym:4O3sym` cambia ≤ 1 ulp de la escala, ≤ 4e-15 relativo por
        componente);
        test `tensors-operators-ops_ddot_full` contra la contracción de tensores completos (84
        entradas, 16·ε(k)); 78/78. Revisión independiente pendiente (junto con el paso 5).
- [x] Verificar que el conjunto de interfaces `dp` generado ⊇ el actual: `ops_signatures.py --check`
      en `check.sh` (2026-10-30).
- [x] Paso 4a (2026-10-30): tests `tensors-operators-ops_<módulo>_full` (include común
      `tests/tensors/operators/ops_full_common.fypp`): cada entrada de la tabla (440 en total, 3D/2D,
      sp/dp) contra la operación sobre los tensores completos en dp; tolerancia 16·ε(k)·Σ|términos|
      (32 en `.transform.` de 4º orden), `=` exacto con el lado izquierdo pre-cargado; una mutación
      por operador detectada; 83/83.
- [x] Paso 4b (decisión del usuario 2026-10-30): cubierto por 4a. Cada entrada 2D se compara con la
      operación sobre los tensores completos 3x3/3x3x3x3 (operandos 2D embebidos, 13/23 nulos
      comprobados en el resultado), equivalente al test de embebido de `dims-2d-3d.md`.
- [x] Paso 5 (2026-10-30): borrados los 16 módulos de compatibilidad (`muscle_tensor_{3d,2d}*.f90`,
      `muscle_tensor_iden_*.f90`). La fachada `muscle_tensors` importa directamente de los módulos de
      familia los tipos con kind, los alias dp históricos (D-003) y los genéricos; los 12 tests que
      los usaban pasan a `use muscle_tensors` (sólo cambia el `use`); 83/83.
- [x] Revisión de la tanda B, tests 4a y paso 5 (2026-10-30), aplicada: un test aún usaba módulos
      borrados (`iden_kinds`, nombre templado) y pasaba por `.mod` obsoletos del build incremental →
      corregido, y `check.sh` borra ahora los `.mod` de módulos que ya no existen (build desde cero:
      83/83); cifra del cambio de bits de D-015.8 corregida. Pendiente menor: medir con `benchmarks/`
      `3D4O3sym:4O3sym` (el original compartía 30 temporales; el generado no factoriza `2*b`). Los
      fallos de tests OpenMP con `MUSCLE_UDTIO=ON` son el bug de gfortran ya conocido (`check.sh` sólo
      corre ahí los tests de tensores).
- [x] Paso 6 (2026-10-30): prototipo GPU (scratchpad `gpu-intel/`) con los operadores de la
      librería vía la fachada (`.ddot.`, `.tdotsym.`, `.tdot.`, `*`, `.transform.`, `+`/`-`, `=`):
      gfortran nvptx y amdgcn compilan y enlazan en sp y dp (los `muscle_tensor_ops_*` salen como
      funciones de dispositivo); Intel UHD 770 en sp: error ≤ 2.7 ε (3D) / 1.9 ε (2D) frente a CPU,
      100k puntos: kernel 0.14 s vs 0.25 s CPU 1 hilo / 0.055 s 32 hilos. Intel dp con emulación:
      el driver descarta los kernels grandes (igual con el prototipo anterior; no es regresión).
      Bug de ifx encontrado y resuelto: el `-` unario se perdía al importar `muscle_tensors`
      (fusión de `operator(-)` de familias + `ops_addition_subtraction`) → módulo
      `muscle_tensor_minus` que fusiona sólo `operator(-)`; verificado con ifx y gfortran.
- [x] Paso 7 (2026-10-30): doc FORD de los específicos de operadores en una línea (−5–10%: FORD
      genera siempre la tabla de argumentos); portada (`my_project.md`) y `pages/tensors.md`
      actualizadas a la API con kinds (ejemplos compilados).
- Pendientes que pasan a fases siguientes: medir `3D4O3sym:4O3sym` con `benchmarks/`; ejecución en
  GPU con FP64 nativo y registros por kernel (Fase 8); opción CMake `MUSCLE_GPU` (Fase 8).

### Fase 5 — Math sin closures
- Rama `fypp/f5-math` desde `fypp/f4-operators` (2026-10-31). Decisiones: D-016. Plan de trabajo:
      0) `test-builder`: golden dp de `muscle_math_{derivatives,operations,spectral_derivs}` y de
      `muscle_math_linsolvers` contra el original, commit aparte; 1) `spectral_derivs` +
      `operations` (eigenvals) + `dsquare` → `.fypp` sp/dp, `pure` + `@:DEVICE()` donde aplique;
      2) diferencias finitas: macros FYPP sin closures (sp/dp, 2D/3D) + API `procedure(...)` vieja
      en sp/dp host-only (D-016); 3) solver LU pequeño (A x = b, ≤ 9x9) sp/dp para kernels sobre la
      base de `muscle_math_lu_inverse`; `linsolvers` (Gauss, Gauss-Seidel, MKL) host-only;
      4) tests sp y de kernels; 5) revisión, docs y cierre.
- [x] Paso 0 (2026-10-31): golden dp de `derivatives`, `operations`, `spectral_derivs`, `linsolvers`
      (`tests/reference/math/`, driver `tests/common/muscle_test_golden.f90`); 87/87.
- [x] Paso 1 (2026-10-31): `spectral_derivs`, `operations` (`eigenvals`) y `dsquare` → `.fypp` sp/dp,
      `pure` + `@:DEVICE()` (compilan para nvptx); rutinas con `procedure(...)` en módulos `*_host`
      (D-016.8); umbrales sp por la regla tol_dp^0.442 (`REPEATED_TOL` 1e-4, `EPS_TOL` 1e-5,
      `EPS_TOL_SECOND` 2e-3), espectrales relativos a max λ² en sp y dp (D-016.6), `eigenvals` sp
      1e-9 (D-016.5), pasos FD por kind en `muscle_kinds` (D-016.7); revisión aplicada; golden sin
      cambios; 87/87. Pendiente (Fase 8): regla de `check.sh` que impida a un kernel llamar a un
      módulo `*_host`.
- [x] Paso 2 (2026-10-31): macros de diferencias finitas sin closures `FD_DERIV`, `FD_GRAD`,
      `FD_JAC`, `FD_HESS` (`include/muscle_fd.fypp`; sp/dp, 2D/3D, aptas para kernel, compilan para
      nvptx); API `procedure(...)` en sp/dp host-only construida con ellas
      (`muscle_math_derivatives_host.fypp`); pasos por kind en `muscle_kinds`; bug del argumento
      nulo corregido (D-016.4, regla del usuario: sólo cambia x = 0; 18 golden actualizados);
      `derivative_escalar_tensor3x3sym` con paso relativo a |x|; test `math-math_fd_macros`;
      helpers `norm_terms`/`family_module` en `muscle_meta`; revisión aplicada; 88/88.
- [x] Paso 3 (2026-10-31): solver LU pequeño `muscle_math_lu_solve` (`LU<n>SOLVE`, n = 2..9, sp/dp,
      `pure` + `@:DEVICE()`, `KAPPA` opcional; frente `lu_solve` para host/tests) sobre la LU compartida
      con `M<n>INV` (mismos bits verificados); `muscle_math_linsolvers` host-only sp/dp (`sgesv`/
      `dgesv`); `USE_MKL` sólo como define de fypp; test `math-math_lu_solve` (sp/dp, NaN, casi
      singular, solvers host en sp); D-016.10; revisión aplicada; 89/89.
- [x] Paso 4 (2026-10-31): tests `math-math_spectral_kinds` (6 rutinas espectrales sp/dp, escalas
      1e-3…2.1e5, autovalores distintos/casi repetidos/dobles/triples rotados, referencia dp
      independiente; falla con los umbrales absolutos antiguos de D-016.6), `math-math_fd_operations_kinds`
      y `math-math_kernels_target` (kernel con eigenvals, descomposición, dEigenvalues, `FD_GRAD` y
      `LU6SOLVE` dentro de `!$omp target`, bit a bit con CPU; host-fallback con `-foffload=disable`).
      Bug sp encontrado y corregido: `spectral_decomposition` perdía la normalización con cortes
      ~1e-11 (|c|² subnormal) → escalado exacto por potencia de 2 (dp mismos bits, golden sin
      cambios); 92/92. Para la Fase 8: `eigenvals` divide por cero (1/0 = Inf, resultado correcto)
      con la raíz mayor exactamente repetida; con trampa de división por cero aborta.
- [x] `muscle_math_derivatives`: macros FYPP que generan la diferencia finita llamando
      directamente a la función concreta (`@:FD_GRAD(dfds, stress_eq_vonmises_dp(y, {x}), sigma, dp,
      3D)`, `include/muscle_fd.fypp`) — paso 2. Eliminar los procedimientos internos de los modelos
      y la API `procedure(...)` pasa a la Fase 6 (D-016.1).
- [x] `muscle_math_spectral_derivs`, `muscle_math_operations` genéricos en kind — paso 1.
- [x] Solver lineal pequeño propio (LU desenrollado ≤ 9x9) para uso en kernel; LAPACK queda
      sólo como opción CPU — paso 3.
- [x] Paso 5 (2026-10-31): revisión final de la fase aplicada (docs coherentes con el código,
      `kspecifics` en todas las plantillas math, doc FORD de argumentos de linsolvers, D-016.11);
      build sólo dp 92/92; paso en x = 0 de las primeras derivadas revisado a √ε (D-016.4).
      **Fase 5 cerrada (2026-10-31).** Pendientes que pasan a otras fases:
      - Fase 7 (antes Fase 6): tangente numérico de los return mappings (`closest_point`, `_time`, `_perzyna`): pasar
        un paso explícito ligado al material o usar el tangente consistente analítico (D-016.4, sp).
      - Fase 7 (antes Fase 6): borrar la API `procedure(...)` (`*_host`) al pasar los modelos a las macros (con
        tolerancia, D-016.9); golden mínimo del primer incremento (deformación, `ep`, `epd` nulos) de
        los return mappings y de `muscle_vp_base`/`muscle_vp_overstress_base` antes de migrarlos.
      - Fase 8: regla de `check.sh` que impida a un cuerpo con `@:DEVICE()` llamar a nombres
        exportados por un módulo `*_host` (ojo: `muscle_math_operations` y
        `muscle_math_derivatives` re-exportan sus `*_host`, así que una regla por `use` daría falsos
        positivos); división por cero de `eigenvals` con raíz mayor repetida; ejecución real en GPU.

### Fase 6 — Modelos genéricos en kind y dimensión
- Rama `fypp/f6-models` desde `fypp/f5-math` (2026-10-31). Decisiones: D-017. Plan de trabajo por
      familia, en el orden hardening → viscoplastic → yield → elasticity → viscoelasticity →
      kinematics: 0) `test-builder`: golden dp de la familia contra el OO original (commit aparte);
      1) macro de contrato de la familia (el piloto, hardening bilinear, la diseña) + modelos de la
      familia en sp/dp con funciones libres `pure` + `@:DEVICE()`; 2) adaptadores OO host-only que
      llaman a las funciones nuevas (solvers/tests/UMAT sin cambios); 3) tests sp; 4) revisión y
      commit. Golden del primer incremento (deformación, `ep`, `epd` nulos) de return mappings y
      leyes viscoplásticas antes de migrarlas (pendiente de la Fase 5).
- [x] Familia hardening (2026-10-31): golden dp (`plasticity-hardening_laws-hard_golden`); patrón
      invertido de contrato (D-017.6: macros chicas en `include/muscle_models_common.fypp`, metadatos
      de familia en `muscle_meta`); bilinear, Ludwik, Swift y Voce modificada en sp/dp con funciones
      libres `pure` + `@:DEVICE()` y genéricos `hard_stress`/`hard_dstress_dep`/`hard_ddstress_ddep`
      (`muscle_hard_laws`); adaptadores OO generados (D-017.5c); Ludwik offset 1e-10 en ambos kinds
      (D-017.4); FD con `scale`/`lower`; dp bit a bit; revisiones aplicadas; 94/94.
- [x] Familia viscoplástica (2026-10-31): golden dp (`plasticity-viscoplastic_laws-vp_golden`, 356
      registros, incluido ep/epd = 0); familias `vp` (JC, JCfull, MRK, NNL, RK, VA) y `os` (Perzyna,
      Chaboche) en sp/dp con composición desacoplada (D-017.6b) y API de parciales + total con dt
      (D-017.8); bugs del original documentados, no corregidos (D-017.7); adaptadores OO; dp bit a
      bit; revisión aplicada; 96/96.
- [x] Familia yield, piloto (2026-11-01): golden dp (2012 registros, 37 estados); familia `yield`
      (D-017.10) con von Mises, Hill48 y CPB06 en sp/dp, evaluación conjunta con `ddf` opcional,
      parámetros derivados privados, `MODEL_FD` sobre 2Osym; correcciones aprobadas D-017.9 (f = 0)
      y D-017.10c/d; revisión aplicada; 99/99.
- [x] Familia yield completa (2026-11-01): Drucker-Prager, CPB06ex2 y Yld2004 (evaluación conjunta
      propia con `ddf` opcional; código común en `muscle_anisotropic_yield`); parámetros derivados
      privados con `yield_<m>_init` en los cuatro anisótropos/DP; decisiones D-017.11 (redondeo del
      reescalado aceptado, umbral sp DP 1e-5, hessiana estable de DP); revisiones aplicadas; 102/102.
- [x] Familia elasticity, piloto (2026-11-01): golden dp (1311 registros, 3D y 2D); familia `elas`
      con `Elasticity_linear` en una plantilla 2D/3D sp/dp (D-017.12d); un solo camino analítico para
      los desviadores (D-017.12a/c, cambios aprobados) y ν = 0.5 documentado (D-017.12b); revisión
      aplicada. Caché del tangente sólo en el adaptador OO (D-017.12e).
- [x] Familia viscoelasticity, piloto (2026-11-01): golden dp (historias de varios pasos, 2188
      registros); familia `visco` con estado explícito (D-017.13d), Maxwell y Kelvin sp/dp,
      `elas_frozen`, status común; D-017.13a/b aplicados (golden verificado con un oráculo
      independiente); revisión aplicada; 107/107.
- [x] Familia kinematics (2026-11-01): golden dp (2071 registros, 20 F + secuencias); familia `kin`
      con estado `kin_state` y status (D-017.14), las 6 cinemáticas en sp/dp (small strain, TL, UL,
      logarítmicas material y espacial, corrotacional con `kin_corot_update`); adaptadores OO con
      los mismos bits; revisión aplicada. Ya no queda ninguna cinemática en `.f90` (sólo bases y
      adaptadores host y la fachada).
- Para la Fase 7 (de la revisión de kinematics): `kin_update` de TL/UL calcula la polar completa en
  cada punto aunque no use R (necesario hoy para el status y los mismos bits): valorar una variante
  sin R tras perfilar; no hay equivalente de `update_incremental` en el kernel (el solver construye
  F_{n+1}); falta el término geométrico/Jaumann de los tangentes (USERMAT de ANSYS, D-017.14b).
- Para la Fase 7 (de la revisión de visco): el estado de tamaño fijo copia 6 + 6·MAX_VISCO_BRANCHES
  reales por punto (`state_np1 = state_n`); con MAX = 10 y dos estados por hilo son ~130 reales
  (riesgo de derrame de registros): considerar un máximo de ramas por material declarado. Unificar
  `STATUS_OK = 0` con los códigos de resultado de los solvers históricos (`STATUS_ELASTIC_CASE = 0`).
Orden sugerido (de menos a más dependencias): hardening → viscoplastic → yield → elasticity
→ viscoelasticity → kinematics.
- [x] Cada modelo concreto pasa a `type` concreto + procedimientos de módulo `pure`, genérico en `k`
      y, cuando aplique, en `d` (patrón invertido, D-017.6).
- [x] `Elasticity_linear`: `stress_3D`/`stress_2D` desde una plantilla (D-017.12d).
- [x] Unificar firmas por familia (contrato D-007): genéricos `<fam>_<func>` iguales en todos los
      modelos de una familia (D-017.5b).
- [x] Estado plástico `plastic_state_<k>` (`muscle_plastic_state`, módulo de soporte sin familia):
      tamaño fijo, pack/unpack `pure` + `@:DEVICE()`, sin status (no valida);
      `Plastic_material_history` queda como adaptador host (2026-11-01).
- [x] Revisión final de la fase (2026-11-01): build sólo dp 111/111; sin modelos en `.f90` salvo
      bases/adaptadores host, fachadas y solvers. **Fase 6 cerrada (2026-11-01).** Pasan a la Fase 7:
      retirar los `Base_*` y los adaptadores OO (D-017.1); tangente numérico y golden del primer
      incremento de los return mappings; borrar la API `procedure(...)` (aún la usan
      `muscle_vp_overstress_base` y los `closest_point*`); composición múltiple en `ModelFamily`
      (yield con endurecimiento isótropo + cinemático) y valores compuestos tensoriales; decidir si se
      unifican el status de los `*_init` (yield/elas no lo devuelven; visco sí) y el pack/unpack de
      kernel (plastic_state sí; visco no). Fase 8: comprobar con el compilador de offload que los
      módulos de kernel compilan aunque hagan `use` de las bases OO (D-017.5c) antes de retirarlas.
- [x] Maxwell: ramas en array fijo (ya era `MAX_VISCO_BRANCHES = 10`, sin allocatables; D-017.13c).

### Fase 7 — Materiales declarados y API por lotes (D-008)
Diseño y alcance: D-018. Pasos (cada uno un commit con el ciclo completo):
1. Golden dp de los cinco solvers OO (`tests/reference/solvers/`): tensión, ep, `strain_p`, status,
   `iters` y tangentes analítico y numérico; carga, descarga, incremento grande y primer incremento.
2. Registro `tools/fypp/muscle_materials.py` + validador + tests Python; CMake genera N módulos.
3. Códigos de status de solver en `muscle_status` y tolerancias del Newton por kind (relativas, D-018.3).
4. Cuerpo `closest_point` + plantilla de material + kernel de punto sp/dp (tangente analítico).
5. `perzyna`, `time`, `closest_point_2` (LU propio) y `cutting_plane`.
6. Materiales con elasticidad visco congelada y con cinemáticas finitas.
7. API por lotes SoA + `idx` + OpenMP CPU (lote == bucle de puntos).
8. API genérica de integración con programas FE (D-018.17): init desde array, API de punto plana,
   rotación del hsv y ejemplos en la documentación.
9. Migración de los tests de solvers, retiro del OO (D-018.8) y cierre.
   Etapas del paso 9 (D-018.18), un commit cada una:
   - E1: separar `solvers_golden_common` (helpers neutros / fábricas OO), ~15 materiales nuevos,
     pruebas físicas de solvers en kernel y benchmark con la API nueva.
   - E2: goldens de solver generados desde el kernel; borrar los tests OO y los 5 solvers OO.
   - E3: una familia por commit: vp+os → hard → yield → visco (+`Viscoelastic_history`) → elas →
     kin (+fachada) → `Plastic_material_history`.
   - E4: borrar math `*_host`; goldens de derivadas con las macros.
   - E5: quitar las macros OO y los helpers `oo_*` de `muscle_meta`; endurecer `check.sh`.
   - E6: docs (README, `pages/`, `my_project.md`, `benchmarks/README`, `.claude/rules`, tabla de
     equivalencias) y cierre de la fase.
- [x] Paso 1 (2026-11-01): golden dp de los cinco solvers (1533 registros, mutaciones de 1 ulp).
- [x] Pasos 2–4, piloto `closest_point` (2026-11-02):
      - registro `muscle_materials.py` (6 materiales: todos los criterios y endurecimientos, elasticidad
        lineal) y generación CMake de `muscle_mat_<name>`;
      - `STATUS_NOT_CONVERGED` y `NEWTON_TOL_<k>`/`TOL_ELASTIC_<k>`;
      - cuerpo `muscle_closest_point_body.fypp` y kernel de punto sp/dp con tangente analítico;
      - hsv de 7 reales; dp bit a bit con el golden actualizado (217 registros, D-018.3/10a/11).
- [x] Paso 5a, `perzyna` y `time` (2026-11-02):
      - cuerpos `muscle_closest_point_{perzyna,time}_body.fypp`, con hsv común y firma uniforme con `dt`;
      - registro con 7 materiales perzyna y 12 time, que cubren todas las leyes os/vp;
      - goldens actualizados según D-018.3/10b/12;
      - `vp_laws_dp.txt`: 18 registros jcfull con `dstress_dep` analítico.
- [x] Paso 5b, `closest_point_2` y `cutting_plane` (2026-11-02):
      - cuerpos propios;
      - `cp2` con `LU7SOLVE`, sin LAPACK (D-018.13);
      - 8 materiales cp2 y 7 cut en el registro.
      Los casos con elasticidad Maxwell (caso 7 de cada uno) pasan al paso 6.
- [x] Paso 6a, elasticidad visco congelada (2026-11-02):
      - los cinco solvers con Maxwell/Kelvin, todos los casos del golden bit a bit;
      - `mat_step` una vez por paso;
      - hsv 13 + 6·`visco_branches`;
      - tests de los caminos de error (D-018.14).
- [x] Paso 6b, cinemáticas finitas (D-018.15) (2026-11-02):
      - `kin_geometric_tangent` y tangentes espaciales corregidos (105 registros kin);
      - sp robusto (Jacobi y polar de Higham);
      - rotación ΔR de `strain_p` y marco rotado para la anisotropía;
      - corrotacional rotada por el FE;
      - 13 materiales kin con closest_point, perzyna, time y cp2;
      - oráculos físicos y golden de regresión (260 registros).
- [x] Paso 7, API por lotes (2026-11-02):
      - `muscle_update_<name>` SoA y macro `BATCH_LOOP` (CPU OpenMP / `target` con `map(present, alloc)`);
      - tests de lote == bucle de kernels de punto (12 materiales, 1 y 4 hilos) y OpenMP;
      - la ruta GPU compila pero no se ejecutó (Fase 8) (D-018.16).
- [x] Paso 8, API genérica de integración con programas FE (D-018.17) (2026-11-02):
      - `muscle_point_<name>`, `mat_rotate_hsv_<name>` y `mat_init_<name>` + `MAT_<NAME>_NPROPS`;
      - parámetros con fuente única en `muscle_meta`, más tabla `muscle_meta.py props`;
      - esqueletos Abaqus/ANSYS/LS-DYNA en `pages/integration.md`. Las notas de VUMAT y LS-DYNA están
        sin verificar contra esos programas.
- [x] Paso 9, retiro del OO (D-018.18), etapas E1–E6 (2026-11-03):
      - E1 697db87; E2 f854819;
      - E3: vp 2803e87, hard cbe2850, yield 4ff2e2f, visco 93863bc, elas 206cec0, kin 90fab26,
        plastic history f080c4d;
      - E4 ce6af37; E5 755d984; E6 9e946b1.
- [x] Revisión final de la fase (2026-11-03): 115/115 en sp;dp y en sólo dp; el registro de usuario
      con `-DMUSCLE_MATERIALS` compila; correcciones menores aplicadas. **Fase 7 cerrada.**
  - Las tareas de la lista original quedaron así:
    - el cuerpo único con renombre en `use` quedó enmendado por D-018.5 (genéricos + tipos);
    - los adaptadores UMAT generados se reemplazaron por la API genérica (D-018.17);
    - los estados SoA de tamaño fijo se resolvieron con hsv por material (D-018.7/14).
  - Pasan a la Fase 8:
    - camino `is_device_ptr` (memoria de dispositivo pura, cuFFT);
    - verificar el modificador `present` en nvhpc; OpenMP 6 depreca `alloc` en favor de `storage`
      (D-018.16f);
    - achicar `visco_state h(6, MAX_VISCO_BRANCHES)` por punto (D-018.14f);
    - ejecutar `BATCH_LOOP` con `target` en hardware real y medir registros/spills de los kernels de
      material;
    - variante TL/UL sin la descomposición polar.
  - Pasan a la Fase 9 o después:
    - materiales 2D (D-018.11g);
    - validar los esqueletos Abaqus/ANSYS/LS-DYNA contra los programas reales (D-018.17d/g, los
      `ustatev` de ANSYS);
    - composición múltiple en `ModelFamily` (endurecimiento cinemático);
    - unificar el status de los `*_init`;
    - unión etiquetada opcional para elegir el modelo en tiempo de ejecución (D-008.4);
    - convergencia de Chaboche n=0.5 + DP en compresión (PROBLEMAS).

### Fase 8 — Build GPU
Alcance y verificación: D-019. Pasos (cada uno un commit con el ciclo completo):
1. Opción CMake `MUSCLE_GPU` + flags (gfortran nvptx/amdgcn, ifx spir64); `@:DEVICE()` y `BATCH_LOOP`
   con el define; los tests enlazan con offload; `check.sh` construye los 3 targets.
2. Tests target vs CPU por kernel en ambos kinds (tensores, math, familias, materiales y lote con
   `BATCH_LOOP`), ejecutados en la Intel (sp; dp donde aguante la emulación).
3. Bugs de kernel: `eigenvals` divide por cero con una raíz repetida, y los que salgan del paso 2.
4. Auditoría de la ruta de kernel (`tools/audit_kernels.sh`: símbolos de runtime host en las imágenes
   de offload y patrones prohibidos).
5. Traspasos de la Fase 7:
   - camino `is_device_ptr`;
   - `present`/`storage` (OpenMP 6);
   - achicar `visco_state` por material;
   - TL/UL sin polar;
   - registros/spills por kernel.
6. Benchmark de N puntos GPU vs CPU (Intel sp) y cierre.
- [x] Opción CMake `MUSCLE_GPU={OFF,NVPTX,AMDGCN,SPIR64}` → flags de offload de gfortran
      (`-foffload=nvptx-none|amdgcn-amdhsa`) e ifx (`-fopenmp-targets=spir64`).
- [x] Instalar los compiladores de offload de gfortran (ya instalados al empezar la Fase 8) (paquetes `gcc-15-offload-nvptx` /
      `gcc-15-offload-amdgcn`): gfortran 15.2 está configurado con offload pero faltan los backends.
- [x] `@:DEVICE()` expande a la directiva correspondiente.
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
- [x] Test/benchmark: N puntos de integración en un `parallel loop` en GPU vs CPU.
- [→] (Pasa a la Fase 9, D-019.2) Prototipo FFT → MUSCLE → iFFT (cuFFT primero) con datos residentes: medir tráfico
      host↔device (sólo escalares) y rendimiento por vóxel; registros/spills por kernel.
- [x] Auditoría automática (`tools/audit_kernels.sh`): grep de patrones prohibidos en la ruta de kernel.
- [x] Pasos 1–6 (2026-11-03):
  - 1: build GPU (45c1336);
  - 2: tests target (460fc99, 7c082dc);
  - 3–4: `eigenvals` y auditoría (70d1051);
  - 5: `_dev`, estado visco por material, TL sin polar y recursos por kernel (b3eb490, e565aa8,
    d6bff03, 7206dd7);
  - 6: benchmark (c9515ef).
- [x] Revisión final de la fase (2026-11-03); correcciones menores aplicadas. **Fase 8 cerrada.**
  Pasan a la Fase 9:
  - el prototipo FFT;
  - la ejecución en NVIDIA/AMD por parte del usuario;
  - volver a habilitar la lista de exclusiones de GPU cuando se actualice IGC;
  - el error interno de ifx en `p%init` (D-019.4f);
  - `storage` de OpenMP 6;
  - UL/SL sin polar (haría falta ΔR sin R);
  - el benchmark de CPU con ifx (`omp_target_alloc`);
  - el fallo en Release de `pz_*`/`tm_jc` en la GPU Intel (D-019.8);
  - reducir el spill.

### Fase 9 — Legibilidad y mantenibilidad (D-020)
Objetivo: que el código se lea como Fortran. **No agrega funcionalidad.** FYPP queda en lo que Fortran
no puede hacer; tipos, firmas, declaraciones, cuerpos y cierres se escriben en Fortran. Criterio:
`.claude/rules/fypp-templates.md` («FYPP mínimo», con las relajaciones de la enmienda de D-020).
Rama `fypp/f9-readability` desde `fypp/f8-gpu`. Alcance: todo `src/`, `include/`, `tools/fypp/`,
tests y la infraestructura de build.

**Tamaño de partida (2026-11-04):**

| Área | Plantillas | Líneas | Con macros |
|---|---|---|---|
| tensors/types, operators | 12 | 2283 | 12 |
| tensors (fachada, io, funciones) | 3 | 446 | 0 |
| math | 7 | 1532 | 5 |
| hard, vp/os, yield (+ plastic_state) | 24 | 2843 | 24 |
| elasticity, visco | 7 | 775 | 7 |
| kinematics | 8 | 1467 | 8 |
| solvers (cuerpos) | 6 | 1054 | 6 |
| materials | 4 | 1174 | 3 |
| tests `.fypp` | 73 | 22075 | — |

Macros más usadas: `@:DEVICE` 152 (se queda), `MODEL_FUNCTION` 115, `MODEL_PUBLIC` y
`MODEL_PARAMS_TYPE` 28 cada una, `YIELD_ZERO_AT_NULL` 10, `FAMILY_USES`/`FAMILY_EXPORTS`/`MODEL_INIT`
7 cada una. `include/*.fypp` suma ~1190 líneas y `muscle_meta.py` 2467 (111 funciones, 5 tests Python).

**Método de cada paso** (el ciclo del proyecto, con una red de seguridad extra):
1. Línea base: guardar el `.f90` generado de las plantillas que se tocan.
2. Reescribir la plantilla a mano.
3. **Comparar el generado antes/después**, normalizado (sin comentarios, `#line`, espacios, ni el orden de
   las declaraciones): se espera el mismo código; lo que cambie se revisa a mano. Además, la tabla de
   símbolos públicos (`nm`) de la librería debe ser idéntica.
4. `tools/check.sh` completo: goldens dp bit a bit, tests target, los tres builds GPU, auditoría.
5. `migration-reviewer` mira legibilidad (¿se entiende leyendo sólo ese archivo?) y deriva de diseño.
6. Aprobación del usuario → commit local. Un commit por unidad (familia, solver, etc.).
Nada de cambios de bits: si el generado difiere y hay cambio numérico, se detiene y se pregunta.

**Etapas:**
- [ ] **E0. Preparación.** Rama; script `tools/gen_snapshot.sh` (guarda y compara el `.f90`
      normalizado y la tabla de símbolos); `tools/fypp_metrics.sh` (nº de macros, líneas de
      `include/` y `muscle_meta.py`, plantillas con macros) para ver el avance de cada etapa.
- [ ] **E1. Piloto: `muscle_hard_swift` y `muscle_hard_bilinear`** (216 líneas en total las 5 leyes).
      Se escriben tipo, funciones y cierres en Fortran con `${k}$`. Preguntas de diseño a resolver
      aquí y a registrar: cómo se declaran a mano los genéricos de familia (`hard_stress`...) y
      `public`; dónde va `@:DEVICE()`; formato de la doc FORD; cómo se escribe el `*_init`.
      **Pausa: se revisa con el usuario cómo quedó y se ajusta la regla antes de seguir.**
- [ ] **E2. Inventario de macros y metadatos.** Tabla con cada macro de `include/*.fypp` y cada
      función de `muscle_meta.py`: usos, veredicto (se queda, se reescribe en Fortran, pasa a ser un
      dato, o la comprueba el compilador o un test). Se registra en DECISIONES.
- [ ] **E3. Familias, una por commit** (de menor a mayor riesgo):
      hard (5, 216 líneas) → elas (2) → vp/os (10, 927) → visco (5, 621) → yield (8, 1607; el más
      grande y con la normalización por exponente y el cero en nulo) → kin (8, 1467).
      Al terminar una familia se borran sus macros de contrato ya sin uso.
- [ ] **E4. Cuerpos de solver** (cp, cp2, perzyna, time, cutting plane; 6 plantillas, 1054 líneas):
      subrutinas Fortran legibles con `${k}$`, un nivel de `#:include`. Se elimina `SOLVE_DUMMIES` y
      lo que arme declaraciones; el tangente compartido queda como texto fijo.
- [ ] **E5. Plantilla de material y APIs** (`muscle_material*.fypp`, 1174 líneas): separar lo que es
      generación por material (justificada: tipos, `mat_init`, lote, punto) del esqueleto, que se
      escribe en Fortran. `mat_init` asigna por nombre y el compilador verifica los nombres.
- [ ] **E6. Tensores, operadores y math** (21 plantillas, 4261 líneas): con el criterio laxo,
      `TENSOR_COMMON_*` se conserva si pasa la prueba de «texto fijo con sustitución de nombres»;
      se revisan `ops_common`, las listas de específicos y las plantillas de math caso a caso.
- [ ] **E7. `muscle_meta.py` e `include/`.** Borrar lo que quedó sin uso, pasar a datos lo que sea
      lista de nombres, y reducir los tests Python a lo que queda. Se mide con `fypp_metrics.sh`.
- [ ] **E8. Tests.** Reemplazar los tests generados por macros por tests escritos a mano cuando sean
      más legibles (73 plantillas, ~22 000 líneas); los drivers de goldens y `tests/gpu/target_cases.py`
      se quedan donde evitan duplicación real.
- [ ] **E9. Infraestructura.** Simplificar `tools/check.sh` (lint, builds GPU) y `cmake/` donde la
      complejidad haya sobrado, sin perder cobertura; documentar en `check.sh --help`.
- [ ] **E10. Documentación y cierre.** `pages/contributing.md` (cómo se escribe un modelo con el
      criterio nuevo, con un ejemplo completo), `pages/tensors.md`, README, reglas de `.claude/` y
      CLAUDE.md; revisión final de la fase; cifras antes/después.

**Riesgos:**
- *Deriva de bits.* Mitigación: el generado normalizado y los goldens; ningún paso se aprueba con
  diferencias sin explicar.
- *Regresión en GPU.* Mitigación: `check.sh` construye siempre nvptx, amdgcn y spir64 y corre los
  tests target en la Intel.
- *Duplicación a mano entre sp y dp.* Se evita con `#:for k`: se escribe una vez, con `${k}$`.
- *Cansancio de reescritura.* Se trabaja por familias pequeñas, con pausa en el piloto.
- *Pérdida de unificación.* Algunas tablas únicas dejan de generarse (p. ej. orden de props); queda
  documentado y verificado por el compilador y por un test (D-020).
**Criterios de salida:** mismos bits en todos los goldens, mismos símbolos públicos, tests target y los
tres builds GPU en verde; ninguna macro prohibida por la regla; un modelo nuevo se escribe y se lee en
un solo archivo; y cifras de «macros, líneas de `include/` y de `muscle_meta.py`» antes y después.

## Pendientes fuera del plan
Ya no pertenecen a una fase. Se retoman a pedido del usuario, y cada uno con su decisión:
- Prototipo FFT → MUSCLE → iFFT (cuFFT primero) con datos residentes, midiendo el tráfico
  host↔device (D-019.2).
- Ejecución en NVIDIA/AMD por parte del usuario; volver a habilitar los tests excluidos de
  `tools/gpu_skip_tests.sh` cuando Intel actualice IGC; error interno de ifx en `p%init` (D-019.4f);
  `storage` de OpenMP 6; benchmark de CPU con ifx (`omp_target_alloc`); fallo en Release de
  `pz_*`/`tm_jc` en la GPU Intel (D-019.8).
- UL y logarítmica espacial sin la descomposición polar (haría falta ΔR sin R).
- Reducir registros/spill de los kernels de material, midiendo en hardware real NVIDIA/AMD con
  `tools/kernel_resources.sh` (D-019.7h). Opciones, de más barata a más cara:
  1. flags: Intel `-ze-opt-large-register-file`, NVIDIA `maxrregcount`/`launch_bounds`, tamaño de
     bloque en AMD;
  2. menos estado vivo, por ejemplo sin tangente cuando no se pide o con el Jacobiano 7×7 de cp2;
  3. recalcular C en lugar de guardarlo;
  4. separar el return mapping del cálculo del tangente;
  5. controlar el inlining de las funciones de punto en amdgcn.
- Materiales 2D (D-018.11g); validar los esqueletos Abaqus/ANSYS/LS-DYNA contra los programas reales
  (D-018.17); composición múltiple en `ModelFamily` (endurecimiento cinemático); unificar el status
  de los `*_init`; unión etiquetada opcional (D-008.4).
- Reescribir los usuarios externos de la API OO (guía: `pages/migration.md`).
- Convergencia de Chaboche n=0.5 + DP en compresión y de `tm_nnl_vm` ×40 (PROBLEMAS).
- Opcionales de la Fase 6: forma estable de `stress_eq` y del gradiente de Drucker-Prager en sp,
  borrar `finite_get_E_log` y `finite_get_b`, medir con benchmarks el `.ddot.` 3D4O3sym:4O3sym.
Hecho en su día (se conserva como historia):
- [x] FORD sobre las plantillas (2026-10-30, adelantado a la Fase 4 a pedido del usuario):
      `my_project.md` añade `fypp` a `fpp_extensions` y `preprocessor: fypp -I include -M tools/fypp
      -m muscle_meta`. FORD 7.0.13 preprocesa los `.fypp` con todas las dimensiones y kinds: la doc
      pasa de 66 a 81 módulos (los 15 migrados faltaban). Ejecutar **desde la raíz del repo**
      (`ford my_project.md`): los `-I`/`-M` son relativos al directorio actual. Ya no hace falta
      `pcpp` (preprocesador por defecto de FORD, que no estaba en `environment.yml`); sin Graphviz
      no hay grafos (como antes).
      Mejoras (2026-10-30): `macro: ENABLE_UDTIO` (documenta la salida formateada), `include:
      ./include`, target CMake `doc` (`cmake --build <dir> --target doc` → `<dir>/doc`, corre desde
      la raíz), paso opcional en `check.sh` (si FORD está instalado, debe correr sin errores) y
      docstrings antiguos que FORD leía como metadatos corregidos. Entorno conda `ford`
      (`environment.yml`, sin `gfortran` para no tapar el compilador del sistema): FORD 7.0.13, fypp
      3.3 y Graphviz → grafos de dependencias, 0 warnings, ~40 s.

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
