# Re-análisis 2: código repetido (2026-11-09, rama `fypp/f9-readability`, HEAD `eee4c95`)

Aspecto 2 de 15 del re-análisis posterior al plan R. Compara con `ANALISIS_REPETICION.md` (2026-11-04).
Método: el mismo `tools/dup_report.py` (bloques de ≥ 6 líneas normalizadas que aparecen más de una vez;
números y kinds normalizados, así que es una **cota superior** de las copias textuales), una variante del
mismo detector para Python, y el recuento de lugares que tocó la prueba del recién llegado del aspecto 1 (ley
Perić, `build-snapshots/newcomer_peric.patch`). Todo VERIFICADO; no se cambió código.

## 1. Cifras

| Zona | 2026-11-04: redundantes / líneas de código | Hoy | Lectura |
|---|---|---|---|
| `src/` + `include/` | 1 122 / 9 028 (12 %) | 1 163 / 9 216 (**12.6 %**): 546 en el mismo archivo, 617 entre archivos | Igual. Lo nuevo (kinemáticas a mano, 3.6) agregó ~40 |
| `tests/` | 6 454 / 29 395 (22 %) | 6 290 / 26 306 (**23.9 %**) | Menos líneas (plan R), proporción algo mayor |
| — tests `.f90` heredados (46 archivos) | — | 4 687 / 9 650 (**49 %**) | 75 % de toda la repetición de tests |
| — tests plantilla (64) | — | 1 603 / 16 656 (**9.6 %**) | Razonable |
| Python (`tools/`, `tests/`: 5 146 líneas) | no medido | 0 | Sin bloques repetidos |

El total de líneas de código de tests es más bajo que en 2026-11-04 por el plan R (22 referencias → 3),
aunque los tests unitarios nuevos agregaron líneas.

## 2. Librería (`src` + `include`): focos

| Bloque | Repetición | Estado |
|---|---|---|
| Firma de 6 líneas ×17 en los módulos de operadores de tensores | 96 | Decidido (R5, D-020): firmas a mano |
| `muscle_kin_laws`: la lista de 26 funciones de `muscle_kin_finite` en el `use` y en el `public` | 26 | Nuevo. Es la fachada (re-exporta): lista de nombres, no lógica |
| Cabecera común de `closest_point` y `closest_point_2` (21 + 15 líneas); `use` y argumentos de los procedimientos internos de `perzyna`/`time` (10 ×4); firma de `SOLVE` (7 ×5) | ~95 | **No se hará** (3.3, decisión del usuario): la firma se escribe a mano en cada cuerpo |
| Evaluación conjunta de yield (5 copias de 36 líneas) | ~145 | **No se hará** (3.2): se lee cada criterio completo. El detector ya no la pone entre los primeros 16 bloques porque las copias difieren en nombres |
| Cuerpos 3D/2D a mano en `elasticity_linear` y `visco_frozen` (15–17 líneas ×2, cuatro bloques) | ~65 | Mantener (R3, estilo de la Fase 9) |
| Tipos de tensores entre familias (`4o`/`4o2sym` 19, `2o`/`2osym` 16, `4o2sym`/`4o3sym` 15) | ~50 | R4 sin hacer: diferencias de forma de `vals` y de nombres; ganancia baja |
| `muscle_math_inverses` (16 ×2), `muscle_kin_finite` (15 ×2) | 31 | Variantes por tamaño o medida, escritas a mano a propósito |

**Conclusión:** la repetición de la librería quedó donde la dejaron las decisiones de la Fase 9 y del bloque 3
(explicitud sobre factorización). No apareció repetición nueva evitable.

## 3. Tests: dónde está

Por archivo (redundantes): `test_muscle_math_spectral_derivs.f90` **1 600** (25 % de toda la repetición de
tests), `test_muscle_yield_cpb06.f90` 486, `test_muscle_tensor_4o2sym_kinds.fypp` 265,
`test_muscle_yield_vonmises.f90` 227, `test_muscle_hard_swift.f90` 208, `test_muscle_hard_ludwik.f90` 208,
`test_muscle_tensor_4o_kinds.fypp` 201, `test_muscle_tensor_inverses.f90` 201.

- **Heredados `.f90` (etiqueta `legacy`)**: casos casi idénticos escritos uno tras otro. La regla vigente es no
  moverlos ni reescribirlos (eran la red de seguridad de la migración). **Hecho nuevo:** desde el plan R hay
  tests unitarios con verdad independiente y una herramienta (`mutation_check.py`) que mide qué detecta cada
  grupo. Se puede medir si un heredado aporta algo que los `*_kinds` no detecten, con el mismo criterio que se
  usó para borrar referencias (aspecto 6).
- **Plantillas de tensores** (`test_muscle_tensor_{2o,2osym,4o,4o2sym,4o3sym}_kinds`): bloques de 20–28
  líneas entre familias; parte ya está en `muscle_test_tensor_cases.fypp`.
- **Materiales**: `test_mat_point` y `test_mat_batch` comparten 37 + 24 líneas (preparación de casos); los
  cinco `test_mat_<solver>` comparten 9 líneas de cabecera (ya factorizadas en `mat_sequence_common.fypp`
  salvo la cabecera).

## 4. Repetición de información: un modelo nuevo (prueba Perić)

Una ley `os` nueva tocó **9 archivos y 15 lugares**: la plantilla; `CMakeLists.txt`; el módulo de familia (5:
`use`, `public` y una `interface` por función); `muscle_meta.py` (4: `OVERSTRESS_MODELS`,
`OVERSTRESS_COMPOSED`, nombre para mostrar, `MODEL_PARAMS`); el registro de materiales; `case_params.py`;
el `solver_case` del test de secuencia; el test de la ley; la lista del README.

- La **lista de modelos de la familia** se repite en el módulo de familia y en `muscle_meta.py`
  (`*_MODELS`, `*_COMPOSED`, nombres). La escritura explícita del módulo de familia está decidida (3.7).
- **`MODEL_PARAMS` repite los componentes del tipo Fortran con su valor por defecto y su doc**: es la
  repetición de información más clara que queda (lo usa la tabla de `props` y `mat_init`). La cubre
  `tests/tools/test_muscle_model_params.py` (nombre, orden, tipo y valor por defecto; la doc no).
- El nombre para mostrar y la lista del README son documentación; nada comprueba que estén al día (el
  README de la prueba hubo que editarlo a mano).

## 5. Propuestas (no aplicadas)

- P1. **Hecho (2026-11-10, re-análisis 6 y D-025).** **Medir los tests heredados con `mutation_check`** (en el aspecto 6): por área, quitar del grupo
  unitario los `legacy` y ver qué mutaciones dejan de detectarse. Si ninguna, borrar o reducir los heredados
  con el criterio de D-023 (el usuario decide); `test_muscle_math_spectral_derivs.f90` sólo es 1 600 líneas
  redundantes.
- P2. Test que compare la lista del README y los nombres para mostrar con `*_MODELS` (o generar la lista del
  README en FORD desde `muscle_meta`). Bajo esfuerzo.
- P3. Extender `test_muscle_model_params.py` a la doc de cada parámetro (o tomarla del `!!` del tipo, sin
  duplicarla). A decidir: hoy la doc vive en dos lugares.
- P4. Factorizar la preparación compartida de `test_mat_point`/`test_mat_batch` (~60 líneas). Ganancia baja.
- Sin propuesta para `src`: lo repetido responde a decisiones tomadas (3.2, 3.3, 3.7, R3, R5).
