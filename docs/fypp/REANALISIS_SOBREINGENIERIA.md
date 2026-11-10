# Re-análisis 3: sobreingeniería (2026-11-10, rama `fypp/f9-readability`, HEAD `f2acf0a`)

Aspecto 3 de 15 del re-análisis posterior al plan R. Compara con `ANALISIS_SOBREINGENIERIA.md` (2026-11-04).
Mismo criterio: un mecanismo es sobreingeniería si su costo (código, conceptos, tiempo, pasos para quien lo
toca) supera lo que aporta **hoy**; no lo es lo que sirve a un objetivo declarado (sp/dp, 2D/3D, GPU).
Método: consumidores medidos con `grep` y un análisis estático (`ast`) de los 173 nombres de nivel superior de
`muscle_meta.py` y de los campos de sus `namedtuple`. VERIFICADO = medido; los veredictos son juicio mío.

## 1. Lo que señaló el análisis anterior

| # (2026-11-04) | Mecanismo | Entonces | Hoy (VERIFICADO) |
|---|---|---|---|
| 1 | Registro de funciones de familia en Python | Único exceso claro (~250–300 líneas) | **Resuelto** (3.1): `muscle_meta.py` 2 287 → 2 102 líneas; 0 nombres muertos; sólo 3 funciones (24 líneas) las usan únicamente los tests (`model_dims`, `model_params`, `literal_value`). Resto: 2 campos de `ModelInit` (`prm_intent`, `prm_doc`) se asignan y **nadie los lee** |
| 7 | `check.sh` sin modo intermedio | Costo recurrente | **Resuelto**: `--quick` ~1 min |
| 8 | Goldens obligatorios por material | Dudoso: 3 pasos extra por material | **Resuelto** (plan R): `solver_case(0, ...)`, sin registros; la prueba Perić lo confirmó |
| 9 | Proceso documental, 384 `D-0xx` en `src` | Dudoso | `D-0xx` en `src` 78 (resuelto). `docs/fypp/` creció a 2 995 líneas (ver 2.2) |
| 2–6, 11 | `MODEL_PARAMS`, generación de materiales, puntos de entrada, macros, `*_tbp`, soporte sp | Se justifican | Sin cambio: siguen teniendo consumidores |
| 10 | Capacidades sin consumidor interno (2D en modelos, UDTIO, esqueletos FE) | Carga a vigilar | Sin cambio (ver 2.3) |

## 2. Hallazgos nuevos

### 2.1 Infraestructura de referencias dimensionada para 8 600 registros que ahora guarda 83 — **dudoso, costo bajo**
Para 3 archivos (129 líneas, 83 registros) quedan: `muscle_test_reference.f90` (311 líneas, 3 categorías de
tolerancia; las tres se siguen usando), `write_reference.sh` + `merge_reference.py` (83: conservan cabeceras y
registros que no cambian al regenerar), dos `regenerate.sh` (91), 4 tests `*_reference_mutation(_strict)` y
los modos `MUSCLE_REFERENCE_STRICT/REPORT`. Es ~600 líneas y 4 tests de control para 83 registros. Funciona y
no estorba, pero ya no guarda proporción. Los tests de mutación de la comparación tienen sentido mientras
haya referencias; la fusión de cabeceras (`merge_reference.py`) es lo más prescindible, porque estos archivos no
se deberían regenerar salvo por una decisión. Detalle: `write_reference.sh` dice «one reference reference
file» (resto del renombre).

### 2.2 El proceso de la migración sigue siendo el proceso por defecto — **sobreingeniería de proceso**
La migración terminó (Fases 0–9), pero lo que se diseñó para ella sigue activo:
- `CLAUDE.md` describe el ciclo por paso «`test-builder`: tests/valores de referencia contra el código
  **original** → `fypp-migrator` → `check.sh` → `migration-reviewer`», y la regla dura «un PR/commit migra un
  módulo».
- Agentes `test-builder` (su descripción todavía habla de capturar *golden* del código original, contrario a
  D-023), `fypp-migrator` y `migration-reviewer`; skill `migrar-modulo` (155 líneas en total).
- `PLAN.md` (709 líneas) es el plan de una migración cerrada; `DECISIONES.md` 1 294 líneas.
Costo: quien desarrolla hoy (un modelo, un solver) recibe instrucciones de migración que no aplican y una que
contradice D-023. No afecta a un colaborador externo (`.claude/` no se versiona), sí al trabajo con agentes.

### 2.3 Capacidades sin consumidor interno — sin cambio, **vigilar**
- 2D en `elas_linear` y `visco_frozen`: ningún material es 2D (los materiales y solvers son 3D).
- UDTIO (`muscle_tensor_io`, 180 líneas): desactivado por defecto con gfortran < 16.3.
- Esqueletos UMAT/USERMAT/umat41 en `pages/integration.md`: «not compiled nor tested». Son documentación de
  la API principal para usuarios FE y nada comprueba que compilen (aspecto 13).

### 2.4 Tests heredados — **posible exceso, a medir**
Los 46 `.f90` heredados (9 650 líneas de código, 49 % en bloques repetidos; aspecto 2) se mantienen por la
regla «no se mueven ni se reescriben». Desde el plan R hay tests unitarios con verdad independiente que pueden
cubrir lo mismo; si cubren todo, los heredados son redundantes. Se mide en el aspecto 6 con `mutation_check`.

### 2.5 Lo nuevo que se justifica
- `tools/mutation_check.py` (314 líneas) + `tests/mutations/*.toml` (188): se usó para todo el plan R y la guía
  lo pide para código nuevo; la prueba Perić lo usó sin problemas.
- Etiquetas de ctest (`kinds`, `legacy`, `reference`, `target`): pocas líneas de CMake y permiten filtrar.
- 81 módulos de material generados (222 264 líneas en `build-check/`, sp y dp): Fortran no puede declarar un
  módulo por entrada; el registro de prueba cubre cada modelo por solver.

## 3. Conclusión

El único exceso claro del análisis anterior (registro de familias en Python) está resuelto, y también los dos
costos de proceso (`--quick`, goldens por material). **No quedan excesos de código relevantes en la
librería.** Lo que queda es de proceso y de tests: el andamiaje de la migración que sigue activo después de
cerrarla (2.2), una infraestructura de referencias sobredimensionada pero barata (2.1) y la posible
redundancia de los tests heredados (2.4).

## 4. Propuestas (no aplicadas)

| Prioridad | Cambio | Efecto | Riesgo |
|---|---|---|---|
| S1 | **Hecho (2026-11-10)**: `CLAUDE.md` con el ciclo de desarrollo, agentes `test-builder` (sin golden), `muscle-developer` y `reviewer` (renombrados), skill `migrar-modulo` retirada (sin reemplazo, decisión del usuario), `PLAN.md` marcado como histórico. **Pasar el proceso a «desarrollo»**: en `CLAUDE.md`, ciclo por paso para modelos/solvers nuevos (tests unitarios + `mutation_check` → implementación → `check.sh` → revisión); quitar «migra un módulo»; actualizar `test-builder` (sin golden, D-023) y renombrar o retirar `fypp-migrator`/`migrar-modulo`; `PLAN.md` marcado como histórico | Instrucciones coherentes con D-023 y con el estado del proyecto | Bajo; requiere decisión (archivos de reglas del usuario) |
| S2 | Quitar `prm_intent`/`prm_doc` de `ModelInit` | 2 campos muertos menos | Nulo |
| S3 | Corregir «one reference reference file»; dejar la infraestructura de referencias como está (o quitar `merge_reference.py` si se acepta que regenerar reescribe el archivo entero) | Menos superficie | Nulo / bajo |
| S4 | Medir los heredados en el aspecto 6 (P1 del aspecto 2) | Decidir con datos | Nulo |
| S5 | Decidir sobre 2D en modelos, UDTIO y esqueletos FE (que compilen en `check.sh` o se marquen como ejemplo) | Evita mantener sin dueño | Decisión del usuario (aspectos 8 y 13) |
