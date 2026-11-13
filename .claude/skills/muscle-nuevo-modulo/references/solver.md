# Solver de return mapping nuevo

Guía: `pages/contributing.md`, *How to add a solver*, y la tabla «Which one to declare a material with» del
`README.md` (para justificar en qué se diferencia el solver nuevo de los cinco existentes).

## Piezas

1. **Cuerpo** `src/solvers/muscle_<solver>_body.fypp`: no es un módulo; lo incluye la plantilla de material dentro
   de su `#:for k` (procedimientos contenidos del módulo de cada material). Copia el cuerpo más parecido:
   - declara a mano los argumentos de `SOLVE` (la misma firma en todos los cuerpos: `mat, elas, strain, dt,
     state_n, state_np1, stress, tangent, status, iters`), con `dt` aunque no se use;
   - sólo llama a genéricos de familia (`elas_*`, `yield_*`, `hard_*`, `vp_*`, `os_*`): así funciona con todo
     modelo registrado de esas familias;
   - la cabecera `#!` lista las variables que usa del contrato de la plantilla (`k`, `MAT`, `MAT_T`, `ELAS_T`,
     `SOLVE`, `ITER`, `TANGENT`, `S2`, `C4`, `T4`, ...; ver los grupos comentados en
     `src/materials/muscle_material.fypp`);
   - criterios relativos por kind (`NEWTON_TOL_<k>`, `TOL_ELASTIC_<k>`, `SY_MIN_<k>`, `DGAMMA_MIN_<k>` de
     `muscle_kinds`), `iters = min(i, iter_nw)`, `STATUS_NOT_CONVERGED` al agotar iteraciones.
2. **`MATERIAL_SOLVERS`** en `tools/fypp/muscle_meta.py`: familias que llama, `nstate` (7 para el estado plástico),
   `composes` si una familia compone otra, `dt=True` si usa el paso de tiempo, `iter_nw` por defecto.
3. **Cadena de solvers** en `src/materials/muscle_material.fypp`: un `#:elif MAT.solver == '<solver>'` con su
   `#:include` literal (fypp no acepta nombres calculados; `tests/tools/test_muscle_materials.py` comprueba que la
   cadena coincide con `MATERIAL_SOLVERS`). Si comparte la tangente consistente de `closest_point`, inclúyela como
   ellos.
4. **Registro de prueba:** al menos un material por modelo de cada familia del solver (lo exige el test de
   registro), con prefijo nuevo en el nombre; agrega el prefijo a `pages/glossary.md`.
5. **Tests** en `tests/materials/`:
   - `test_muscle_material_<solver>.fypp`: secuencia de deformación con `solver_case(0, ...)` por material (sp vs
     dp, tangente vs diferencias finitas en cada paso plástico);
   - casos en `test_muscle_material_return_mapping.fypp`: retorno radial de von Mises con endurecimiento lineal en
     forma cerrada, consistencia f ≈ 0 en la convergencia, paso elástico bajo el umbral, convergencia cuadrática
     si el solver la promete;
   - `tests/mutations/solvers.toml`: agrega el cuerpo a `sources` si el glob no lo toma.
6. **GPU:** el cuerpo corre en el dispositivo; agrega un material del solver a los casos de
   `tests/gpu/target_cases.py` (`MATERIAL_CASES`) para que tenga test target. Corre `tools/check.sh` completo.

## Cuidado

- Una referencia guardada sólo se justifica para lo que ninguna verdad independiente detecta (el camino de Newton
  con amortiguamiento o búsqueda lineal); requiere decisión del usuario (D-023).
- Si el solver nuevo comparte residuo con otro, la tangente puede ser la misma: compruébalo con diferencias
  finitas, no lo supongas (al jacobiano de `closest_point_2` le falta el factor 2 de los cortes y no se vio hasta el plan R).
