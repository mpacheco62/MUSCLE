# Re-análisis 5: ciclo de desarrollo (2026-11-10, rama `fypp/f9-readability`, HEAD `e8d7690`)

Aspecto 5 de 15 del re-análisis posterior al plan R. Compara con la sección 2.5 de `ESTADO_OPTIMO.md`
(2026-11-08). Método: cada operación sobre un build Debug limpio (`build-maint-Debug`, gfortran 15.2, Ninja,
`-j8`, `sp;dp`, 169 tests), con los archivos restaurados después de cada prueba; el ciclo «modelo nuevo» aplica
el parche de la prueba Perić (`build-snapshots/newcomer_peric.patch`) y lo revierte. Todo VERIFICADO.

## 1. Tiempos

| Operación | 2026-11-08 | Hoy |
|---|---|---|
| Build limpio Debug / Release (con tests y los 81 materiales) | — | 71 s / 68 s |
| Build sin cambios | — | 0.3 s |
| `ctest` completo, CPU | 9.6 s (181 tests) | **8.9 s** (169 tests) |
| Tocar una plantilla sin cambiar su salida | ~1.3 s | 1.3 s (fypp regenera, nada recompila) |
| Cambiar el cuerpo de una ley (fórmula) | — | **1.1 s** (sólo su módulo; el `.mod` no cambia y nada más recompila) |
| Agregar algo no re-exportado a una ley | — | 1.1 s |
| Cambiar la **interfaz** de una ley de endurecimiento (argumento nuevo) | «recompila 76–80 materiales» | **6 s**: 85 objetos, 76 materiales |
| Agregar una constante pública a `muscle_kinds` | 479 pasos, ~60 s | 457 pasos, 64 s (reconstrucción casi total) |
| Cambiar `include/muscle_common.fypp` (aunque sea un comentario `#!`) | recompila 142 de 144 objetos | **reconstrucción total**: 1 062 pasos, 59 s |
| Cambiar `muscle_meta.py` (aunque sea un comentario) | relanza CMake, 142/144 | **reconstrucción total**: relanza CMake, 1 101 pasos, 58 s |
| **Modelo nuevo completo** (parche Perić: ley, familia, `muscle_meta.py`, material, tests) | — | build 60 s + `ctest` 10 s ≈ **70 s** |
| `mutation_check.py hard` (30 mutaciones, build dp propio) | — | 66 s (~2 s por mutación), `CRITERION: OK` |
| `check.sh --quick` | 57 s (medido al crearlo) | 75 s |
| `check.sh` completo (UDTIO, FORD, 3 builds GPU, tests en la Intel, auditoría) | decenas de minutos | ~13 min (corridas del plan R) |

## 2. Lectura

- **El bucle interno es rápido**: editar la fórmula de un modelo y correr sus tests cuesta ~1 s de build más
  el test (`ctest -R`); todos los tests, 9 s. El aislamiento por módulo funciona: un cambio de cuerpo no
  propaga porque gfortran no reescribe un `.mod` igual y Ninja no recompila a los que dependen de él.
- **El cuello de botella es `muscle_meta.py`**: cualquier cambio (incluso un comentario) relanza CMake y
  reconstruye todo (~60 s), porque todas las plantillas declaran `muscle_meta.py` y `muscle_common.fypp` como
  dependencia (`MUSCLE_FYPP_DEPENDS`). Un modelo nuevo siempre lo toca (`MODEL_PARAMS`, `*_MODELS`), así que su
  primer build es siempre total. Empeoró algo respecto del 2026-11-08 (antes 142/144 objetos de la
  librería; hoy también los tests).
- **`check.sh --quick` subió de 57 a 75 s** desde su creación (más tests unitarios del plan R y el lint).
  Sigue siendo el ciclo de verificación diario razonable.
- **El chequeo de mutaciones es caro si se pide todo**: ~2 s por mutación; un área completa con
  `--samples all` (hasta 1 286 en yield) es de 30–45 min. La guía lo pide para código nuevo: conviene limitarlo
  al archivo nuevo con `--only` (la prueba Perić: 24 mutaciones, ~1 min), lo que la guía no explica (aspecto 1,
  C1).
- `muscle_kinds` sigue siendo una cascada casi total (1.2 descartado por el usuario: 89 de 133 módulos usan
  alguna constante).

## 3. Propuestas (no aplicadas)

| Prioridad | Cambio | Efecto | Riesgo |
|---|---|---|---|
| K1 | **Hecho (2026-11-10)**, sección «Development cycle» de `contributing.md`. `contributing.md`: tiempos esperados y el ciclo recomendado (`ctest -R` mientras se edita, `--quick` antes de cada commit, completo al cerrar una etapa; `mutation_check --only` sobre el archivo nuevo) | Expectativas claras; evita correr un área completa | Nulo (junto con C1) |
| K2 | Dependencias más finas de fypp: que sólo las plantillas que leen el registro de materiales o los parámetros dependan de `muscle_meta.py` (o separar los datos que cambian con cada modelo en un archivo propio) | Un modelo nuevo dejaría de reconstruir los tensores y la math (~60 s → estimado 20–30 s) | Medio: una dependencia faltante deja generados viejos sin aviso; requiere un test de dependencias. **Beneficio dudoso para 60 s** |
| K3 | Nada que hacer en el bucle interno ni en `ctest` | — | — |

Recomendación: K1 sí (barato, junto con C1); K2 no por ahora: el costo total (~1 min) no justifica el riesgo de
dependencias incompletas.
