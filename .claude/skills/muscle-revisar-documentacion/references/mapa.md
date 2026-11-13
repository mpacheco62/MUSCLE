# Qué cambió → qué documentación revisar

| Cambio | Revisar |
|---|---|
| Modelo nuevo de una familia | FORD del módulo (fórmula, parámetros con unidades, referencia); lista de modelos de `README.md` y de `my_project.md` («Key Features»); tabla por familia de `pages/contributing.md` si cambia algo de ella; `pages/glossary.md` si hay prefijo o término nuevo; `tools/fypp/muscle_materials.py` si se registra un material de prueba |
| Material nuevo o solver nuevo | `pages/integration.md` (API de punto y de lote), `pages/architecture.md`, `pages/contributing.md` (*How to add a solver*, qué solver declarar), `README.md` |
| Cambio de API pública (nombre, firma, módulo, constante) | FORD de lo cambiado; `pages/migration.md` (equivalencias); `pages/integration.md`; `pages/conventions.md`; `pages/math.md` o `pages/tensors.md` según el área; toda mención en `README.md`/`pages/` (`doc_drift.py --only names`); decisión en `DECISIONES.md` (sin alias: D-027) |
| Algoritmo cambiado (misma API) | FORD del procedimiento (algoritmo, precisión, casos borde, cuándo falla); `pages/math.md` si es de math; nota de reemplazo en la decisión anterior y decisión nueva; `PROBLEMAS.md` (filas que el cambio corrige o crea) |
| Umbral, tolerancia o constante de `muscle_kinds` | FORD de la constante; `pages/conventions.md`; tolerancias en `pages/tensors.md` y en las cabeceras de los tests que la citan |
| Regla nueva (lint, GPU, estilo) | `pages/contributing.md` (checklists, *Development cycle*); si existen localmente, `.claude/rules/*.md` y `CLAUDE.md`; las skills de `.claude/skills/` que citan la regla |
| Error de compilador o de GPU rodeado | `PROBLEMAS.md` (fila del problema y registro de soluciones provisionales vigentes, con qué lo verifica y cuándo retirarla); `tools/gpu_skip_tests.sh` si cambia lo excluido |
| Test o verificación nuevos | cabecera del test (qué verdad compara); `pages/contributing.md` (*Tests*); `CLAUDE.md` local (línea base de tests) |
| Herramienta nueva (`tools/`, scripts de skills) | `pages/contributing.md`; `CLAUDE.md` local (*Build y tests*) |
| Fin de una etapa | `docs/fypp/ESTADO_OPTIMO.md` (§7 pendientes), `PROBLEMAS.md` completo, `README.md` (toolchain probado), la línea base de tests |

## Qué no se cambia

- Los documentos históricos de `docs/fypp/` (`PLAN.md`, análisis y re-análisis, informes de fases): se les agrega
  una nota «superado por D-0xx» si confunden, pero no se reescriben.
- `pages/migration.md` nombra la API retirada a propósito (columna «antes»).
- El código: si la documentación revela un error del código, se reporta al usuario; no se «documenta» el error.
