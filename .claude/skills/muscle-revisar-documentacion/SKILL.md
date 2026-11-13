---
name: muscle-revisar-documentacion
description: Revisa que la documentación de MUSCLE (librería Fortran 2008 + FYPP de modelos constitutivos, repo MUSCLE_FYPP) esté al día con lo implementado - instala FORD si falta (con permiso), genera la documentación FORD y la valida, detecta deriva entre docs y código (nombres, rutas, enlaces y decisiones D-0xx que ya no existen, módulos o procedimientos públicos sin `!!`, modelos sin mencionar) y revisa a mano lo que un script no ve (que los textos digan lo que el código hace hoy, tablas de pages/, README, decisiones y problemas). Úsala siempre que en MUSCLE se pida revisar, actualizar, generar o validar la documentación, antes de un PR o de cerrar una etapa, después de cambiar una API, un modelo, un algoritmo o una regla (aunque no diga "documentación", p. ej. "¿está al día el README?", "genera la doc con FORD", "revisa que pages/ refleje el cambio", "¿quedó algo desactualizado tras el refactor?").
---

# Revisar la documentación de MUSCLE frente a lo implementado

El código cambia y la documentación se queda atrás sin que nadie lo note: un nombre retirado sigue en una página,
una tabla omite el modelo nuevo, un comentario FORD describe el algoritmo anterior. Esta skill lo detecta con dos
scripts y una revisión dirigida, y corrige la documentación (no el código). La fuente de verdad es el código y
`docs/fypp/DECISIONES.md`; si la documentación y el código discrepan, se corrige la documentación salvo que el
código tenga un error, que se reporta al usuario en lugar de «documentarlo».

## 0. Alcance

- **Cambio acotado** (lo normal, antes de un PR): `git diff --stat dev...HEAD` (o el rango que diga el usuario) da
  los archivos tocados; la revisión manual se concentra en ellos y en lo que los cita. Los scripts corren siempre
  sobre todo el repo (son baratos).
- **Revisión completa** (cierre de etapa): todo `pages/`, `README.md`, `my_project.md` y las skills.

## 1. FORD: generar y validar

`.claude/skills/muscle-revisar-documentacion/scripts/ford_build.sh <repo> [--warn] [OUTPUT_DIR]` busca FORD en el
`PATH`, en el entorno conda `ford` (el de `environment.yml`) o en `build-check/.ford-venv`, genera la documentación
desde la raíz (`my_project.md` preprocesa las plantillas con fypp) y resume los avisos sin el ruido conocido.

- **Si no encuentra FORD, pide permiso al usuario antes de instalar** y luego corre con `--install` (crea el
  entorno conda desde `environment.yml`, o un venv con `pip install ford fypp`; necesita red). Nunca instales en
  el Python del sistema.
- Debe terminar con código 0 y sin errores. Tarda ~1.5 min.
- Con `--warn`: hoy quedan ~91 interfaces genéricas, ~42 procedimientos ligados a tipo (`*_tbp` de la API
  histórica) y ~6 procedimientos sin `!!` (línea base 2026-11-13). Lo **nuevo** de tu cambio no debe sumar
  avisos: documenta lo público que agregaste (interfaz genérica incluida).
- El script imprime el `PATH` con el que `tools/check.sh` también genera la documentación (si FORD no está en el
  `PATH`, `check.sh` la salta y dice «FORD not installed»).
- Abre unas páginas de `<OUTPUT_DIR>/index.html` (la del módulo cambiado, `page/`): fórmulas `\( \)` que no se
  renderizan, listas rotas, tablas, enlaces.

## 2. Deriva automática

`python3 .claude/skills/muscle-revisar-documentacion/scripts/doc_drift.py <repo> [--only <chequeo>]`:

| Chequeo | Qué detecta |
|---|---|
| `undocumented` | módulos, tipos y procedimientos públicos de `src/` sin su línea `!!` |
| `names` | nombres en backticks de `pages/`, `README.md`, `my_project.md` y las skills que no están en el código (renombrados o retirados) |
| `paths` | rutas citadas que no existen |
| `links` | enlaces `](page/<x>.html)` a páginas que no existen |
| `decisions` | referencias `D-0NN` a decisiones que no están en `DECISIONES.md` |
| `models` | modelos de una familia que la documentación de usuario nunca menciona |

Objetivo: 0 hallazgos. Cada uno se corrige o, si es intencional (un nombre de ejemplo o uno retirado que
`pages/migration.md` documenta a propósito), se agrega a `EXAMPLES` o `RETIRED` del script con un comentario. Los
`docs/fypp/` históricos (planes, análisis) no se revisan por nombres: describen el pasado a propósito.

## 3. Revisión dirigida (lo que un script no ve) → `references/mapa.md`

Para cada cosa que cambió, el mapa dice qué documentos revisar. Lee lo que dicen y compáralo con el código:
- **FORD del cambio**: que el `!!` del módulo y de cada procedimiento describa el algoritmo **actual**, sus
  argumentos con unidades y dominio, el `status`, los kinds, las referencias; que no cite rutinas, umbrales o
  decisiones reemplazadas.
- **Páginas y tablas**: `README.md` y `my_project.md` (listas de modelos y capacidades), `pages/contributing.md`
  (tabla por familia, checklists), `pages/math.md` («necesito X → módulo Y»), `pages/conventions.md`,
  `pages/integration.md` y `pages/migration.md` (API pública), `pages/glossary.md` (términos y prefijos).
- **Decisiones y estado**: `docs/fypp/DECISIONES.md` (decisión nueva o nota de reemplazo en la que cambió),
  `PROBLEMAS.md` (estado de las filas afectadas), `ESTADO_OPTIMO.md` §7 (pendientes); si existe en la copia
  local, la línea base de `CLAUDE.md`.
- **Skills** (`.claude/skills/`): reglas, rutas o decisiones que citan.

## 4. Corregir y cerrar

- Corrige la documentación en el mismo PR que el cambio de código (regla de `pages/contributing.md`). Si el
  cambio de código ya está hecho, un commit `docs(...)` aparte.
- Verifica: `doc_drift.py` en 0, `ford_build.sh` con código 0, y `tools/check.sh --lint-only` (líneas ≤ 132
  columnas en las fuentes, SPDX).
- Informe breve: hallazgos de cada chequeo (corregidos o justificados), avisos de FORD frente a la línea base,
  lo revisado a mano y lo que queda para el usuario (p. ej. documentación que exige decidir algo del código).
