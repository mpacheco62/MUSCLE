# Re-análisis 1: comprensión (2026-11-09, rama `fypp/f9-readability`, HEAD `c06e900`)

Aspecto 1 de 15 del re-análisis posterior al plan R (lista en `ESTADO_OPTIMO.md`). Compara con
`ANALISIS_COMPRENSION.md` (2026-11-04), hecho antes de los bloques 0–3 y del plan R. Método: las mismas cifras
de entonces, una **prueba del recién llegado** nueva (otra familia: ley de sobretensión Perić en `os`, que
pasa por un solver con referencias) y una **prueba de lectura** (8 preguntas sobre el código respondidas sólo
con la documentación pública). Ambos agentes tuvieron prohibido leer `CLAUDE.md`, `.claude/` y `docs/fypp/`;
los dos avisan que el arnés les inyectó `.claude/rules/` a mitad de camino y que no lo usaron. Lo VERIFICADO
lo comprobé yo en el repo; no se cambió código.

## 1. Cifras (VERIFICADO)

| Medida | 2026-11-04 | Hoy | Lectura |
|---|---|---|---|
| Plantillas de material: líneas / `#:if` de variante | 1 147 / 83 | 1 147 / 83 | Sin cambio: sigue siendo la zona más difícil (prueba de lectura) |
| `muscle_meta.py` | 2 287 líneas, sin índice | 2 102, con índice de secciones | Mejor |
| Referencias `D-0xx` en `src/`+`include/` | 384 | 78 | Meta (≤ ~100) cumplida |
| TODO/FIXME en `src/`, `include/`, `tools/fypp/` | — | 0 | — |
| Páginas públicas (`pages/` + README) | sin glosario ni arquitectura | 10 archivos, 1 330 líneas (glosario, arquitectura, convenciones) | Mejor |
| Tests: plantillas / `.f90` / Python | 74 / 49 / — | 64 / 46 / 10 | Menos archivos |
| Test más grande | `test_yield_kinds` 1 057 | `test_math_spectral_derivs.f90` 1 354 (legacy), `test_yield_kinds` 1 185 | **Los monolitos crecieron** con el plan R |
| Procedimiento más largo en `src` | — | `UPDATE` del kernel 136 líneas, `yld2004` 133, cuerpo `closest_point_perzyna` 114 | Aceptable para el dominio |
| Restos de «golden» | 933 menciones | 13, todas en `tests/` (comentarios «Cases of test_*_golden», función `_golden` de `props_samples.py`) más la nota histórica del glosario | Casi limpio |

## 2. Prueba del recién llegado (VERIFICADO: parche en `build-snapshots/newcomer_peric.patch`)

Tarea: ley Perić `g = (1/η)[(σeq/σy)^(1/m) − 1]` en la familia `os`, con material de prueba, tests y
`ctest` completo. **Resultado: completo** (antes: parcial, el 2026-11-04; lograda con la guía corregida el
2026-11-08). 169/169, `check.sh --quick` OK, material `pz_peric_vm_swift` con `solver_case(0, ...)` (sin
registros), `mutation_check` del archivo nuevo: 23/24 detectadas (la que sobrevive, `r <= 1` → `r < 1`, es
equivalente: g(1) = 0 por ambas ramas). ~30 comandos, 9 archivos, 1 build fallido (error propio).

Fricciones, por impacto (las tres primeras VERIFICADAS en `pages/contributing.md`):
1. **Nombres de archivo y módulo de la familia `os`**: la guía dice `muscle_<prefix>_<name>.fypp`
   (`contributing.md:166`); las leyes `os` viven en `muscle_vp_<name>.fypp` y el módulo de familia es
   `muscle_vp_overstress_laws`. Falta una tabla familia → prefijo → archivo → módulo de familia.
2. **Listas `*_COMPOSED`** (`OVERSTRESS_COMPOSED`, `VISCOPLASTIC_COMPOSED` en `muscle_meta.py:106,115`): la
   guía no las menciona; el recién llegado adivinó bien.
3. **Área de mutaciones**: la guía dice `tools/mutation_check.py <family>`; el área de `os` es `vp`. No dice
   cómo limitarlo al archivo nuevo (`--list` + `--only`).
4. Qué solvers usan cada familia: la guía sólo lo dice para `hard`; para `os` hubo que leer
   `MATERIAL_SOLVERS`.
5. `solver_case`: la guía no da el orden de los campos (los comentarios de `mat_sequence_common.fypp` sí).
6. El ejemplo de test sirve para `hard`; para `os` copió los de Perzyna/Chaboche. Los ids de `record()`
   (`10*i + n`) se agotaron en 9 leyes.
7. No dice si una ley `os` necesita validador `MODEL_INITS` (Perzyna no tiene) ni si un modelo nuevo va en
   `tests/gpu/target_cases.py`.
8. Agregar el modelo a la lista del README no está en la lista de verificación.
9. El lint no aplica el límite de 132 columnas a las líneas de plantilla.

## 3. Prueba de lectura

Respondió bien las 8 preguntas (recorrido de la API de punto hasta la tangente, registro, prefijos, kinds y
2D/3D, GPU, unitario vs referencia, `.approx.`), en 2–9 archivos por pregunta. Las páginas que más sirvieron:
`architecture.md`, `glossary.md`, `contributing.md`. Hallazgos (VERIFICADOS):

1. **Dos formatos de historia que se contradicen.** `muscle_plastic_state` documenta
   `plastic_state_pack/unpack` como «el estado en el array de historia FE» con 13 componentes
   (tensión, deformación plástica, εp); los materiales empaquetan con su propio `PACK` (kernel, l. 315):
   deformación plástica y εp (7, más el estado viscoelástico), sin tensión. El pack del módulo sólo lo usa su
   test (`tests/plasticity/test_plastic_state_kinds.fypp`); el tipo sí lo usan los materiales.
2. **«plan R»** aparece en `glossary.md:65` y `conventions.md:94` sin definirse; `mutation_check.py --help` y
   algunos comentarios citan `docs/fypp/` y `D-0xx`, que un externo no tiene.
3. **Tolerancia de `.approx.`**: `conventions.md:90` la da como 1e-12 (dp) sin decir que 4O2sym (1e-8) y las
   identidades escaladas (1e-10) tienen otra por defecto; el trato de Inf/NaN sólo está en el código.
4. «Only three remain» (`contributing.md:217`) frente a 7 tests con etiqueta `reference` (3 + 4 de
   mutación): correcto pero confunde.
5. No hay página que diga cuándo elegir `cp` o `cp2` (sólo que existen).
6. Lo más difícil de leer: `muscle_material.fypp` (~25 `#:set` antes del primer Fortran), las cuatro variantes
   de `UPDATE` en el kernel, `FIN`/`XS` definidos distinto en `point` y `batch`, `muscle_meta.py`, y
   `case_params.py` con conjuntos codificados en cadenas y `solver_case` posicional.

## 4. Conclusión

La comprensión **mejoró de forma medible**: la prueba del recién llegado pasó de parcial a completa en una
familia distinta de la que se usó para corregir la guía, y la prueba de lectura respondió todo con la
documentación pública. Lo que queda son **huecos de la guía para familias distintas de `hard`** (la guía se
escribió con `hard` como ejemplo) y algunas contradicciones puntuales. La generación de materiales sigue
siendo la parte opaca (83 `#:if`, igual que antes); la alternativa D sigue archivada.

## 5. Propuestas (no aplicadas; esperan decisión)

**Bajo esfuerzo (≈ 0.5 día en total)**
- C1. **Hecho (2026-11-10)**: tabla por familia, `--list`/`--only` del chequeo de mutaciones, README en la lista de verificación. Propuesta original: `contributing.md`: tabla por familia (prefijo, archivo, módulo de familia, listas de `muscle_meta.py`
  incl. `*_COMPOSED`, área de mutaciones, solvers que la usan, test donde va la ley); `--only` para el chequeo
  de mutaciones; orden de campos de `solver_case`; README y `target_cases.py` en la lista de verificación.
- C2. **Hecho (2026-11-09, D-024, opción (a) del usuario):** se borraron `plastic_state_pack/unpack` y
  `PLASTIC_STATE_SIZE`; queda el tipo `plastic_state_<k>`. Las otras opciones eran documentarlo como formato
  alternativo o que los materiales lo usaran (cambio de layout de `hsv`).
- C3. **Hecho (2026-11-10).** Glosario: entrada «plan R» (o quitar la mención); `conventions.md`: tolerancias por familia y trato de
  Inf/NaN de `.approx.`; aclarar «three reference tests (+4 mutation checks)».
- C4. Restos de «golden» en comentarios de tests y `_golden` de `props_samples.py`.
- C5. **Hecho (2026-11-10)**, tabla en el README. Una línea en `glossary.md` o `architecture.md` sobre cuándo usar `cp` o `cp2`.

**Esfuerzo medio / a decidir**
- C6. Partir los monolitos de tests que crecieron (`test_yield_kinds` 1 185, `test_vp_laws_kinds` 951) o
  dejarlos (el recién llegado agregó su test sin problema, como en 3.5).
- C7. Lint del límite de columnas en plantillas.
- C8. Legibilidad de la plantilla de material (`#:set` agrupados y comentados, `FIN`/`XS` con un solo
  significado): sólo comentarios y nombres, snapshot idéntico.
