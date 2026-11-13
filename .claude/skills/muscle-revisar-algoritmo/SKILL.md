---
name: muscle-revisar-algoritmo
description: Revisa un algoritmo numérico de MUSCLE (librería Fortran 2008 + FYPP de modelos constitutivos, repo MUSCLE_FYPP) - una rutina de math (autovalores, inversas, sistemas, funciones de matrices, derivadas espectrales), un return mapping, una ley o un kernel de material - en precisión (contra una verdad independiente de alta precisión), robustez (casos borde, escalas extremas, excepciones de punto flotante), calidad de sus tests (chequeo de mutaciones), velocidad en CPU y en GPU, tamaño del kernel en la GPU (memoria privada, spill, stack calls, FP64 emulado) y riesgos de dispositivo, compara alternativas y entrega un informe con lo que conviene mejorar. Úsala siempre que en MUSCLE se pida revisar, auditar, medir, perfilar, optimizar, comparar o validar un algoritmo o una rutina numérica (aunque no diga "revisar", p. ej. "¿está bien implementado el Newton de closest point?", "¿por qué es lento este kernel en GPU?", "compara dos formas de calcular la polar", "¿conviene cambiar el algoritmo de autovalores?", "mide cuánto cuesta esta ley").
---

# Revisar un algoritmo numérico de MUSCLE

Esta skill recoge el método con que el proyecto revisó sus algoritmos (plan espectral, `docs/fypp/PLAN_ESPECTRAL.md`
y `docs/fypp/ESPECTRAL_F0.md`; errores de GPU de D-030 y D-031): medir antes de opinar, contra una verdad que no sale
del código revisado, y separar lo que es del algoritmo de lo que es del compilador o del hardware. La fuente de
verdad es el repositorio: si algo de aquí contradice `pages/contributing.md`, `pages/math.md`,
`docs/fypp/DECISIONES.md` o, si existen en tu copia, `CLAUDE.md` y `.claude/rules/`, mandan ellos.

**Revisar no es cambiar.** La skill produce mediciones y un informe con recomendaciones; cambiar el algoritmo, un
umbral o el comportamiento numérico de algo existente necesita una decisión del usuario (registrada en
`docs/fypp/DECISIONES.md`). Los problemas marcados «Marcado, sin corregir» en `docs/fypp/PROBLEMAS.md` no se tocan.

## 0. Delimitar la revisión

1. Trabaja en el repo de MUSCLE que contiene esta skill; para experimentar con variantes usa un **worktree** aparte
   (`git worktree add --detach <dir> HEAD`), nunca el árbol del usuario. Si el scratchpad está en un tmpfs (`/tmp`),
   vigila el espacio (`df -h /tmp`): los builds de GPU de varias variantes lo llenan y el compilador falla con
   errores crípticos.
2. Identifica la rutina, su firma y **sus consumidores** (`grep -rn <nombre> src`): la revisión se mide en los
   kernels que la usan, no sólo aislada. Lee su documentación FORD, su decisión (`D-0xx`) y las filas de
   `PROBLEMAS.md` que la nombran.
3. Pregunta al usuario qué le importa si no está claro (precisión, CPU, GPU, tamaño del kernel) y el presupuesto
   aceptable (p. ej. «el kernel no crece más de 1.2x»): las recomendaciones dependen de eso.

## 1. Precisión y robustez → `references/precision.md`

Verdad independiente en alta precisión (cuádruple o `mpmath`), métricas en unidades de ε del kind con la escala
correcta, un catálogo de casos difíciles (brechas pequeñas, casi degenerados, escalas extremas, ceros, NaN/Inf) y las
excepciones de punto flotante con trampas. Distingue el error del algoritmo del límite del problema (condición):
una pérdida ~ε/brecha que también tendría cualquier algoritmo no es un defecto.

## 2. Tests y mutaciones

1. Lista los tests que ejercitan la rutina y qué verdad compara cada uno (forma cerrada, propiedad, otra ruta,
   cuádruple). Un test que compara contra la salida guardada del propio código no cuenta como verdad.
2. Chequeo de mutaciones sobre las líneas de la rutina:
   `python3 tools/mutation_check.py <área> --samples all --only <ids>` (los ids con
   `.claude/skills/muscle-nuevo-modulo/scripts/mutation_ids.sh <repo> <área> <archivo>`, o filtrando
   `--list --samples all` por rango de líneas). No lo corras a la vez que otra cosa que lea `src/`: modifica el
   árbol temporalmente.
3. Clasifica cada sobreviviente con una razón **comprobada**: equivalente exacta (no puede cambiar ningún
   resultado), equivalente por debajo del redondeo (cambia ~ε), sólo en empates, rama inalcanzable, o **hueco de
   test**. Cuidado con las equivalencias «matemáticas» que no lo son en coma flotante: la invariancia ante un
   desplazamiento sólo vale si el desplazamiento es del tamaño de lo que se descompone (una media 1e6 veces mayor
   que el desviador sí cambia el resultado), y una rama «inalcanzable» puede ser la que propaga un NaN. Si dudas,
   compila el mutante y pruébalo con un caso dirigido.
4. Las que sólo mata `math_reference` u otra referencia mínima son esperables (plan R, D-023): dilo, no las
   borres.

## 3. Velocidad en CPU y en GPU → `references/rendimiento.md`

Costo por llamada (micro-banco) y costo del kernel completo de cada consumidor, en CPU (un hilo, Release, fijado a
un núcleo, variantes intercaladas, mínimo de 7 rondas, con un caso de control que no use la rutina) y en la GPU
Intel (Release, el banco de lote `bench_batch_gpu_<material>_<kind>` con su verificación GPU = CPU). Las corridas de
GPU van **de a una**: dos procesos en la GPU a la vez producen pérdidas del dispositivo falsas.

## 4. Tamaño del kernel y riesgos de dispositivo → `references/kernel_gpu.md`

`scripts/igc_kernel_report.sh <ejecutable> [regex]` (esta skill) da por kernel la memoria privada, el spill, las
stack calls, las llamadas y las llamadas a la emulación FP64, y los `alloca` dinámicos de la imagen;
`tools/kernel_resources.sh` da GRF y spill por material (SPIR64, AMDGCN). La referencia lista los riesgos que el
proyecto ya pagó (D-019.5c/5f/8, D-030, D-031) y las variables de IGC que sirven para diagnosticar, no para
producción.

## 5. Comparar alternativas

Si la revisión sugiere otro algoritmo, mídelo antes de recomendarlo, con el mismo arnés:
1. escribe los candidatos como procedimientos `pure`, `@:DEVICE()`, en un archivo de `benchmarks/` (precedente:
   `benchmarks/math/spectral_candidates.fypp`), incluido por el banco de precisión;
2. para el costo en los kernels reales, en un worktree reemplaza el cuerpo de la rutina por el candidato con un
   define de FYPP (`-DMUSCLE_FYPP_DEFINES=...`, precedente: `SPECTRAL_PROTO` de F0) y corre los mismos bancos y la
   suite de tests con cada variante: qué tests fallan y por qué (diferencias legítimas, como otra base de un
   autoespacio repetido, o errores);
3. resume cada candidato en precisión, costo CPU/GPU, tamaño del kernel y tests, frente al presupuesto del usuario.

## 6. Informe → `references/informe.md`

Un informe con tablas (precisión, costo, GPU, kernel, mutaciones) y recomendaciones priorizadas, cada una con su
evidencia y su costo, separando lo comprobado de lo que es hipótesis. Si el usuario lo quiere guardado, va a
`docs/fypp/REVISION_<rutina>.md` (español) y los defectos nuevos a `PROBLEMAS.md`; nada de cambios en `src/` sin
decisión.
