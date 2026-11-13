---
name: muscle-nuevo-modulo
description: Crea un módulo nuevo en MUSCLE (librería Fortran 2008 + FYPP de modelos constitutivos, repo MUSCLE_FYPP) siguiendo las reglas del proyecto - un modelo nuevo de una familia existente (ley de endurecimiento, ley viscoplástica o de sobretensión, criterio de fluencia, elasticidad, viscoelasticidad, cinemática), un solver de return mapping, una familia de modelos nueva o un módulo de soporte (math, tensores). Úsala siempre que en MUSCLE se pida agregar, implementar, portar o escribir una ley, criterio, modelo, solver, familia o módulo nuevo (aunque no diga "módulo", p. ej. "agrega la ley de Perić", "implementa Hosford", "quiero un solver de Euler implícito", "pasa este modelo de Abaqus a MUSCLE"), y también para revisar si un módulo nuevo cumple las reglas.
---

# Crear un módulo nuevo en MUSCLE

Esta skill encadena lo que el proyecto aprendió en sus análisis (comprensión, repetición, mantenibilidad,
calidad de tests, convenciones) para que un módulo nuevo entre ya cumpliendo las reglas, no después de una
ronda de correcciones. La fuente de verdad es el repositorio: si algo de aquí contradice `pages/contributing.md`,
`pages/conventions.md`, `docs/fypp/DECISIONES.md` o, si existen en tu copia (no se versionan), `CLAUDE.md` y
`.claude/rules/`, mandan ellos (y conviene
avisar al usuario de la diferencia).

## 0. Ubicarse y decidir el caso

1. Trabaja en el repo de MUSCLE que contiene esta skill (`.claude/skills/muscle-nuevo-modulo/`; los scripts se
   llaman desde la raíz como `.claude/skills/muscle-nuevo-modulo/scripts/<script>`), en una rama nueva desde `dev`. Lee
   `CLAUDE.md` y `pages/contributing.md` (secciones *How to add a model*, *How to add a solver*, *Checklists*,
   *Tests*, *Development cycle*). No hagas `git push` nunca; commits locales sólo si el usuario los pide, sin
   autoría de Claude en el mensaje.
2. Identifica el caso y lee **sólo** su referencia:

| Caso | Ejemplo | Referencia |
|---|---|---|
| Modelo de una familia existente (`hard`, `vp`, `os`, `yield`, `elas`, `visco`, `kin`) | ley de Perić, criterio de Hosford | `references/modelo.md` |
| Solver de return mapping | Euler implícito con subincrementos | `references/solver.md` |
| Familia de modelos nueva | daño, endurecimiento cinemático | `references/familia.md` |
| Módulo de soporte (math, tensores, utilidades de kernel) | descomposición polar, nueva inversa | `references/soporte.md` |

3. **Antes de escribir nada que no sea la física del módulo, busca si ya existe.** Corre
   `.claude/skills/muscle-nuevo-modulo/scripts/api_index.py <repo> [palabra ...]` (genéricos, tipos y constantes de toda la librería con su
   documentación; busca por concepto: `eigen`, `deviatoric`, `inverse`, `rotation`, `tolerance`) y lee
   `references/ubicacion.md` (bloques que no se reescriben y dónde va cada cosa nueva). Reescribir el desviador, un
   determinante o una tolerancia que ya está en `muscle_kinds` crea dos versiones que divergen; poner álgebra general
   dentro de un modelo la esconde de los demás modelos y de sus tests. Lo que falte y no sea del modelo va a su área
   (`src/math`, `src/tensors`, un módulo de soporte de la familia): eso es un segundo módulo, de soporte, con su
   propia referencia.
4. Si falta una decisión de diseño (un parámetro con dos convenciones posibles, un cambio de API pública, un
   valor por defecto discutible, validar o no parámetros), **pregunta al usuario con opciones y una
   recomendación** antes de escribir código, y registra la respuesta en `docs/fypp/DECISIONES.md` (siguiente
   `D-0xx`). Los analistas del proyecto encontraron que lo que más cuesta después son las decisiones tomadas
   en silencio.

## 1. Tests primero, con verdad independiente

Escribe los tests antes que el módulo (o con el agente `test-builder`, si existe en tu copia de `.claude/agents/`). Un test útil compara contra algo que
**no sale del código que prueba**: forma cerrada, propiedad (simetría, homogeneidad, normalización, objetividad,
`A A⁻¹ = I`, límite isótropo), derivadas contra diferencias finitas (`tests/common/muscle_fd.fypp`), la misma
cantidad por otra ruta, o un valor publicado. Copiar la salida actual como «esperado» no detecta errores: el
re-análisis de tests encontró que esos eran justo los huecos (p. ej. un factor de normalización mal escrito que
sólo reescala, invisible a diferencias finitas y a homogeneidad).

- Ambos kinds: plantilla `.fypp` con `#:for k in REAL_KINDS`, tolerancias de `muscle_kinds` (`TOL_APPROX_<k>`)
  o justificadas en el test; entradas representables en `sp`.
- Casos borde que el código promete: valores nulos, pisos y recortes, entradas inválidas (y su `status`),
  identidad, escalas extremas (homogeneidad con 2^±20 es exacta).
- Varía los parámetros de prueba: un parámetro que en todos los casos vale 1 no prueba nada de su factor.
- Nombres: el test de `src/<ruta>/muscle_<x>.fypp` es `tests/<ruta>/test_muscle_<x>[_<tema>].fypp` (D-025),
  registrado con `add_muscle_fypp_test(<nombre> <archivo>)` en el `CMakeLists.txt` espejo.
- **No crees tests de referencia** (valores guardados en `tests/reference/`) salvo decisión registrada (D-023):
  sólo existen donde no hay verdad independiente (el camino de Newton de dos solvers).

## 2. Escribir el módulo

- **FYPP mínimo (D-020):** la plantilla se lee como Fortran. Tipos, firmas con `intent`, cuerpos y cierres se
  escriben a mano; FYPP sólo instancia por kind/dimensión, cambia interruptores de build y lista nombres. Nada de
  macros que armen declaraciones desde Python. Copia el modelo vecino de la familia y cambia la física.
- **Directivas FYPP (D-026):** `#:if`/`#:for`/`#:set` indentados al nivel del código que controlan, sin sangrar
  más el contenido; `#:include`, `#:def` y la guarda de archivo en columna 0.
- **Kernel GPU:** todo procedimiento por punto es `pure`, con `implicit none` y `@:DEVICE()` justo después, sin
  `class`, `allocatable`, punteros, I/O, `stop`, LAPACK, `save`, recursión ni módulos IEEE; errores por un entero
  `status` (`STATUS_*` de `muscle_status`). Además:
  - ni `scale()` ni `exponent()`: `fp_scale` y `fp_exponent` de `muscle_math_fp` (mismos bits). En la GPU Intel
    ifx implementa `exponent()` con una rutina que no existe con la FP64 emulada (D-019.5b) y `scale()` con
    operaciones FP64 emuladas, cuyas llamadas IGC compila mal en EU fusionadas (D-031); el lint
    `scale`/`exponent` de `check.sh` las rechaza;
  - ningún arreglo temporal de tamaño conocido sólo en ejecución: reducciones con `dim=` o expresiones de arrays
    sobre argumentos de forma asumida (`maxval(sum(abs(A), dim=2))`). En el dispositivo son un `alloca` dinámico
    que IGC no acota y el kernel desborda la memoria privada (D-030). Bucles explícitos; la norma infinito ya está
    en `NORM_INF`.
- **Kinds:** literales con kind (`1.0_${k}$`), nada de `real(8)`/`dble`/`1d0`; umbrales que dependen de la
  precisión, por kind en `muscle_kinds` (o como constante con nombre y documentación en el módulo si es propia de
  él). Ningún número mágico sin nombre.
- **Nombres:** el módulo de un modelo es `muscle_<prefijo>_<modelo>`, el de la familia `muscle_<prefijo>_laws`
  (D-027); tipo de parámetros `<prefijo>_<modelo>_<k>`; específicos `<prefijo>_<función>_<modelo>[_<dim>]_<k>`;
  `prm` primer argumento; `res` como resultado; `*_init` termina en `status`. Revisa la tabla por familia de
  `contributing.md` para el resto de los lugares (listas de `muscle_meta.py`, solvers, área de mutaciones).
- **Documentación FORD:** `!!` en la línea siguiente a cada declaración (módulo, tipo, procedimiento, argumento),
  fórmulas en `\( \)`/`\[ \]`, unidades y dominio de cada parámetro; referencia bibliográfica del modelo.
  Comentarios en inglés que explican el porqué, no repiten el código. Cita una decisión `D-0xx` sólo si explica
  algo no obvio.
- **Estilo:** 4 espacios, líneas ≤ 132 columnas (comentarios incluidos), cabecera SPDX, `use ..., only :`.
- **No toques** el orden de almacenamiento Voigt de ningún tipo ni el comportamiento numérico de lo existente sin
  decisión registrada; no corrijas de paso los problemas marcados en `docs/fypp/PROBLEMAS.md` como «Marcado, sin
  corregir»: el usuario decidió dejarlos así.

## 3. Verificar

En este orden (tiempos aproximados en `contributing.md`, *Development cycle*):

1. `cmake --build build --target build_all_tests` y `ctest --test-dir build -R <área>` mientras editas.
2. **Chequeo de mutaciones de lo nuevo**: `.claude/skills/muscle-nuevo-modulo/scripts/mutation_ids.sh <repo> <área> <archivo>` lista los ids de las
   mutaciones del archivo nuevo; córrelos con
   `/usr/bin/python3 tools/mutation_check.py <área> --samples all --only <ids>` (necesita Python ≥ 3.11: el
   Python del entorno `ford` no sirve). Una mutación que cambia un resultado y sobrevive es un test que falta;
   si es equivalente (no puede cambiar ningún resultado), explica por qué en el informe. Si el archivo no está en
   ningún área de `tests/mutations/*.toml`, agrégalo al área que corresponda.
3. **Ubicación y reutilización**: `.claude/skills/muscle-nuevo-modulo/scripts/check_new_module.py <repo> [archivos]` (por defecto, los cambiados en
   `src/` e `include/`). Señala: módulo con nombre distinto de su archivo o fuera de la carpeta de su familia;
   funciones fuera del contrato de la familia (las que no llevan el nombre del modelo son sospechosas de estar en el
   módulo equivocado); patrones de código que ya existe (traza o desviador por componentes, determinante o
   autovalores escritos a mano, umbrales sin nombre); bloques repetidos con el resto de la librería. Es una ayuda,
   no un veredicto: cada hallazgo se corrige o se justifica en el informe.
4. Si el cambio no debía alterar el código generado de lo existente (por ejemplo, sólo agregar un modelo a la
   lista de un módulo de familia), `tools/gen_snapshot.sh save` antes y `compare` después.
5. `tools/check.sh --quick`: lint (kinds, literales, GPU, 132 columnas, SPDX), build sin **ningún aviso del
   compilador**, ctest. Un aviso nuevo hace fallar el chequeo: arréglalo, no lo silencies.
6. `tools/check.sh` completo al cerrar (≈ 25 min: builds GPU, tests en la GPU Intel, auditoría). Para correr en
   GPU a mano: `systemd-run --user --scope -p MemoryMax=4G -p MemorySwapMax=0` + `timeout`; nunca
   `ulimit -s unlimited`.
7. Revisión independiente del diff (agente `reviewer`, si existe en tu copia) antes de pedir aprobación; pídele también que mire si
   algo reimplementa lo que ya existe o está en el módulo equivocado.

## 4. Cerrar

- Documentación de usuario: la lista de modelos del `README.md`, la tabla por familia de `contributing.md` si
  cambia algo de ella, `pages/glossary.md` si aparece un término o prefijo nuevo.
- `docs/fypp/`: decisión nueva en `DECISIONES.md` si la hubo; cualquier anomalía encontrada en `PROBLEMAS.md`.
- Informe breve al usuario: qué se creó (archivos), qué verdad comprueba cada test, resultado de ctest y del
  chequeo de mutaciones (detectadas / sobrevivientes con su explicación), los hallazgos de `check_new_module.py`
  (corregidos o justificados), avisos (deben ser 0), y las decisiones
  que tomó o que quedan abiertas. Commit local sólo si lo pide, en estilo Conventional Commits
  (`feat(hard): add the ... law`), un paso verificable por commit.
