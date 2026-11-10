# Estado óptimo para continuar el desarrollo de MUSCLE

Fecha: 2026-11-08. Rama: `fypp/f9-readability` (Fases 0–9 cerradas).
Fuentes: `ANALISIS_COMPRENSION.md`, `ANALISIS_REPETICION.md`, `ANALISIS_SOBREINGENIERIA.md`,
`ANALISIS_MANTENIBILIDAD.md` y las mediciones del ciclo de desarrollo de la sección 2.5.

**Objetivo:** dejar la librería lista para desarrollar de forma eficiente: que se entienda, que haya poco que
repetir, que no cargue con mecanismos que ya no pagan su costo y que cambiarla sea rápido y seguro.
**Fuera de alcance** (etapa posterior, por decisión del usuario): validación física, robustez numérica,
rendimiento, distribución, licencia y concurrencia.

Convención: VERIFICADO = medido en esta rama; ESTIMADO = juicio o cálculo propio. Esfuerzos en días-persona,
ESTIMADOS.

---

## 1. Resumen

El código de modelos quedó en buen estado tras la Fase 9: se lee como Fortran, un modelo vive en un archivo
y la red de seguridad (goldens dp bit a bit, 181 tests, snapshot del generado) es fuerte. Lo que frena el
desarrollo futuro no está en los modelos sino **alrededor**:

1. **Conocimiento frágil:** las reglas del proyecto no están en git y la rama vigente está muy lejos de
   `main` y `dev`.
2. **Ciclo de verificación lento o incompleto:** no hay modo rápido de `check.sh`, no hay CI, la cobertura no
   se puede medir y los goldens sólo valen en Debug.
3. **Guía de colaboración incompleta:** agregar un material exige pasos que no están escritos.
4. **Restos del diseño anterior:** el registro de funciones de familia en Python duplica firmas que ya están
   en Fortran.
5. **Repetición acotada:** dos focos claros en `src`; el grueso está en los tests heredados.

Corregir los bloques 0–3 de la sección 3 cuesta del orden de **8 a 11 días-persona** (ESTIMADO) y no toca la
numérica.

---

## 2. Diagnóstico por aspecto (cifras clave)

### 2.1 Comprensión
- Modelos legibles: 0 llamadas a macros que armen declaraciones (eran 207).
- Generación de materiales opaca: 4 plantillas, 1 147 líneas, 83 `#:if` de variante.
- Prueba del recién llegado (agregar una ley y un material siguiendo `contributing.md`): resultado parcial.
  Faltan en la guía `golden_params.py`, un material por solver y los goldens; los mensajes de error no
  nombran el archivo a editar.
- Sin glosario (`cp`, `cp2`, `pz`, `tm`, `cut`, `vm`, `hsv`, `SoA`) ni mapa de arquitectura.
- 384 referencias `D-0xx` dentro de `src/`.

### 2.2 Repetición
- `src`+`include`: 12 % de líneas en bloques repetidos (cota superior); `tests`: 22 %; `tools`: 0 %.
- Focos evitables en `src`: evaluación conjunta de yield (5 copias de 36 líneas) y preámbulo/cierre de los
  cuerpos de solver (~120 líneas).
- Información repetida: un modelo nuevo se declara en 6 lugares.

### 2.3 Sobreingeniería
- Único exceso claro: el registro de funciones de familia en Python (`MODEL_FAMILIES`, 222 líneas, y 5
  clases). Producción lee sólo prefijo, módulo, módulo de familia y modelos; el resto lo leen sólo los
  generadores de tests.
- Todo lo demás tiene consumidores: 0 funciones muertas en `muscle_meta.py`, las 12 macros de `include/`
  tienen 2–401 usos, los 9 puntos de entrada por material y los 16 métodos históricos de tensores se usan.
- El costo pesa más en el proceso: goldens obligatorios por material, `check.sh` sin modo intermedio.

### 2.4 Mantenibilidad
- `.claude/rules/` (7 reglas, 476 líneas) y `CLAUDE.md` están en `.gitignore`.
- `fypp/f9-readability` va 282 commits por delante de `main` y 146 de `dev`; `main` es el código anterior.
- Sin CI. Entorno sin fijar; el README dice «GFortran 9+» sin verificar.
- Goldens dependientes de flags: con `Profile_unittests` fallan `kin_golden` y `math_golden_derivatives`.
- Target de cobertura roto con lcov 2.x.
- Soluciones provisionales de compiladores (≈ 8) sin criterio de retiro.
- 52 avisos «Unused variable» en el build.

### 2.5 Ciclo de desarrollo (VERIFICADO, gfortran Debug, `-j8`)

| Operación | Costo |
|---|---|
| `ctest` completo (181 tests, CPU) | 9.6 s |
| Tocar una plantilla sin cambiar su salida | ~1.3 s (fypp regenera, nada recompila) |
| Agregar una constante pública a `muscle_kinds` | **reconstrucción total: 479 pasos, ~60 s** |
| Cambiar `include/muscle_common.fypp` o `muscle_meta.py` | recompila 142 de 144 objetos de la librería (y `muscle_meta.py` relanza CMake) |
| Cambiar la interfaz de una ley de endurecimiento | recompila los 76–80 módulos de material |
| Builds GPU de `check.sh` | 17 s (nvptx), 97 s (amdgcn), 175 s (spir64) más tests target en la Intel |
| `check.sh` completo | del orden de decenas de minutos (no cronometrado) |

Puntos de concentración: `muscle_tensors` (872 `use`), `muscle_yield_criteria` y `muscle_elasticity_laws`
(564 cada uno), `muscle_hard_laws` (536), `muscle_kinds` (439). `muscle_kinds` mezcla las definiciones de
kind (estables) con las tolerancias por kind (cambiaron en 13 de 103 commits de la migración).

---

## 3. Qué corregir, en bloques

Cada ítem: **qué**, por qué, esfuerzo, riesgo y criterio de terminado. «Decisión» marca lo que requiere al
usuario.

### Bloque 0. Proteger la base (≈ 1 día)

| ID | Qué | Por qué | Esfuerzo | Riesgo | Terminado cuando |
|---|---|---|---|---|---|
| 0.1 | ~~**Versionar las reglas del proyecto**~~ **No se hará (decisión del usuario, 2026-11-08).** (`.claude/rules/`, o una copia en `docs/dev/`), dejando fuera sólo los settings personales. **Decisión:** qué partes | Hoy un clon nuevo no tiene las reglas que aplican el lint y la revisión | 0.25 | Nulo | Las reglas están en git y `pages/`/`docs/` apuntan a una ruta versionada |
| 0.2 | ~~**Integrar la cadena de ramas**~~ **No se hará (decisión del usuario, 2026-11-08).** `f0…f9` en `dev` y decir en el README cuál es la rama vigente. **Decisión:** cuándo; el push lo hace el usuario | Lo que ve un externo es el código anterior; integrar más tarde cuesta más | 0.25 (más el push) | Bajo | `dev` contiene la Fase 9 |
| 0.3 | **Hecho (2026-11-08).** **Declarar versiones probadas** (gfortran 15.2, ifx 2026.2, fypp ≥ 3.1, CMake ≥ 3.20) y corregir «GFortran 9+»; guía breve del toolchain de offload | Una versión distinta rompe el build o los goldens sin aviso | 0.5 | Nulo | README con versiones verificadas |

### Bloque 1. Ciclo de desarrollo rápido y fiable (≈ 2.5 días; revisado con el usuario el 2026-11-08)

Orden de ejecución: 1.1 → 1.6 → 1.3 → 1.4. Un commit local por ítem; cada uno se verifica con `check.sh`
completo y, los que no deben cambiar la librería (1.1, 1.6), también con `tools/gen_snapshot.sh`.

| ID | Qué | Por qué | Esfuerzo | Riesgo | Terminado cuando |
|---|---|---|---|---|---|
| 1.1 | **Hecho (2026-11-08): 57 s** frente a decenas de minutos del completo. **`check.sh --quick`**: lint + build CPU `sp;dp` Debug + `ctest`; sin build UDTIO, sin FORD, sin GPU, sin auditoría. El completo queda para cerrar etapas | Hoy la elección es «sólo lint» o decenas de minutos (`ctest` solo: 9.6 s) | 0.5 | Muy bajo | `--quick` existe, está en `--help` y su tiempo queda medido |
| 1.6 | **Hecho (2026-11-08); corregido el 2026-11-10:** la medición de entonces fue sobre un build incremental (la librería no se recompiló); un build limpio tenía **514 avisos** (393 argumentos sin usar, 62 comparaciones exactas, el resto en tests). Bloque de avisos del re-análisis 4 (`1da8f79`, `e12d1c2`, `6f6cfba`): 0 avisos en build limpio con `-Wno-unused-dummy-argument` y, pendiente de retirar, `-Wno-compare-reals`; `check.sh` falla con cualquier aviso. Texto original: **Limpiar los 52 avisos**: todos están en `tests/materials/test_mat_kinematics.fypp` (variables sin usar); la librería no tiene avisos | Tapan avisos nuevos | < 0.5 | Nulo | Build Debug sin avisos; el test comprueba lo mismo |
| 1.3 | **Hecho (2026-11-08, D-021):** 186/186 en Debug, Release y `Profile_unittests`; modo estricto pasa bit a bit en Debug; con FMA quedan 6 problemas de diseño de tests en PROBLEMAS. **Goldens con tolerancia por categoría** (decisión del usuario: los goldens eran para la migración; ahora son tests de regresión de funcionalidad, no de bits). Ver detalle abajo | Fallan con otros flags/compiladores (2 falsos fallos con `Profile_unittests`) y bloquean optimizaciones legítimas | ~1 | Bajo | Pasan en Debug, `Profile_unittests` y Release; un cambio deliberado mayor que la tolerancia hace fallar el test; modo estricto disponible |
| 1.4 | **Hecho (2026-11-08):** librería (`src/`+`include/`) 98.2 % líneas, **81.7 % funciones** (3 750 de 4 588 procedimientos generados), 91.9 % ramas. Las líneas están infladas: gcov atribuye el código generado a la línea de la plantilla (una línea de un `#:for` cuenta cubierta si se ejecuta cualquier instancia); el dato representativo es el de funciones. **Arreglar la cobertura** (lcov 2.0 falla con los `#line`: `--ignore-errors inconsistent`, `--rc geninfo_unexecuted_blocks=1`, o gcovr) y medirla una vez. Después de 1.3 | No se sabe qué cubren los tests | 0.5 | Nulo | Informe generado; cifra anotada aquí y en `PROBLEMAS.md` |
| 1.2 | ~~Separar `muscle_kinds`~~ **No por ahora**: 89 de los 133 módulos que lo usan traen alguna constante (incluidos los 80 materiales); la cascada bajaría de 133 a ~89 | Beneficio menor que el costo de tocar 133 `use` | — | — | — |
| 1.5 | ~~CI mínima~~ **Pospuesta**: sólo se activa con el push de una rama que la contenga, y las ramas no se integran por ahora (0.2) | — | — | — | Se retoma al decidir la integración |

**Detalle de 1.3 (goldens con tolerancia):**
- **Regla nueva por defecto:** cada golden compara con tolerancia relativa por categoría, en todas las
  configuraciones. Valores de partida (dp; en sp se escalan con su epsilon y su `NEWTON_TOL`), a calibrar
  midiendo las diferencias reales entre Debug, `Profile_unittests` y Release:

  | Categoría | Tolerancia relativa de partida | Por qué |
  |---|---|---|
  | Álgebra directa (tensores, leyes, cinemática) | ~1e-13 | Sólo redondeo |
  | Diferencias finitas | ~1e-9 | El redondeo se amplifica por ε/h |
  | Solvers iterativos (tensión, estado, tangente) | ~10 × `NEWTON_TOL` (≈ 1e-9) | Otro redondeo puede cambiar la trayectoria del Newton |
  | Número de iteraciones | ±1 | Idem |
  | `status` | exacto | Es lógica, no aritmética |

- **Modo estricto opcional:** una variable de entorno reactiva la comparación bit a bit en la configuración de
  referencia (gfortran Debug), para cuando se quiera certeza total.
- **Refactors que no deben cambiar nada:** se verifican con `tools/gen_snapshot.sh` (código generado igual),
  que es mejor herramienta para ese caso que el bit a bit.
- **Se toca:** la comparación de `tests/materials/mat_golden_common.fypp` y las de cinemática, math, visco y vp;
  constantes de tolerancia por categoría en `tests/common/`; las referencias **no** se regeneran.
- **Se documenta:** decisión nueva en `DECISIONES.md` (D-021); la regla de tests (`.claude/rules/testing.md`) y
  la regla dura de `CLAUDE.md` «no cambiar el comportamiento numérico en dp», que pasan de «mismos bits» a
  «dentro de la tolerancia de su categoría».
- **Fuera de alcance:** ampliar tests de propiedades (validación física, etapa posterior).
- **Nombre:** se mantiene «golden» por ahora (decisión del usuario, 2026-11-08). Se corrigen los 81 mensajes
  «bit-for-bit» de los tests («outside the tolerance of its category (rel. X)»; «bit-for-bit» sólo en modo
  estricto) y se define el término en la regla de tests: *golden = valores de referencia guardados en
  `tests/reference/` desde una versión validada, comparados con la tolerancia de su categoría*. El renombre
  queda para 3.8.

### Bloque 2. Que se entienda sin preguntar (≈ 2–3 días)

| ID | Qué | Por qué | Esfuerzo | Riesgo | Terminado cuando |
|---|---|---|---|---|---|
| 2.1 | **Hecho (2026-11-08):** listas de verificación y página pública `pages/conventions.md` en FORD (el usuario descartó versionar `.claude/`). **`contributing.md` completo**: listas de verificación para modelo, material y solver (CMake, `golden_params.py`, un material por solver, goldens, `write_golden.sh`, tests), el `use` de la familia en su propia línea, el flujo de build (Ninja, `build_all_tests`, `ctest -R`, `--quick`); sin remisiones a archivos no versionados | La prueba del recién llegado falló justo ahí | 0.5 | Nulo | Repetir la prueba del recién llegado sin pasos no documentados |
| 2.2 | **Hecho (2026-11-08).** **Mensajes que dicen qué archivo editar** (`tests/materials/props_samples.py:60`, `tests/tools/test_muscle_materials.py:47`, `mat_golden_common.fypp:138`) | Hoy nombran el síntoma | 0.25 | Nulo | Cada error de registro incompleto nombra el archivo y la entrada |
| 2.3 | **Hecho (2026-11-08).** **Glosario** en `pages/` | Nombres crípticos del registro | 0.25 | Nulo | Página publicada en FORD |
| 2.4 | **Hecho (2026-11-08), con extracto (decisión del usuario: sin archivo generado versionado).** **`pages/architecture.md`**: diagrama registro → plantilla de material → solver → modelos → API de punto y de lote, con el `.f90` generado de un material de ejemplo, regenerado y comparado en `check.sh` | Falta el mapa | 0.5–1 | Nulo | Página publicada; el ejemplo no puede quedar viejo |
| 2.5 | **Hecho (2026-11-08).** **Índice de secciones en `muscle_meta.py`** (qué consume cada bloque) | Sin guía de navegación en 2 287 líneas | 0.25 | Nulo | Índice al inicio del archivo |
| 2.6 | **Hecho (2026-11-08):** 14 filas en `PROBLEMAS.md`. **Registro de soluciones provisionales**: compilador, versión, archivo, test que la verifica y condición de retiro | Se acumulan sin fecha de salida | 0.25 | Nulo | Tabla en `PROBLEMAS.md` o en las reglas GPU |
| 2.7 | **Hecho (2026-11-08): 384 → 77**, snapshot idéntico; cabeceras de cp2 y cutting plane completas. **Recortar las referencias `D-0xx`** en `src/` a las que explican un porqué no obvio; completar las cabeceras de `closest_point_2_body` y `cutting_plane_body` | Ruido para quien no tiene `DECISIONES.md` | 0.5 | Nulo | ≤ ~100 referencias; cabeceras completas |
| 2.8 | **Hecho (2026-11-08):** fortls 3.2.2 (entorno `ford`) funciona sobre el código generado, no sobre las plantillas; `.fortlsrc` documentado en `contributing.md` y `tools/generated.sh` da el `.f90` de una plantilla. **Soporte de editor para `.fypp`**: investigar y documentar cómo navegar las plantillas (fortls sobre el generado, `compile_commands.json`, comando para abrir el `.f90` de una plantilla) | Sin «ir a la definición» en las plantillas | 0.5 | Nulo | Guía de configuración probada |

### Bloque 3. Menos repetición y menos andamiaje (≈ 3–4 días; cada ítem se verifica con el snapshot)

| ID | Qué | Por qué | Esfuerzo | Riesgo | Terminado cuando |
|---|---|---|---|---|---|
| 3.1 | **Hecho (2026-11-08):** `muscle_meta.py` −226 líneas (sin firmas, `ModelArg` ni campos sin uso); las firmas que necesita el test generado viven en una tabla de 18 líneas en `tests/gpu/target_cases.py`, con un test que la compara con el registro; Fortran generado de los tests idéntico. **Reducir el registro de funciones de familia** a lo que produce nombres (prefijo, módulo, modelos) y reescribir `tests/gpu/target_cases.py` y `test_family_target.fypp` con una tabla pequeña o a mano | Duplica en Python firmas que ya están en Fortran (~250–300 líneas) | 1–1.5 | Bajo (sólo datos y tests; el generado de `src` no cambia) | `MODEL_FAMILIES` sin `args`, `funcs`, `joint` ni `composition`; tests target iguales |
| 3.2 | **No se hará (decisión del usuario, 2026-11-08): repetición aceptada para que cada criterio se lea completo.** **Evaluación conjunta de yield** como texto fijo con sustitución de nombres | 5 copias de 36 líneas | 0.5 | Bajo | Una sola copia; snapshot idéntico |
| 3.3 | **No se hará (2026-11-08): la repetición es sobre todo la firma de `SOLVE`, escrita a mano en cada cuerpo por D-020.** **Preámbulo y cierre de los cuerpos de solver** como texto fijo compartido | ~120 líneas repetidas | 0.5 | Bajo | Snapshot idéntico |
| 3.4 | **No se hará (2026-11-08): la prueba del recién llegado regeneró los goldens sin problema con la guía.** **Golden de un material nuevo**: script que lo genere, o escape explícito y listado. **Decisión:** cuál | Hoy son 3 pasos manuales | 0.5 | Bajo | Agregar un material es un comando más |
| 3.5 | **No se hará (2026-11-08): la guía corregida basta; el recién llegado agregó su test sin trabarse.** **Un test pequeño por ley** como patrón copiable, en lugar de ampliar los monolitos | La prueba del recién llegado no supo dónde ponerlo | 0.5 | Nulo | Patrón documentado en `contributing.md` |
| 3.6 | **Hecho (internas, 2026-11-08):** `muscle_kin_finite` (929 líneas) seguía generando 24 funciones desde tablas de Python (`MEASURES`, `CONV`, ...), una violación de D-020 que las métricas no veían (no usaba macros `MODEL_*`): reescrito a mano, snapshot idéntico salvo orden; párrafo «cómo agregar un modelo» en los 7 módulos de familia. `*_init` sin `status`: el usuario decidió agregarlo a todos (D-022, hecho). **Consistencia de convenciones**: inventario de desviaciones entre familias (sufijo `_3D_` en yield sí y en hard no, nombres de `*_init`, estilo de tests) y corrección de las internas. **Decisión:** las que cambian API pública | Menos carga mental | 0.5–1 | Bajo (internas) / requiere decisión (públicas) | Inventario hecho; internas corregidas |
| 3.7 | **Decidido (usuario, 2026-11-08): se mantiene la escritura explícita.** **Lista de modelos en los módulos de familia**: explicitud (hoy) o lista única con `#:for`. **Decisión** | Información repetida en 3 lugares por familia | 0.25 | Bajo | Decisión registrada en D-020 |
| 3.8 | **Reemplazado por el bloque R (2026-11-09).** Antes: **Más adelante (decisión del usuario, 2026-11-08): renombrar «golden» → «reference» y reordenar los tests.** Hoy: 28 archivos con «golden» en el nombre, 933 menciones en `tests/`/`src/`/`tools/`, 22 tests de ctest; tests `.f90` heredados y plantillas mezclados, monolitos de hasta 1 057 líneas. Un commit sólo de renombre y otro de reorganización, separados del cambio de tolerancias (1.3) | Nombre que describe lo que hace; estructura de tests navegable | 1–2 | Bajo (los nombres de ctest cambian: avisar en `CLAUDE.md`) | Sin «golden» en nombres; estructura de `tests/` documentada |

### Bloque R. Tests de referencia → tests unitarios (decidido con el usuario, 2026-11-09)

Reemplaza y amplía el ítem 3.8. **Objetivo:** reducir al mínimo los tests de referencia (hoy «golden»: 22 tests,
~8 600 registros en `tests/reference/`) y convertirlos en tests unitarios con una verdad independiente (forma
cerrada, propiedades, otra ruta de cálculo). **Decisiones del usuario:** (1) se acepta un mínimo de referencias donde
no hay verdad independiente (solvers sobre secuencias de carga, criterios anisótropos); (2) **criterio de
mutaciones**: una referencia sólo se borra si el test unitario que la reemplaza detecta las mismas mutaciones
deliberadas del código (cambiar un coeficiente, intercambiar componentes, alterar una fórmula); (3) se hace todo el
plan. Estructura de `tests/`: se mantiene el espejo de `src/`, con etiquetas de ctest por tipo (`kinds`, `legacy`,
`reference`, `target`); los `.f90` heredados no se mueven (etiqueta `legacy` + nota en la guía); los monolitos no se
parten por ahora; se renombra todo «golden» → «reference» (archivos, API `reference_init/reference_check`,
`REFERENCE_TOL_*`, `MUSCLE_REFERENCE_STRICT`, ctest), con nota histórica en el glosario.

| Etapa | Qué | Esfuerzo |
|---|---|---|
| R0 | Inventario: por test de referencia, qué comprueba, qué tests unitarios ya lo cubren, qué falta y qué referencia queda; un juego de mutaciones por área | 0.5 |
| R1 | Operadores de tensores: contra el cálculo con matrices completas (ampliar `ops_full_common.fypp`); borrar las 6 referencias | 0.5–1 |
| R2 | Math: sistemas lineales (residuo, A·A⁻¹ = I), operaciones (invariantes de autovalores), derivadas (polinomios, analítica vs diferencias), espectral (V Λ Vᵀ = A, ortogonalidad) | 1–1.5 |
| R3 | Elasticidad (forma cerrada, embebido 2D/3D), endurecimiento y leyes vp (fórmulas cerradas, derivadas vs diferencias), viscoelasticidad (relajación de Maxwell, fluencia de Kelvin) | 1.5–2 |
| R4 | Cinemáticas: estiramiento uniaxial, corte simple, rotación rígida en forma cerrada; objetividad | 1 |
| R5 | Fluencia: von Mises y Hill cerrados; homogeneidad, gradiente/hessiano vs diferencias, reducción isótropa; referencias mínimas para CPB06/Yld2004 anisótropos | 1–1.5 |
| R6 | Solvers: retorno radial cerrado para von Mises lineal; consistencia f ≈ 0, paso elástico, tangente vs diferencias; referencias mínimas por solver | 1.5–2 |
| R7 | Renombre «golden» → «reference» de lo que quede, arreglos de estructura (mover `tests/test_muscle_math_operations.f90` a `tests/math/`, `muscle_test_solver_golden.f90` a `tests/materials/`, unir `reference/solvers` y `reference/materials`, documentar `tests/gpu/`), etiquetas de ctest, guía | 1 |

Un commit por etapa (tests nuevos + borrado de referencias juntos, tras pasar el criterio de mutaciones);
`check.sh --quick` en cada uno y `check.sh` completo al final. Meta: ~6–8 tests de referencia con unos pocos cientos
de registros.

**Resultado (2026-11-09, hecho, D-023).** Chequeo con **todas** las mutaciones de cada área (`tools/mutation_check.py
<area> --samples all`), no con una muestra:

| Etapa | Área | Mutaciones que detectaba la referencia | Unitarios antes → después | Qué quedó |
|---|---|---|---|---|
| R0 | herramienta `tools/mutation_check.py`, `tests/mutations/*.toml` | — | — | — |
| R1 | operadores de tensores (+ 8 a mano en `muscle_meta.py`) | 25 | 25 → 25 | nada (6 tests borrados; `tools-ops_inventory` reemplaza la cobertura) |
| R2 | math | 655 | 633 → 649 | `math_reference`: 6 registros (empates exactos; rama de `eigenvals`) |
| R3 | hard / elas / vp / visco | 48 / 152 / 322 / 194 | 48 / 152 / 304 / 187 → todas | nada (test nuevo `vp_limits`) |
| R4 | kin / mat_kinematics | 245 / 237 | 166 → 245 / 236 → 237 | nada (test nuevo `kin_finite`, verdad en precisión cuádruple) |
| R5 | yield | 1 286 | 1 281 → 1 286 | nada |
| R6 | solvers | 130 | 102 → 117 (+ referencia reducida: 130) | `closest_point_time` 50 y `closest_point_perzyna` 27 registros: camino de Newton |
| R7 | renombre y orden | — | — | etiquetas `kinds/legacy/reference/target`, guía (Tests), D-023 |

De 22 tests y ~8 600 registros quedan **3 tests y 83 registros**. Desvíos del plan: `tests/test_muscle_math_operations.f90`
no se movió sino que se borró (no estaba registrado, no compilaba y lo cubre `tensors-types-3d2osym`);
`reference/materials` desapareció con R4, así que no hubo nada que unir; `muscle_test_solver_golden.f90` pasó a
`tests/materials/muscle_test_solver_sequence.f90`. Hallazgos (PROBLEMAS): `check_approx` aceptaba Inf (corregido),
el mismo hueco en `.approx.` (corregido a pedido del usuario), y marcados sin corregir: la rama `|1 - d| <= EPS` de
`eigenvals`, el jacobiano de `closest_point_2` sin factor 2 de corte y `closest_point_time` + JC sin convergencia con
poca sobretensión.

### Bloque 4. No ahora (sólo si aparece la necesidad)
- **Alternativa D** (tipos vacíos en materiales; parche en `build-snapshots/d_visco_kin.patch`), **B** o **C**:
  beneficio modesto frente al costo; D además reintroduce generación desde metadatos.
- **Tests `.f90` heredados a tablas de casos** (49 archivos, mucha repetición interna): son la red de seguridad
  de la migración; esperar a tener cobertura medida (1.4).
- **Cuerpos 3D/2D escritos a mano** y **andamio de firmas de operadores**: decisiones de estilo de la Fase 9;
  se mantienen.
- **Carga sin uso interno** (2D en materiales, UDTIO, esqueletos FE): decidir en la etapa posterior.

---

## 4. Orden recomendado

1. **Bloque 0** completo (base y conocimiento).
2. **1.1, 1.2, 1.6** (ciclo rápido) y **2.1, 2.2** (lo que falló en la prueba del recién llegado).
3. **1.3, 1.4, 1.5** (red de seguridad portable).
4. **3.1, 3.2, 3.3, 3.4, 3.5** (refactors acotados).
5. **2.3–2.8, 3.6, 3.7** (documentación y consistencia).
6. **Repetir la prueba del recién llegado** y las métricas (sección 5) para confirmar el estado.

## 5. Cómo saber que se llegó

| Criterio | Hoy (VERIFICADO) | Meta |
|---|---|---|
| Reglas del proyecto en git | no | descartado (0.1) |
| Rama vigente integrada en `dev` | 146 commits por delante | descartado por ahora (0.2) |
| CI | no | pospuesta (1.5) |
| Prueba del recién llegado | parcial, 3 pasos no documentados | **lograda (2026-11-08)**: ley + 3 materiales + 45 registros golden siguiendo sólo la guía, 186/186; fricciones menores corregidas en la guía |
| Verificación diaria | sólo lint o `check.sh` completo | `--quick` en pocos minutos (medir al implementarlo) |
| Constante de tolerancia nueva | reconstrucción total (~60 s) | sin meta por ahora (1.2 pospuesto) |
| Avisos del build Debug | 52 (medido mal: 514 en build limpio) | **0 en build limpio** (2026-11-10), comprobado por `check.sh`; `-Wno-compare-reals` pendiente de retirar |
| Cobertura | no medible | medida: 81.7 % de funciones de la librería (98.2 % líneas, inflado por las plantillas) |
| Goldens fuera de Debug | 2 falsos fallos (comparación bit a bit) | 0, con tolerancia por categoría y modo estricto opcional |
| Tests de referencia | 22 tests, ~8 600 registros | **3 tests, 83 registros** (plan R, criterio de mutaciones) |
| Tests heredados `.f90` (`legacy`) | 38 tests, 9 650 líneas | **0** (D-025, 2026-11-10: los unitarios detectan todas sus mutaciones; re-análisis 6) |
| `muscle_meta.py` | 2 287 líneas | **2 093 tras 3.1** (sin el registro de firmas) |
| Bloques repetidos en `src` (`tools/dup_report.py`) | 12 % | ≤ 10 % |
| Macros prohibidas (`tools/fypp_metrics.sh`) | 0 | 0 |

## 6. Decisiones pendientes del usuario
1. Qué partes de `.claude/` se versionan (0.1).
2. Cuándo se integran las ramas en `dev` (0.2).
3. Tolerancias de los goldens fuera de la configuración de referencia (1.3).
4. Dónde corre la CI y con qué compilador (1.5).
5. Script o escape para el golden de un material nuevo (3.4).
6. Cambios de convención que afectan la API pública (3.6).
7. Explicitud o lista única en los módulos de familia (3.7).
