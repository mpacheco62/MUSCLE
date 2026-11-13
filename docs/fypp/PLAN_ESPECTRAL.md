# Plan C: una sola descomposición espectral (2026-11-11)

Pedido del usuario tras `REANALISIS_MATH.md`: hacer la opción C, «unificar los algoritmos espectrales», con un
plan que se pueda desarrollar sin mayores inconvenientes. Este documento es ese plan: qué significa unificar, qué
se mide antes de elegir el algoritmo, en qué pasos se hace, cómo se verifica cada uno y qué decisiones quedan
para el usuario. No se cambió código.

## 0. Criterio del usuario (2026-11-11) y orden resultante

> «Lo ideal es que el kernel no aumente mucho más, ya que en un futuro optimizaremos estos casos para correrlos en
> GPU. Pero un punto importante es que quede ordenada esa carpeta y se entienda a primera vista qué hace qué. En
> caso que sea muy complejo unificar alguna o aumente mucho el kernel, podemos dejarla separada, pero ordenada.»

Consecuencias para el plan:

1. **Prioridad 1: el orden de la carpeta** (opción B de `REANALISIS_MATH.md`), con resultados idénticos bit a bit.
   Se hace **primero** (fase B, abajo): es de bajo riesgo, deja la carpeta legible en seguida y las fases de C
   trabajan ya sobre la estructura final.
2. **Prioridad 2: unificar sólo lo que salga barato.** Presupuesto R7: el lote GPU y la llamada en CPU de un
   material anisótropo y de uno `kin_ml` no crecen más de ~1.2x (valor de partida; el usuario puede ajustarlo
   con los números de F0). El candidato J (Jacobi para todo) probablemente no lo cumple en los criterios
   anisótropos; K o KJ sí podrían. **Si un consumidor no cabe en el presupuesto, conserva su algoritmo**, en el
   mismo módulo ordenado y con la diferencia documentada en una línea («por qué hay dos»).
3. Las fases F2 (cinemática) y F3 (criterios anisótropos) son independientes: cada una se unifica o no según su
   propia medida.

### Fase B: ordenar `src/math/` (antes de F0) — **hecha** (D-029)

Distribución final (cada módulo con una cabecera «qué hay aquí y cuándo usarlo»):

| Módulo | Contenido | Capa |
|---|---|---|
| `muscle_math_inverses` (include `muscle_math_lu_bodies.fypp`, renombre de `muscle_math_lu_inverse.fypp`) | `M33INV` … `M99INV`, `COND_INF`, `ADJ_KAPPA_MAX` | bajo tensores |
| `muscle_math_lu_solve` | `LU<n>SOLVE`, `lu_solve` | bajo tensores |
| `muscle_math_fp` (nuevo) | `QUIET_NAN`, `fp_exponent`: reemplazos de IEEE/intrínsecos que no existen en el dispositivo | bajo tensores |
| `muscle_math_eigen` (nuevo) | `eigenvals`, `spectral_decomposition`, `eigenprojections`, `eigenvals_refined` (Jacobi privado); luego el núcleo de C | sobre tensores |
| `muscle_math_spectral_derivs` | `dEigenvalues_dTensor`, `d2Eigenvalues_dTensor2`, `spectral_gradient`, `spectral_hessian`, `dsquare` (absorbe `muscle_math_derivatives`) | sobre tensores |
| `muscle_math_tensor_functions` (renombre de `muscle_math_operations`) | `polar_rotation`, `log_divided_difference` | sobre tensores |
| `muscle_math_linsolvers` | host, sin cambios por ahora (retirar o no: decisión aparte) | host |

Más una página `pages/math.md` con la tabla «necesito X → módulo Y» y el dibujo de las capas, enlazada desde
`architecture.md`. Los tests se renombran para espejar los módulos (D-025) y las áreas de mutaciones se
actualizan. Sin alias (D-027). Verificación: `gen_snapshot` (sólo cambian nombres de módulos y símbolos), ctest,
referencias estrictas, `check.sh` completo. Commits: uno por movimiento (fp, eigen, spectral_derivs +
derivatives, tensor_functions, include + estilo, tests + página).

## 1. Qué significa «unificar»

Hoy hay **cinco caminos** para el mismo problema (autovalores y autovectores de un tensor simétrico 3x3), con
**dos estrategias distintas** para los pares de autovalores cercanos:

| Rutina | Módulo | Método | Pares cercanos | Consumidores |
|---|---|---|---|---|
| `eigenvals` | operations | forma cerrada trigonométrica | umbral absoluto 1e-10 sobre el desviador (triple raíz), rama `abs(1-d) <= EPS` | f de CPB06, CPB06ex2, Yld2004; dentro de las otras cuatro |
| `spectral_decomposition` | spectral_derivs | `eigenvals` + producto vectorial de filas (Kopp) + **una** rotación plana | rotación plana del par restante | gradiente/hessiana de CPB06, CPB06ex2, Yld2004 |
| `eigenprojections` | spectral_derivs | `eigenvals` + Cayley-Hamilton fuera de la banda, **Jacobi cíclico** dentro | banda `JACOBI_GAP_LO/HI` (dp [1e-12, 1e-2); sp siempre Jacobi) | cinemática finita (log, polar en dp) |
| `eigenvals_refined` | spectral_derivs | lo mismo, sin proyecciones | ídem | cinemática (`update_norot`) |
| `dEigenvalues_dTensor` | spectral_derivs | Cayley-Hamilton, reparto isótropo en raíces repetidas | umbral relativo `EPS_TOL` | `eigenprojections` (fuera de la banda) |

**Objetivo de C:** un único núcleo `(λ, V) = descomposición(S)`, robusto y preciso para cualquier separación entre
autovalores, del que sale todo lo demás:

- autovalores precisos: `λ`;
- autoproyecciones: `Pₐ = vₐ ⊗ vₐ`;
- gradiente y hessiana espectrales: `spectral_gradient(V, …)` y `spectral_hessian(λ, V, …)`, que no cambian;
- medidas logarítmicas y sus tangentes, desde `(λ, P)`, como hoy.

Desaparecen la banda de Jacobi, la doble estrategia y la dependencia de las rutinas vectoriales de los umbrales
absolutos de `eigenvals`.

**Fuera del núcleo, con decisión propia (fases F4 y F5):**
- `eigenvals` como camino barato «sólo valores» para la f de los criterios anisótropos;
- `dEigenvalues_dTensor` y `d2Eigenvalues_dTensor2` como API Cayley-Hamilton (sin vectores).

## 2. Qué toca y qué no (impacto medido)

| Elemento | Impacto | Fuente |
|---|---|---|
| Referencias de solvers (77 registros, `closest_point_time`/`perzyna`) | **ninguno**: sólo von Mises y Drucker-Prager | `tests/reference/solvers/*.txt` |
| `math_reference` (6 registros) | cambian los de los empates exactos de `spectral_decomposition` (qué base de un autoespacio repetido elige); los de `eigenvals` sólo si F4 cambia `eigenvals` | cabecera de `math_reference_dp.txt`, D-023 |
| Materiales `*_cpb06*`, `*_yld2004*` (cp, cp2, pz, tm, cut) y `kin_*` finitos | resultados cambian al nivel del redondeo o de eps/gap; los tests son por propiedad o con tolerancia; **cambian iteraciones de Newton** posiblemente en ±1 | `tests/materials/*`, `muscle_materials.py` |
| Tests de yield (`yield_laws_kinds`: gradiente y hessiana contra diferencias finitas) | deben seguir pasando; la f no cambia mientras no se toque `eigenvals` | `tests/plasticity/yield_criteria` |
| Tests de cinemática (verdad en cuádruple) | deben seguir pasando; es la red principal | `tests/kinematics` |
| GPU Intel (IGC) | **riesgo alto**: compilador frágil con estas rutinas (`scale()` D-019.6; `closest_point` con más de un work item D-019.5c) | PROBLEMAS |
| Problemas «Marcado, sin corregir» | sólo `eigenvals` (rama `abs(1-d)`); no se toca salvo que el usuario lo autorice en F4 | PROBLEMAS.md |
| API pública | `spectral_decomposition`, `eigenprojections`, `eigenvals_refined`, `spectral_gradient`, `spectral_hessian`: **se conservan las firmas**; cambia la implementación. `JACOBI_GAP_*` de `muscle_kinds` se retiran | — |

## 3. Requisitos del núcleo (criterios de aceptación medibles)

Con `S = Q diag(λ) Qᵀ` conocido en precisión cuádruple, `ε = epsilon(1.0_k)` y `c` una constante chica fijada en
F0:

| Id | Requisito | Medida |
|---|---|---|
| R1 | autovalores con error absoluto ≤ c·ε·‖S‖ para **toda** separación (0, 1e-16 … 1) | máximo sobre la batería |
| R2 | `‖VᵀV − I‖ ≤ c·ε` y `‖V Λ Vᵀ − S‖ ≤ c·ε·‖S‖` | ídem |
| R3 | cada proyección con error ≤ c·ε/gapₐ (relativo); la suma de las de un par repetido, ≤ c·ε | ídem |
| R4 | ninguna excepción de punto flotante bajo `-ffpe-trap=zero,overflow,invalid`: tensor nulo, isótropo, subnormal, escalas 1e±30 en dp y 1e±18 en sp, pares a 1 ulp | test con trampas |
| R5 | continuidad: sin saltos de bandas o ramas mayores que c·ε | barrido de gap 1e-16 → 1 |
| R6 | kernel: `pure`, `@:DEVICE()`, sin llamadas de runtime (`audit_kernels.sh`), pasa los tests target en la GPU Intel | `check.sh` completo |
| R7 | costo ≤ presupuesto del usuario frente a hoy, en CPU por llamada y en el lote GPU de un material anisótropo y de uno `kin_ml` | banco de F0 y benchmarks |

R1–R3 son los que hoy cumplen `eigenprojections` (cinemática) y, salvo los umbrales de `eigenvals`,
`spectral_decomposition`. Ningún algoritmo cumple R1 con error *relativo* por autovalor: para C ≈ I la parte
desviadora de 1e-11 se conoce a ~ε/1e-11. Eso es inherente, también a Jacobi.

## 4. Candidatos (F0 los mide antes de elegir)

| Id | Algoritmo | A favor | En contra |
|---|---|---|---|
| **J** | Jacobi cíclico (el `jacobi_eigen` actual, ya probado en cuádruple y con D-028), **optimizado**: rotaciones con actualización explícita de 2 filas/columnas en lugar de dos `matmul` 3x3 por rotación (~10x menos operaciones), barridos fijos (sin divergencia en GPU) o parada por umbral | el más preciso y robusto; sin umbrales absolutos; ya es el camino de sp en la cinemática | el más caro: 4–6 barridos de 3 rotaciones |
| **K** | Kopp mejorado: `spectral_decomposition` actual (vector más separado por producto vectorial + rotación plana del par) con autovalores recalculados como cocientes de Rayleigh `vₐᵀ S vₐ` en lugar de `eigenvals`, y la rotación plana sin desborde (H3 de REANALISIS_MATH) | el más barato, sin iteración | depende de `eigenvals` para elegir el vector más separado; precisión a medir cerca de C = I |
| **KJ** | K seguido de un barrido de Jacobi de limpieza | casi el costo de K y la precisión de J | dos algoritmos en uno (más código) |

F0 entrega una tabla: precisión (R1–R5) por candidato, kind, separación y escala; costo en CPU por llamada; y un
prototipo en la GPU Intel. La elección es la **decisión D-032**.

## 5. Fases

> **D-032 (2026-11-12):** K, implementado como `spectral_decomposition` (sin `sym_eigen`). Como el núcleo conserva
> el nombre, F1 y F3 se hacen en un paso junto con F4 (b) (la f de los criterios usa los λ del núcleo); después F2
> (cinemática), y F5 mantiene las derivadas de Cayley-Hamilton. Donde abajo dice `sym_eigen`, léase
> `spectral_decomposition`.
>
> **F1 + F3 + F4 (b) hechos (2026-11-12):** el cuerpo de `spectral_decomposition` es K y la f de CPB06, CPB06ex2 e
> Yld2004 usa sus λ (CPB06 construye Σ con `yield_lin_apply`, como el gradiente). Los registros de empates de
> `math_reference` se conservan: K elige la misma base que la rutina anterior y siguen pasando (F3.3 no hace
> falta). Tests nuevos: corte diminuto sobre una presión, par casi repetido arriba y abajo, acoplamiento casi
> `tiny` (H3), media grande y propagación de NaN. Iteraciones de Newton idénticas en los 4 materiales
> anisótropos; GPU Intel 48/48.
>
> **F2 hecho (2026-11-12):** `eigenprojections`/`eigenvals_refined` sobre `spectral_decomposition`; sin franja de
> Jacobi (D-018.19a reemplazada). Ver D-032 para el cambio numérico de la deformación de Hencky. Para F6: el nombre
> `eigenvals_refined` ya no describe la rutina (no refina nada): candidata a retiro o alias si se limpia la API.

Cada fase es uno o más commits verificables, con todos los tests en verde. Ninguna empieza sin cerrar la
anterior.

### F0. Banco de pruebas y elección (sin tocar `src`) — **hecha** (`ESPECTRAL_F0.md`, 2026-11-12)

1. **Arnés de verdad en cuádruple** reutilizable, en `tests/common/muscle_test_spectral.fypp`. Generador de
   `S = Q diag(λ) Qᵀ`:
   - familias de separación: distintos, par arriba, par abajo, par al medio, triple, a 1 ulp, 1e-16 … 1e-1;
   - escalas: tensiones 1e-3 … 2e5 y C = I + δ con δ = 1e-11 … 1;
   - rotaciones: 2 o más, incluida la identidad (tensor diagonal: empates exactos);
   - medidas de R1–R5.
   Sale de generalizar la verdad de `test_muscle_math_spectral_derivs_projections`.
2. **Prototipos J, K y KJ** dentro del banco (`benchmarks/math/bench_spectral.fypp`, nuevo), no en `src`.
3. **Microbanco de costo en CPU:** ns por llamada en sp y dp, Release, frente a `spectral_decomposition` y
   `eigenprojections` de hoy.
4. **Humo en GPU Intel:** el prototipo elegido dentro de un `target`, con límites de memoria y tiempo, para
   detectar pronto si IGC lo rechaza.
5. **Informe** `docs/fypp/ESPECTRAL_F0.md` con la tabla y la recomendación. Luego la **decisión D-032** del
   usuario: algoritmo, presupuesto de costo (R7) y constantes c de los criterios.

**Commit:** `test(math): quad-precision spectral harness and candidate benchmark`.

### F1. El núcleo, sin consumidores

1. Nuevo genérico público con el algoritmo elegido: `sym_eigen(s, lam, V)` (nombre a confirmar), en `sp` y `dp`.
   Va en un módulo nuevo `muscle_math_eigen`, que es el primer paso de la opción B; la decisión se toma con
   D-032. Si se prefiere no crear el módulo todavía, va en `muscle_math_spectral_derivs`.
2. Tests unitarios con el arnés de F0, que cubren R1–R5 bajo `-ffpe-trap`. Se suma a
   `test_muscle_math_target.fypp`: CPU = GPU bit a bit en el host, `check_target` en el dispositivo.
3. Chequeo de mutaciones del área `math` sobre el núcleo (`--only`): toda sobreviviente se explica o se cubre.
4. `gen_snapshot`: sólo se agregan símbolos; nada existente cambia.
5. `check.sh` completo, porque es la primera vez que el núcleo pasa por IGC.

**Commit:** `feat(math): sym_eigen, one spectral decomposition for every eigenvalue gap (D-032)`.

### F2. La cinemática pasa al núcleo

1. `eigenprojections` y `eigenvals_refined` se reescriben sobre `sym_eigen` (`Pₐ = vₐ ⊗ vₐ`). Desaparecen
   `jacobi_band`, `jacobi_eigen` y el uso de `dEigenvalues_dTensor` en la cinemática.
2. `muscle_kinds` deja de tener `JACOBI_GAP_LO/HI` (D-018.19a queda reemplazada). Sin alias, como en D-027.
3. Verificación:
   - tests de `tests/kinematics` (verdad en cuádruple), sp y dp;
   - el caso de PROBLEMAS «F = I + 1e-6, error 0.67» y C = I + 1e-11;
   - los bordes de la antigua banda: el salto ≤ 1.1e-11 debe desaparecer;
   - tests `kin_*` de materiales;
   - tests target.
4. Registro de cambios numéricos: tabla de las diferencias máximas, frente al código anterior, de σ, tangente y
   R en los tests de cinemática (script en el scratchpad, no en el repo).
5. `check.sh` completo.

**Commit:** `refactor(kin)!: eigenprojections on sym_eigen; retire the Jacobi band (D-032)`.

### F3. Los criterios anisótropos pasan al núcleo

1. `spectral_decomposition` queda como envoltorio de `sym_eigen` (misma firma) o se reemplaza en las llamadas de
   CPB06, CPB06ex2 y Yld2004; se decide en D-032. Se conserva el escalado por potencia de 2 si el algoritmo
   elegido lo necesita (D-016.11), con `fp_scale` (D-031; `scale()` está prohibido en kernels).
2. Verificación:
   - `yield_laws_kinds`: gradiente y hessiana contra diferencias finitas, en sp y dp;
   - tests de materiales anisótropos de los 5 solvers;
   - **tabla de iteraciones de Newton** antes y después, por material y paso, para el informe;
   - tests target y GPU Intel. Es la fase de mayor riesgo con IGC.
3. `math_reference`: los registros de empates de `spectral_decomposition` dejan de valer, porque la base elegida
   en un autoespacio repetido no es única. Se reemplazan por tests unitarios de invariantes:
   - `V Λ Vᵀ`;
   - ortogonalidad;
   - suma de proyecciones del autoespacio.

   Se corre `mutation_check math --samples all` con el criterio de D-023: el registro sólo se borra si los
   unitarios detectan lo mismo.
4. Benchmarks GPU en Release, de los casos `cut_yld2004_ludwik`, `pz_cpb06ex2_bilinear` y `tm_jc_cpb06_ludwik`
   (`benchmarks/gpu/bench_cases.py`), contra la línea base de F0. Deben cumplir R7.
5. `check.sh` completo.

**Commit:** `refactor(yield)!: anisotropic criteria on sym_eigen (D-032)`.

### F4. El camino «sólo valores» (decisión)

Hoy la f de los criterios anisótropos usa `eigenvals` (forma cerrada, barata) y su gradiente usará los λ del
núcleo. La diferencia es de redondeo, pero sigue habiendo dos fuentes de autovalores. Opciones:

- **(a)** Dejar `eigenvals` como está: es lo barato, y conserva el problema marcado y el umbral absoluto 1e-10
  documentado en D-017.9.
- **(b)** La f usa los λ del núcleo. Hay una sola fuente y desaparecen el problema marcado y el umbral absoluto.
  Cuesta una descomposición más por evaluación de f, y cambia f al nivel del redondeo.
- **(c)** Como (a), y además corregir en `eigenvals` la rama `abs(1-d)` y el umbral absoluto.

Cualquier opción que no sea (a) toca el problema «Marcado, sin corregir», así que necesita autorización explícita
del usuario. Si se elige (b) o (c), el registro de `eigenvals` de `math_reference` se resuelve como en F3.3.

### F5. `dEigenvalues_dTensor` y `d2Eigenvalues_dTensor2` (decisión)

Sin consumidores en la librería después de F2: `d2` nunca los tuvo y `d1` los pierde con la cinemática. Opciones:
- mantenerlas como API Cayley-Hamilton documentada, sin vectores, con teoría en `Theory/`;
- reimplementarlas sobre el núcleo, `∂λₐ/∂S = Pₐ` (cambia resultados);
- retirarlas junto con `dsquare`, si nadie más lo usa: Drucker-Prager sí lo usa.

### F6. Limpieza y documentación

1. Borrar el código muerto:
   - restos de `spectral_decomposition`, de Jacobi y de `cross`, si ya no se usan;
   - constantes `EPS_TOL*`, si F5 retira las rutinas.
2. Áreas de mutaciones (`tests/mutations/math.toml`, `kin`, `yield`): fuentes al día y una corrida completa.
3. Documentación:
   - `DECISIONES.md`: D-032, más una nota de reemplazo en D-018.19a, D-016.11 y D-028;
   - `PROBLEMAS.md`: H3 de REANALISIS_MATH, la fila de la banda de Jacobi, la de los empates con FMA y la fila
     marcada si se resolvió en F4;
   - `pages/`, el README y la skill (`ubicacion.md`: «Autovalores, autovectores» apunta al núcleo);
   - `REANALISIS_MATH.md` (H2 resuelto).

### F7. Cierre

`check.sh` completo, benchmarks CPU y GPU en Release contra la línea base de F0, revisión independiente del diff
total de C (agente `reviewer`) e informe final al usuario: precisión, costo, iteraciones y decisiones tomadas.

## 6. Verificación común a todas las fases

1. `cmake --build build --target build_all_tests && ctest` mientras se edita.
2. Los tests espectrales y de cinemática con `-ffpe-trap=zero,overflow,invalid`. Los de math ya lo tienen; en F2
   se agrega a los de cinemática que no lo tengan.
3. `tools/mutation_check.py <área> --samples all --only <ids>` sobre lo nuevo, en las áreas `math`, `kin` y
   `yield`.
4. `tools/gen_snapshot.sh` donde no deba cambiar nada (F1, y los pasos de sólo limpieza de F6).
5. `tools/check.sh --quick` en cada paso; el completo al cerrar cada fase, con la GPU Intel.
6. Agente `reviewer` sobre el diff de cada fase antes de pedir el commit.
7. Ningún test se debilita ni se borra para que pase (`testing.md`). Si una tolerancia de material no alcanza,
   se mide, se explica y se decide.

## 7. Riesgos y mitigaciones

| Riesgo | Mitigación |
|---|---|
| IGC (GPU Intel) rechaza o compila mal el núcleo, como con `scale()` o `closest_point` | humo en GPU ya en F0.4; el núcleo pasa por los tests target en F1, antes de que lo use algún material; la implementación anterior queda en el historial para volver por fase |
| El costo en GPU sube (registros, iteraciones) | barridos fijos, sin divergencia; rotaciones explícitas, sin `matmul`; R7 con presupuesto acordado; si no se cumple, se usa KJ o K |
| Cambian las iteraciones de Newton de materiales anisótropos | tabla antes y después en F3; las referencias de solvers no se ven afectadas; los tests de secuencia son por propiedad |
| Autoespacios repetidos: la base no es única | sólo se prueban invariantes (proyecciones del autoespacio, `V Λ Vᵀ`), nunca componentes sueltas de V; `spectral_hessian` ya usa el límite de θ |
| La cinemática dp fuera de la banda cambia de Cayley-Hamilton al núcleo | sus tests tienen verdad en cuádruple; la tabla de diferencias de F2.4 documenta el cambio, que debería ser una mejora |
| Tocar el problema marcado de `eigenvals` sin permiso | aislado en F4, que necesita la autorización explícita del usuario |
| sp: escalas extremas (subnormales, desborde) | R4 con trampas; escalado por potencia de 2 si el algoritmo lo necesita |

## 8. Decisiones del usuario

- **Tomadas (2026-11-11):** orden de la carpeta primero (fase B); presupuesto R7 ~1.2x; unificar sólo lo que
  quepa en él.
- **Antes de la fase B:** confirmar la distribución de módulos de la sección 0.
- **Al cerrar F0 (D-032):** algoritmo (J, K o KJ), nombre del núcleo y constantes c de los criterios.
- **F3:** si `spectral_decomposition` queda como envoltorio o desaparece.
- **F4:** camino «sólo valores», (a), (b) o (c), con la autorización sobre el problema marcado.
- **F5:** destino de `dEigenvalues_dTensor`, `d2Eigenvalues_dTensor2` y `dsquare`.
