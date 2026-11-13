# Plan C, F0: medición de los candidatos espectrales (2026-11-12)

Fase F0 de `PLAN_ESPECTRAL.md`: precisión y costo de las descomposiciones espectrales actuales y de los candidatos
J, Jf, K y KJ, para decidir qué se unifica (decisión D-032). El usuario pidió el criterio original (el kernel no
debe crecer más de ~1.2x) **y también los casos hasta 2x**. No se cambió `src`.

## Qué se midió y cómo

- **Candidatos** (`benchmarks/math/spectral_candidates.fypp`, kernel: `pure`, `@:DEVICE()`, sin `matmul`):
  - **J**: Jacobi cíclico sobre S con rotaciones explícitas, salida temprana cuando un barrido no rota nada (como
    máximo 6 barridos).
  - **Jf**: el mismo con **4 barridos fijos**, sin salida dependiente de los datos (igual trabajo en cada hilo de
    GPU). Con 3 barridos fijos algunos cortes puros rotados no convergen en dp (fallan 6 tests); con 4 pasa todo.
  - **K**: Kopp mejorado: desviador escalado por potencia de 2 (exacto), autovalores en forma cerrada (Kopp 2008),
    vector del autovalor más separado por producto vectorial, par restante por una rotación plana sin desborde
    (D-028) y autovalores finales como cocientes de Rayleigh.
  - **KJ**: K y después, como mucho, 2 barridos de Jacobi sobre Vᵀ D V.
- **Precisión** (`benchmarks/math/bench_spectral.fypp`): S = Q diag(λ) Qᵀ en cuádruple, redondeado al kind. La verdad
  es la descomposición en cuádruple de la S redondeada, así que se mide sólo el error del algoritmo. Se usan:
  - 3408 casos por kind: tensiones a escala 1e-3, 1 y 2e5; C = I + δ/2 μ con δ = 1 … 1e-11; pares de autovalores
    con brecha 1e-1 … 1e-15 y 0, arriba y abajo; raíz triple; tres rotaciones; 3000 tensores aleatorios, la mitad en
    corte puro;
  - escalas extremas (dp 1e±30, sp 1e±18) aparte;
  - banderas de punto flotante con `ieee_exceptions`.
- **Costo por llamada**: un hilo, Release, sobre los mismos casos.
- **Costo del kernel en CPU**: en un worktree se reemplazó la implementación de `spectral_decomposition`,
  `eigenprojections` y `eigenvals_refined` por cada candidato (define de FYPP) y se midió el kernel de punto
  completo, con tangente:
  - materiales con criterio anisótropo: `bench_batch_gpu_<material>_<kind>` en CPU, 1 hilo, 20000 puntos;
  - cinemática finita: `benchmarks/math/bench_kin_kernel.fypp`, nuevo, con F = R(I + e) y puntos elásticos y
    plásticos;
  - `cut_vm_bilinear`, que no usa nada espectral, como control.
  En WSL2 una CPU virtual cae en núcleos P o E del anfitrión: las primeras corridas variaban hasta 3.6x (el control
  daba 1.9x). El método final alterna las cinco variantes caso por caso, 7 rondas, fijadas con `taskset`, y toma el
  mínimo. El control queda en 1.00 y mediana/mínimo ≤ 1.02.
- **Costo del kernel en GPU**: Intel UHD 770, ifx Release SPIR64, lote `cut_yld2004_ludwik` sp (el único material
  anisótropo que corre bien en Release en esa GPU). 1e5 y 1e6 puntos, 2 rondas, con límite de memoria y tiempo. La
  cinemática finita no se pudo medir en GPU: el banco de lote sólo admite materiales con entrada en deformación.
- **Tests**: la suite completa (133) en Debug con cada candidato en lugar de las rutinas actuales.

## 1. Precisión (máximos sobre los casos de escala física; ε = epsilon del kind; |S| = max |λ|)

| | e_λ /(ε\|S\|) | ortogonalidad /ε | e_P·brecha/ε | ln C (abs/ε) | ln C rel. al desviador | C = I + 1e-11 | escalas extremas | excepciones FP |
|---|---|---|---|---|---|---|---|---|
| actual criterios anisótropos, dp | 3.0e4 | 4.2 | 4.5e4 | 3.0e4 | **1.0 (se pierde)** | **falla** | **9e15 en 1e-30** | 0 |
| actual criterios anisótropos, sp | 1.9 | 3.8 | 3.4 | 2.3 | 3.6e-5 | bien | **1.7e7** | **6** (isótropo rotado; H3) + 49 extremas |
| actual cinemática, dp | 3.0e4 | – | 3.0e4 (55 aleatorios) | 3.0e4 | **1.0** | **falla** | **6e15 en 1e-30** | 0 |
| actual cinemática, sp | 5.9 | – | 7.8 | 4.3 | 8e-5 | bien | bien | 48 extremas |
| J / Jf, dp | 2.0 | 2.5 | 3.2 | 1.4 | 6e-5 | bien | bien | 0 |
| J / Jf, sp | 1.8 | 2.4 | 2.7 | 1.2 | 2.7e-5 | bien | bien | 0 |
| K, dp | 3.0 | 4.2 | 3.8 | 2.6 | 2.2e-5 | bien | bien | 0 |
| K, sp | 2.8 | 4.0 | 3.4 | 2.5 | 1.1e-5 | bien | bien | 0 |
| KJ, dp | 4.2 | 4.2 | 3.8 | 3.5 | 2.2e-5 | bien | bien | 0 |

- Los errores actuales de 3e4–4.5e4 en dp vienen todos de C = I + 1e-11: el umbral absoluto 1e-10 de `eigenvals`
  (D-017.9) devuelve una raíz triple y ln C pierde todo el desviador. La banda de Jacobi de la cinemática no lo
  evita, porque decide con esos mismos autovalores. Pasa igual con tensiones de escala 1e-30.
- Todos los candidatos cumplen R1–R5 con constantes de 2 a 4 ε, en todas las brechas, incluidas 1e-15 y 0. Ninguno
  levanta excepciones, ni en escalas extremas.
- KJ no mejora a K: los barridos de limpieza no reducen el error medido.

## 2. Costo por llamada (ns, un hilo; casos de arriba, la mayoría tensores aleatorios)

| | sp | dp | dp relativo a la actual de los criterios |
|---|---|---|---|
| `eigenvals` (sólo valores) | 93 | 101 | 0.61 |
| actual criterios anisótropos (`spectral_decomposition`) | 152 | 165 | 1.00 |
| actual cinemática (`eigenprojections`) | 523 | 193 | 1.17 |
| J (media de 3 barridos en dp, máximo 4) | 155 | 224 | 1.36 |
| K | 108 | 116 | **0.70** |
| KJ | 130 | 135 | 0.82 |

## 3. Costo del kernel completo en CPU (mínimo de 7, relativo a hoy)

| material | kind | hoy (ns/punto) | J | Jf | K | KJ |
|---|---|---|---|---|---|---|
| cut_vm_bilinear (control) | sp / dp | 145 / 416 | 1.01 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 1.01 / 1.00 |
| cut_yld2004_ludwik | sp / dp | 4868 / 7851 | 0.94 / 1.01 | 0.94 / 1.01 | 0.93 / 0.90 | 1.01 / 0.94 |
| pz_cha_yld2004_ludwik | sp / dp | 22411 / 32984 | 0.95 / 1.00 | 0.96 / 1.00 | 0.94 / 0.91 | 1.00 / 0.94 |
| pz_cpb06ex2_bilinear | sp / dp | 12396 / 18623 | 0.93 / 1.01 | 0.94 / 1.01 | 0.92 / 0.89 | 0.98 / 0.92 |
| tm_jc_cpb06_ludwik | sp / dp | 13327 / 18284 | 0.96 / 1.02 | 0.96 / 1.02 | 0.94 / 0.91 | 0.98 / 0.93 |
| kin_ml_vm_bilinear | sp / dp | 10506 / 8423 | 0.31 / 0.44 | 0.31 / 0.44 | 0.34 / 0.48 | 0.35 / 0.49 |
| kin_sl_vm_bilinear | sp / dp | 14208 / 10010 | 0.29 / 0.41 | 0.30 / 0.42 | 0.31 / 0.43 | 0.33 / 0.45 |
| kin_ul_vm_bilinear | sp / dp | 7724 / 5857 | 0.39 / 0.46 | 0.40 / 0.47 | 0.42 / 0.50 | 0.43 / 0.52 |

- En CPU **ningún candidato encarece ningún kernel** (máximo 1.02x).
- La cinemática finita se acelera 2 a 3.4 veces con todos. La rutina actual hace siempre Jacobi en sp, y en dp usa
  rotaciones con dos `matmul` 3x3 dentro de la banda; además el gradiente de los logaritmos la llama varias veces.

## 4. Costo del kernel en GPU (Intel UHD 770, `cut_yld2004_ludwik` sp, 1e6 puntos)

| | tiempo (s) | relativo a hoy | verificación GPU = CPU |
|---|---|---|---|
| hoy | 0.459 | 1.00 | correcta |
| K | 0.530 | **1.15** | correcta |
| KJ | 0.570 | 1.24 | correcta |
| J | 0.780 | 1.70 | **errores enormes** (salida 1) |
| Jf | 0.790 | 1.72 | **errores enormes** (salida 1) |

- J y Jf dan resultados incorrectos en la GPU Intel; en CPU son correctos (tests y verificación del banco). Es la
  fragilidad ya conocida de IGC (D-019.5c, `scale()`). Usarlos exigiría reducir el caso y quizá cambiar la
  formulación, sin garantía.
- K encarece este kernel un 15 % en la GPU Intel, frente a una reducción del 7 % en CPU.

## 5. Tests con cada candidato en lugar de las rutinas actuales (Debug, 133)

| | fallan | por qué |
|---|---|---|
| K | 1: `yield_kinds` | comparaciones exactas de dos caminos de Yld2004/CPB06 (evaluación conjunta y procedimiento de punto contra las funciones): uno toma los autovalores de `eigenvals` y el otro de la descomposición. Es la doble fuente de autovalores de F4, no un error del candidato. |
| J, Jf (4), Jf (5) | 2: `yield_kinds` (igual) y `math_reference` | los empates exactos eligen otra base del autoespacio repetido (previsto en el plan, F3.3). |
| Jf (3 barridos) | 6 | no converge en cortes puros rotados (descartado). |

Con K incluso `math_reference` pasa: en los empates exactos elige la misma base que la rutina actual.

## 6. Qué cabe en cada presupuesto

| Consumidor | ≤ 1.2x (criterio original) | ≤ 2x (pedido adicional) |
|---|---|---|
| Criterios anisótropos (CPU) | J, Jf, K, KJ (0.89–1.02x) | todos |
| Criterios anisótropos (GPU Intel) | **K (1.15x)** | K, KJ (1.24x); J/Jf (1.7x) entran por costo pero **dan resultados incorrectos** en IGC |
| Cinemática finita (CPU) | todos (0.29–0.52x: más rápidos) | todos |
| Cinemática finita (GPU) | no medido (el banco de lote no admite entrada F); en CPU la mejora es de 2 a 3.4x | – |

## 7. Recomendación para D-032

**K como núcleo único** (`sym_eigen`, nombre a confirmar), para los criterios anisótropos y la cinemática:
- **Presupuesto:** es el único que cabe en ≤ 1.2x en todo lo medido. En GPU Intel sube 1.15x el kernel anisótropo; en
  CPU lo baja un 6–11 %; la cinemática pasa a 0.34–0.50x.
- **Precisión:** del orden de J (2–4 ε), sin umbrales absolutos. Corrige los fallos actuales: C = I + 1e-11 en dp,
  escalas extremas y las excepciones de H3 en sp.
- **Tests:** pasa la suite salvo la doble fuente de autovalores de F4; ni siquiera cambian los empates de
  `math_reference`.
- **Simplificación:** desaparecen la banda de Jacobi (`JACOBI_GAP_*`), `jacobi_eigen` y las proyecciones de
  Cayley-Hamilton en la cinemática.

**Con presupuesto ≤ 2x**, las alternativas son:
- **KJ:** 1.24x en GPU, sin ganancia de precisión medible frente a K. No hay razón para preferirlo.
- **J/Jf:** son los más simples y precisos, pero en GPU Intel cuestan 1.7x y **fallan**. Sólo serían opción si se
  diagnostica y evita el problema de IGC, o si el objetivo deja de ser la GPU Intel (en NVIDIA/AMD no se midió).

**Decisiones para el usuario:**
1. Algoritmo: K (recomendado), KJ, J/Jf (con investigación de IGC), o no unificar.
2. Alcance: los dos consumidores (recomendado) o sólo uno (por ejemplo, sólo la cinemática, que gana 2–3x sin
   costo en GPU anisótropo).
3. F4, que se vuelve necesaria con cualquier opción: una sola fuente de autovalores para la f de los criterios
   anisótropos. Hay que elegir:
   - (b) la f usa los λ del núcleo: con K cuesta ~116 ns contra 101 de `eigenvals` en dp, por evaluación de f;
   - mantener `eigenvals` y que la evaluación conjunta y el procedimiento de punto también lo usen para los
     valores.
   Con (b) los criterios dejan de usar `eigenvals`, y con eso la rama marcada sin corregir deja de afectarlos;
   `eigenvals` sigue existiendo para quien sólo quiera valores.

## 8. Re-medición tras D-030 y D-031 (2026-11-12)

Con `fp_scale` (D-031) y `NORM_INF` (D-030) en el código, se repitieron las mediciones sobre el mismo arnés y se
sumaron dos candidatos:
- **Jnb**: Jacobi de 4 barridos sin ramas (todas las rotaciones, también las de elementos ya nulos).
- **Js**: Jacobi escalar de 4 barridos, rotaciones desenrolladas sobre escalares, escalado por potencia de 2.

La precisión de J, Jf, K y KJ no cambia (`fp_scale` da los mismos bits que `scale()`).

**Precisión de los nuevos:**
- Js: como J (2–3 ε).
- Jnb: falla en escala 1e-18 (error 1e10 ε): descartado.

**Tests (CPU, 133):**

| Variante | Fallan |
|---|---|
| actual | 0 |
| K, KJ | `yield_kinds` (F4) |
| J, Jf, Js | `yield_kinds` y `math_reference` (empates) |
| Jnb | los dos anteriores y `math_spectral_kinds` |

**Kernel en CPU** (mínimo de 7 rondas, relativo a hoy, sp / dp):

| | criterios anisótropos (4 materiales) | cinemática finita (3) |
|---|---|---|
| K | 0,90–0,92 / 0,87–0,91 | 0,31–0,42 / 0,43–0,50 |
| Js | 1,15–1,21 / 1,08–1,11 | 0,39–0,50 / 0,52–0,58 |
| J / Jf | 0,91–0,94 / 1,00–1,01 | 0,29–0,40 / 0,42–0,47 |
| KJ | 0,96–0,98 / 0,90–0,94 | 0,32–0,43 / 0,44–0,51 |
| Jnb | 3,4–5,9 / 1,2–1,4 | 0,83–0,95 / 0,85–0,95 |

**Kernel en GPU** (Intel UHD 770, Release, `cut_yld2004_ludwik` sp, 1e6 puntos, mínimo de 2 rondas):

| | tiempo (s) | relativo a hoy | fallos GPU vs CPU (muestra de 1000) |
|---|---|---|---|
| hoy | 0,145 | 1,00 | 0 |
| Js | 0,135 | 0,93 | 0 |
| K | 0,165 | 1,14 | 0 |
| KJ | 0,231 | 1,60 | 0 |
| J / Jf | 0,655 / 0,663 | 4,5 | 0 |
| Jnb | 1,089 | 7,5 | 0 |

- **Todos son correctos en GPU:** el fallo de J/Jf de la sección 4 era el error de IGC con la emulación FP64
  (`docs/fypp/igc_mre/`), que `fp_scale` saca del kernel.
- **El kernel de hoy es 3,2 veces más rápido que en la sección 4** (0,145 s contra 0,459 s), al no llamar a la
  emulación. Por eso los cocientes de J/Jf suben: 4,5x contra 1,7x.
- **Los otros tres materiales anisótropos** (`pz_cha_yld2004_ludwik`, `pz_cpb06ex2_bilinear`, `tm_jc_cpb06_ludwik`)
  siguen sin poder medirse en Release: IGC entrega la imagen sin sus kernels de lote, también con la rutina actual
  (D-019.8).

**Decisión (D-032):** K. Es el único que queda en ≤ 1,2x en todo lo medido y pasa `math_reference`. Js sería la
opción si la prioridad fuera la GPU Intel (0,93x, el más preciso), pero encarece la CPU en sp hasta 1,21x y cambia los
empates de `math_reference`.
