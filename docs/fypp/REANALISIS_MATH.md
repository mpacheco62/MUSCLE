# Re-análisis de `src/math/` (2026-11-11)

Pedido del usuario: la carpeta «está desordenada, hay varios archivos y algunos se pisan respecto de lo que
hacen». Análisis de sólo lectura sobre `fypp/f9-readability` en `29a5d31`: contenido, consumidores, capas,
solapamientos y una propuesta de reorganización. No se cambió código.

## 1. Inventario

| Archivo | Líneas | Contenido | Kernel | Consumidores en `src` (fuera de math) |
|---|---|---|---|---|
| `muscle_math_inverses` | 300 | `M33INV`, `M44INV` (adjunta + LU si κ > `ADJ_KAPPA_MAX`), `M55INV`, `M66INV`, `M99INV` (LU), `COND_INF`, `QUIET_NAN`, `ADJ_KAPPA_MAX` | sí | los 5 tipos de tensores (`.inv.`, `inverse`), `kin_finite` |
| `muscle_math_lu_inverse.fypp` | 171 | **include, no módulo**: macros `lu_factor`, `lu_subst`, `lu_inverse_body`, `lu_solve_body` | — | `inverses` y `lu_solve` |
| `muscle_math_lu_solve` | 134 | `LU<n>SOLVE` (n = 2..9) y `lu_solve` (frente con `select case`, sólo host) | sí | `closest_point_2` (`LU<n>SOLVE`); `lu_solve` sólo tests |
| `muscle_math_linsolvers` | 276 | `gauss_solve`, `gauss_seidel`, `mklsolve` (LAPACK `?gesv`), tamaño libre, `allocatable` | **no** (host-only) | **ninguno** (sólo tests) |
| `muscle_math_operations` | 172 | `eigenvals` (forma cerrada), `polar_rotation`, `log_divided_difference` | sí | `kin_finite`, `yield_cpb06`, `cpb06ex2`, `yld2004` |
| `muscle_math_derivatives` | 67 | `dsquare` (∂A²/∂A) | sí | `yield_druckerprager`; dentro de math, `d2Eigenvalues_dTensor2` |
| `muscle_math_spectral_derivs` | 624 | `dEigenvalues_dTensor`, `d2Eigenvalues_dTensor2`, `spectral_decomposition`, `eigenprojections`, `eigenvals_refined`, Jacobi (privado), `spectral_gradient`, `spectral_hessian`, `fp_exponent`, `cross`/`sym_outer` (privados) | sí | `kin_finite` (`dEigenvalues_dTensor`, `eigenprojections`, `eigenvals_refined`); yield CPB06/CPB06ex2/Yld2004 (`spectral_decomposition`, `spectral_gradient`, `spectral_hessian`, `fp_exponent`) |

Tests (`tests/math/`, 3 962 líneas): un archivo por módulo salvo `spectral_derivs` (4: `_kinds`, `_cases`,
`_eigenvals`, `_projections`), más `fd_macros`, `math_reference` y `math_target`. Áreas de mutaciones: `linalg`
(inverses, lu_inverse, lu_solve) y `math` (derivatives, operations, spectral_derivs, linsolvers).

## 2. Capas: la carpeta tiene dos niveles que no se ven en los nombres

```
  modelos (yield, kin, solvers)
        │
  math «sobre tensores»: operations, derivatives, spectral_derivs   (usan ten_3D2Osym, .ddot., iden_2O ...)
        │
  tensores (muscle_tensor_*)
        │
  math «bajo tensores»: inverses (+ lu_inverse.fypp), lu_solve      (matrices reales; los tensores los usan)

  host, al costado: linsolvers (allocatable, LAPACK)
```

Consecuencia práctica: nada de lo que usan los tensores (`QUIET_NAN`, inversas) puede ir a un módulo de math
que dependa de tensores. Cualquier reorganización tiene que respetar ese corte, y hoy sólo lo sabe quien lee
los `use`.

## 3. Hallazgos

### H1. Los nombres de los módulos no dicen lo que contienen (alto, legibilidad)

- `muscle_math_operations` es un cajón: autovalores, rotación polar y cociente de logaritmos. Su documentación
  de módulo sigue diciendo sólo «Eigenvalues of symmetric second-order tensors» (no se actualizó en D-028).
- `muscle_math_spectral_derivs` contiene derivadas espectrales, pero también la **descomposición** espectral
  (`spectral_decomposition`, `eigenprojections`, `eigenvals_refined`, Jacobi), una utilidad de bits IEEE
  (`fp_exponent`) y álgebra de vectores (`cross`, `sym_outer`). Es el 36 % de la carpeta.
- `muscle_math_derivatives` tiene una sola función (`dsquare`), y la derivada que más se usa
  (`dEigenvalues_dTensor`) no está ahí sino en `spectral_derivs`.
- `muscle_math_lu_inverse.fypp` tiene nombre de módulo pero es un include, y contiene también el cuerpo del
  solver (`lu_solve_body`), no sólo la inversa.

### H2. Lo espectral está partido en dos módulos y hay tres caminos a las autoproyecciones (alto, «se pisan»)

> **Resuelto (D-029 y D-032, 2026-11-12):** lo espectral está en `muscle_math_eigen` (D-029) y hay un solo camino a
> autovalores, autovectores y autoproyecciones: `spectral_decomposition` (K), que usan los criterios anisótropos y la
> cinemática (`eigenprojections`). `eigenvals` (sólo valores) y las derivadas de Cayley-Hamilton quedan como API pública.

| Necesidad | Rutina | Módulo | Método | Quién la usa |
|---|---|---|---|---|
| autovalores | `eigenvals` | operations | forma cerrada (Cayley-Hamilton, trigonométrica) | yield, kin |
| autovalores precisos con pares cercanos | `eigenvals_refined` | spectral_derivs | forma cerrada + Jacobi en la banda | kin |
| autovectores | `spectral_decomposition` | spectral_derivs | producto vectorial (Kopp) + **una** rotación plana | yield anisótropos |
| autoproyecciones | `eigenprojections` | spectral_derivs | Cayley-Hamilton fuera de la banda, Jacobi dentro | kin |
| autoproyecciones (= ∂λ/∂S) | `dEigenvalues_dTensor` | spectral_derivs | Cayley-Hamilton con reparto isótropo en raíces repetidas | kin, y `eigenprojections` por dentro |

Cada una existe por una razón documentada (los yield necesitan vectores para `spectral_gradient`/`hessian`;
la cinemática necesita proyecciones precisas cerca de C = I), pero quien busca «autovalores» los encuentra en un
módulo y «autovalores refinados» en otro, y dos estrategias distintas para el mismo problema (pares cercanos):
rotación plana única en `spectral_decomposition`, Jacobi cíclico en `eigenprojections`. Unificar los algoritmos
**cambia resultados** (decisión aparte, no recomendada ahora); juntar las rutinas en un módulo no.

### H3. La rotación plana de `spectral_decomposition` repite la fórmula que D-028 corrigió en Jacobi (medio, corrección posible)

> **Resuelto (D-032, 2026-11-12):** reproducido (el banco de F0 da 6 excepciones en sp con un isótropo rotado, y el
> test `diag(4, 2, 1)` con acoplamiento casi `tiny` aborta bajo `-ffpe-trap=overflow` con la fórmula anterior) y
> corregido: K usa la tangente sin desborde de D-028.

`spectral_decomposition` calcula `theta = (r - p)/(2 q)` con sólo `abs(q) > 0` y luego `theta*theta`
(`muscle_math_spectral_derivs.fypp`, rotación del bloque 2x2). Es el mismo patrón que desbordaba en
`jacobi_eigen` (D-028): con `q` mucho menor que `r - p` (≲ 1e-19 relativo en sp), `theta*theta` desborda; sin
trampa el resultado sigue siendo correcto (`tg` → 0), pero levanta la excepción de overflow, y con
`-ffpe-trap=overflow` el programa aborta. **No reproducido** (los tests actuales no generan ese `q`):
PLAUSIBLE. Corrección análoga a D-028 (`t = q/(r - p)` sin formar theta, o la rama `abs(2q) <= 1e-15 |r - p|`),
con un test que la provoque bajo `-ffpe-trap`. Sería una corrección de bug con su propio paso; cambiaría a lo
sumo el último bit de `tg` en esos casos.

### H4. Utilidades de kernel viviendo en módulos ajenos (medio)

- `QUIET_NAN` (NaN sin módulos IEEE) está en `inverses`; la usan los 5 tipos de tensores para la inversa
  singular.
- `fp_exponent` (el `exponent()` intrínseco desde los bits, para la GPU Intel) está en `spectral_derivs`; la
  usan los tres yield anisótropos.
- Ambas son lo mismo en espíritu: reemplazos de intrínsecos/IEEE que no existen en el dispositivo. Hoy un yield
  nuevo que necesite `fp_exponent` importa un módulo de derivadas espectrales.
- `cross` y `sym_outer` son privadas de `spectral_derivs`; `sym_outer(x, y)` es una operación de tensores
  (diádica simétrica de dos vectores → `ten_3D2Osym`), candidata a `muscle_tensor_functions` si un segundo
  módulo la necesita (hoy no).

### H5. Tres solvers lineales, uno solo usado (medio, API)

- `LU<n>SOLVE` (kernel, tamaño fijo): lo usa `closest_point_2`.
- `lu_solve` (frente genérico, host): sólo tests.
- `muscle_math_linsolvers` entero (`gauss_solve`, `gauss_seidel`, `mklsolve`): **sin consumidores en la
  librería**. `gauss_solve` repite el algoritmo de `lu_solve` con `allocatable`; `gauss_seidel` se documenta
  como «útil para estudiar convergencia»; `mklsolve` es el único punto de la librería que llama a LAPACK (y su
  referencia dependía de la LAPACK enlazada, PROBLEMAS). Es el único módulo `host-only` de math, el único con
  `allocatable` y el que aporta dos de las comparaciones de `-Wcompare-reals` pendientes. D-016 lo dejó «host-only en
  sp/dp», sin decidir si sigue siendo API pública. Retirarlo es un cambio de API (decisión del usuario).

### H6. `d2Eigenvalues_dTensor2` sin consumidores en la librería (bajo, API)

Sólo la usan tests. La que usan los modelos es `spectral_hessian` (vía autovectores). Es la única razón por la
que `spectral_derivs` depende de `dsquare`. Mantener (API pública documentada, teoría en `Theory/`) o retirar:
decisión.

### H7. Umbrales por kind repartidos con criterios distintos (bajo)

`EIG_EPS` y `LOG_SERIES` (operations), `EPS_TOL`, `EPS_TOL_SECOND`, `REPEATED_TOL` (spectral_derivs) y
`ADJ_KAPPA_MAX` (inverses) son `#:set` + constantes del módulo; `JACOBI_GAP_LO/HI` está en `muscle_kinds`. Todos
tienen nombre y documentación (cumplen `precision-kinds.md`), pero el lector no sabe dónde buscar. Criterio
sugerido: en el módulo si sólo lo usa ese módulo (todos menos `JACOBI_GAP`, que movimos por historia), y una
línea en la doc del módulo que los liste.

### H8. Estilo heredado (bajo)

- `eigenvals_3D2Osym` con sangría de 2, `.le.`, `**2.0_k`, líneas en blanco dobles y `  contains` desalineado
  (`muscle_math_operations`).
- `inverses` en MAYÚSCULAS (`PURE SUBROUTINE`, `REAL(...)`, `INTENT`), distinto del resto de la librería; los
  nombres `M33INV`... son API histórica, el estilo de las palabras clave no.

### H9. Tests que no espejan los módulos (bajo)

- Los tests de `eigenvals` (en `operations`) se llaman `test_muscle_math_spectral_derivs_eigenvals`.
- `test_muscle_fd_macros.fypp` prueba `tests/common/muscle_fd.fypp`, no algo de `src/math`.
- Nombres de ctest mezclados: `math-muscle_math_derivatives`, `math-muscle_math_inverses` frente a
  `math-math_lu_solve`, `math-math_spectral_kinds`.

## Estado (2026-11-11)

Fase B hecha (D-029): H1, H4, H8 y H9 resueltos; H2 queda ordenado (todo lo espectral en `muscle_math_eigen`) y su
unificación sigue en `PLAN_ESPECTRAL.md`; H3 (desborde posible de la rotación plana), H5 (`linsolvers`), H6
(`d2Eigenvalues_dTensor2`) y H7 (umbrales) siguen abiertos.

## 4. Propuesta

### A. Arreglos sin tocar la API (recomendado, primero)

1. Doc de módulo de `operations` al día; estilo de `eigenvals` e `inverses` (commit de formato separado,
   `gen_snapshot` = mismo código salvo comentarios/espacios).
2. H3: test con `-ffpe-trap` que provoque el desborde y la corrección análoga a D-028 (paso propio, D nueva).
3. Tests: nombres de archivo y de ctest coherentes con D-025.

### B. Reordenar por tema (renombres sin alias, como D-027; resultados idénticos bit a bit)

| Módulo propuesto | Contenido | Capa |
|---|---|---|
| `muscle_math_inverses` (+ `muscle_math_lu_bodies.fypp`, renombre del include) | inversas, `COND_INF` | bajo tensores |
| `muscle_math_lu_solve` | sin cambios | bajo tensores |
| `muscle_math_fp` (nuevo) | `QUIET_NAN`, `fp_exponent` (reemplazos de IEEE/intrínsecos para el dispositivo) | bajo tensores |
| `muscle_math_eigen` (nuevo) | `eigenvals`, `eigenvals_refined`, `eigenprojections`, `spectral_decomposition`, Jacobi | sobre tensores |
| `muscle_math_spectral_derivs` | `dEigenvalues_dTensor`, `d2Eigenvalues_dTensor2`, `spectral_gradient`, `spectral_hessian`, `dsquare` (absorbe `derivatives`) | sobre tensores |
| `muscle_math_tensor_functions` (renombre de `operations`) | `polar_rotation`, `log_divided_difference` | sobre tensores |
| `muscle_math_linsolvers` | retirar, o dejar como está (host) | host |

Costo: cambian los `use` de ~12 módulos de `src` y de los tests, `tests/mutations/*.toml`, la doc FORD, la
skill (`ubicacion.md`) y `pages/`. Sin cambio numérico (`gen_snapshot`: sólo nombres de módulo y símbolos).
Una alternativa más chica es B parcial: sólo `muscle_math_eigen` + `muscle_math_fp` y la doc al día, que
resuelve lo que más «se pisa» (H2, H4) sin tocar el resto.

### C. Unificar los algoritmos espectrales (no recomendado ahora)

Un solo método robusto para pares cercanos (Jacobi para todo, o Kopp + rotación para todo) eliminaría la
duplicación de H2 de raíz, pero cambia resultados de yield o de cinemática y su costo en GPU. Sólo si aparece un
motivo (precisión o rendimiento medidos).

## 5. Decisiones abiertas

1. ¿A? (doc, estilo, H3, tests).
2. ¿B completa, B parcial o nada?
3. `muscle_math_linsolvers`: ¿retirar o mantener?
4. `d2Eigenvalues_dTensor2`: ¿mantener o retirar?
