# Re-análisis 6: calidad de los tests (2026-11-10, rama `fypp/f9-readability`, HEAD `8390913`)

Aspecto 6 de 15 del re-análisis posterior al plan R (nuevo; no tiene análisis anterior). Preguntas: ¿qué cambios
del código no detecta ningún test?, y ¿aportan algo los 38 tests heredados (`legacy`) que no aporten los demás?
Método: `tools/mutation_check.py --samples all` en **todas** las áreas (≈ 4 300 mutaciones, build dp Debug, 6
carriles en worktrees desechables, ~45 min de pared). Para cada área con heredados, el grupo «referencia» fueron
los heredados y el grupo «unitario» todos los demás tests del área: `CRITERION: OK` significa que los heredados
no detectan nada que los otros no detecten (el criterio de D-023). Las mutaciones que ningún grupo detectó las
clasificaron tres agentes leyendo el código y los tests (EQUIVALENTE, SÓLO-SP, BAJO TOLERANCIA, HUECO, DUDOSO);
revisé a mano los huecos principales. Informes completos en `build-snapshots/re6/*.md` (no versionado).

Cambios de herramienta hechos para medir (sin commit aún): `mutation_check.py --reference/--unit` (sobrescriben
los grupos del área) y tres áreas nuevas que el plan R no cubría: `types` (tipos de tensores y
`muscle_tensor_common`), `linalg` (inversas, LU) y `materials` (plantillas de material).

## 1. Resultados por área

| Área | Mutaciones | Detectan los heredados | …y también los demás | Sólo los demás | Nadie | No compila |
|---|---|---|---|---|---|---|
| `tensors` (operadores) | 27 | 20 | 20 | 5 | 1 | 1 |
| `types` (nueva) | 172 | 143 | 143 | 20 | 9 | 0 |
| `linalg` (nueva) | 538 | 476 | 476 | 43 | 19 | 0 |
| `math` | 720 | 571 | 571 | 94 | 26 | 29 |
| `elas` | 163 | 85 | 85 | 78 | 0 | 0 |
| `hard` | 48 | 46 | 46 | 2 | 0 | 0 |
| `vp` | 407 | 61 | 61 | 317 | 29 | 0 |
| `visco` (sin heredados) | 229 | — | — | 214 | 7 | 8 |
| `mat_kin` (kin + material cinemático) | 427 | 105 | **104** | 168 | 139 | 15 |
| `yield` | 1 361 | 1 254 | **1 241** | 32 | 54 | 21 |
| `solvers` (grupo «referencia» = las 2 referencias que quedan) | 200 | 80 | 67 | 55 | 42 | 23 |
| `materials` (nueva; referencias de solvers) | 28 | 9 | 9 | 13 | 2 | 4 |

## 2. Los tests heredados

**En 8 de 10 áreas los heredados no aportan nada** que no detecten los demás (`tensors`, `types`, `linalg`,
`math`, `elas`, `hard`, `vp`; `visco` no tiene). Con el criterio de D-023 podrían retirarse. Excepciones:

- **`yield`, 13 mutaciones** que sólo detectan `test_muscle_yield_cpb06.f90` y `_cpb06ex2.f90`:
  - 12 en el **factor de normalización B de CPB06** (`muscle_yield_cpb06.fypp:134–136`). Un B erróneo sólo
    reescala f; diferencias finitas, homogeneidad y sp/dp no lo ven, y el único chequeo de valor absoluto
    (CPB06 isótropo = von Mises) usa C = I. El heredado compara con un valor de mpmath. Falta: **normalización
    uniaxial f(e₁⊗e₁) = 1** para cada conjunto (Cazacu 2006, ec. 12).
  - 1 en CPB06ex2 (`:163`, C′ toma la columna de C). Falta: **simetría de intercambio**
    f(C,k,C′,k′) = f(C′,k′,C,k), o C′ = 0 ⇒ CPB06.
- **`mat_kin`, 1 mutación**: `kin_pull_back_tangent` (`muscle_kin_finite.fypp:757`) **no tiene test unitario**
  (`test_kin_kinds` la importa y no la llama); sólo la cubre el ida y vuelta de `test_muscle_kinematics.f90:388`.
  Falta: pull-back(push-forward(C)) = C con J ≠ 1.
- `solvers`: las 13 que sólo detectan las referencias de `closest_point_time/perzyna` son las conocidas del plan
  R (camino de Newton, D-023).

Con esos tres chequeos agregados, **los 38 heredados (9 650 líneas, 49 % repetidas; aspecto 2) quedarían
redundantes** según el criterio del usuario.

## 3. Lo que nadie detecta (328 mutaciones)

| Clase | Cantidad | Comentario |
|---|---|---|
| Equivalentes | 141 | Comparación en un borde donde ambas ramas dan lo mismo; `tiny(2.0) = tiny(1.0)`; triángulo inferior de hessianos simétricos que nadie lee; ramas Φ = 0 reverificadas por el llamador; `<=`/`<` en `huge()` (todas las de `types`) |
| Sólo-sp | 97 | `polar_rotation` de `kin_finite` existe sólo en sp: el build dp de mutaciones no puede verlas (**límite de la herramienta**) |
| Bajo tolerancia | 30 | Criterios secundarios de convergencia, pivote en empates, ~1 ulp |
| **Huecos** | **53** | Sección 4 |
| Dudosas | 7 | Recorte Δγ < 0 de `closest_point` (¿alcanzable?), `phi == 1.0` en Yld2004 |

## 4. Huecos, por importancia (cada uno con el chequeo propuesto; VERIFICADOS los marcados ✔)

1. ✔ **Normalización de CPB06** (12) y **c11 ≠ 1**: ningún conjunto de prueba de CPB06 tiene c11 ≠ 1, así que
   `dev(1)*C1(1,1)` → `/C1(1,1)` sobrevive. Agregar un conjunto con c11 = 1.2 y la normalización uniaxial.
2. ✔ **`vp_newton_dstress_dep` de Johnson-Cook** (9): pública, sin test unitario; sólo la cubre la referencia de
   `closest_point_time`. Chequeo: igual a `vp_dstress_dep` fuera del piso; en el piso, la fórmula de la derecha.
3. ✔ **`kin_pull_back_tangent`** (1, sección 2).
4. ✔ **`kin_spin_tangent` en F = I** (`:1133`): el retorno temprano no se comprueba; con la mutación el término
   geométrico sería todo unos en la primera iteración. Chequeo: término geométrico = 0 y `kin_log_tangent` = A
   en F = I.
5. ✔ **Material con deformación inválida** (`material_kernel:151–152`): se comprueba el status, no que tensión
   y tangente salgan en cero. Chequeo: F = diag(−1, 1, 1) ⇒ status, `iters` = 0, tensión y tangente nulas.
6. **Tolerancia de las derivadas espectrales con un par casi nulo junto a uno dominante**
   (`spectral_derivs:132, 211, 392`): estados uniaxiales de CPB06/Yld2004. Los tests usan funciones polinómicas
   (la diferencia dividida se cancela limpia). Chequeo: φ = Σ exp(λ) en λ = (10, 3e-16, 1e-16): hessiano
   continuo con el de (10, 0, 0). (A confirmar corriendo el mutante: depende del redondeo.)
7. **Drucker-Prager `q/stress_norm`** (3): fallaría para |σ| ≲ 3e-5 en dp o hidrostáticos ruidosos grandes.
   Chequeo: homogeneidad exacta con λ = 2^±20 en dp e hidrostático no exacto en binario (p = 7.1e4).
8. **Ramas recortadas de MRK y RK** y derivadas nulas para ε̇ ≤ 1e-12 que promete la doc (`test_vp_limits`).
9. **`eig_c`/`proj_c` de `kin_finite_update`** (12): con el orden roto en la banda de Jacobi nadie lo ve
   (salida pública opcional). Chequeo: descendentes y proyectores = qₐqₐᵀ.
10. **Matrices no cuadradas** en `gauss_solve`, `gauss_seidel`, `mklsolve`, `lu_solve` (7): nadie prueba el
    rechazo (`info = −1`).
11. **Curvatura de CPB (`a <= 2`)** con un valor principal nulo y 2 < a ≤ 3 (el test salta ese caso con k ≠ 0).
12. Menores: `stress_hist` por defecto de `visco_frozen` (6, la librería lo sobrescribe), criterio `<` de
    Gauss-Seidel (1).

Además (aspecto 4): `fd_check_<m>` de los materiales Hill48 en `test_mat_return_mapping` se genera y no se
llama; **no es un hueco** (la tangente de esos materiales la comprueban los tests de secuencia en cada paso
plástico), sólo código muerto.

### 4.1 Huecos 1–5 cerrados (2026-11-10, Q1)
Tests nuevos con verdad independiente (sp y dp); cada mutación del hueco se volvió a correr con
`mutation_check --only` y quedó **detectada por los unitarios** (sin contar los heredados):
- `test_yield_kinds`: normalización uniaxial f(s e_x⊗e_x) = s en todos los conjuntos de CPB06 (Cazacu 2006,
  ec. 12; todos resultaron normalizados), conjunto 6 con c11 = 1.2 en todos los chequeos existentes, y
  CPB06ex2 con la simetría de intercambio (C, k) ↔ (C′, k′) y C′ = 0 ⇒ CPB06 (f, gradiente, hessiano):
  15/15 mutaciones (las 13 que sólo detectaban los heredados y las 2 de c11).
- `test_vp_limits`: `vp_newton_dstress_dep` de JC en el piso de tasa (forma cerrada y diferencia finita
  unilateral del flujo) y fuera de él (forma cerrada e igual a `vp_dstress_dep`), 4 juegos de parámetros: 9/9.
- `test_kin_finite`: pull-back y push-forward contra la fórmula de índices (verdad en precisión cuádruple, J ≠ 1)
  e ida y vuelta; con F = I, término geométrico nulo y tangentes logarítmicas = A: 2/2 (incluida la que sólo
  detectaba el heredado `kinematics-kinematics`).
- `test_mat_kinematics`: F inválida (det F ≤ 0) en los 11 materiales con F ⇒ `STATUS_INVALID_DEFORMATION`,
  `iters` = 0, tensión y tangente nulas (antes ningún test de material comprobaba ese status): 2/2.
Con esto **los heredados ya no detectan nada propio** en ninguna área (criterio de D-023 cumplido para Q2).
Queda sin test: F_n inválido en las cinemáticas espaciales (segunda llamada a `kin_update` del kernel).

## 5. Conclusión

- La batería es fuerte: de ~4 300 mutaciones, los tests detectan todas salvo ~50 huecos reales (≈ 1 %); el
  resto de lo que sobrevive es equivalente, sólo-sp o por debajo de las tolerancias.
- **Los tests heredados son casi totalmente redundantes**: aportan 14 detecciones en total, todas cubribles con
  3 chequeos con verdad independiente.
- Los huecos se concentran en **funciones públicas sin test propio** (`vp_newton_dstress_dep`,
  `kin_pull_back_tangent`), **casos borde** (F = I, F inválida, matrices no cuadradas, par casi nulo) y
  **parámetros de prueba poco variados** (c11 = 1 en todos los CPB06).
- Límite de la herramienta: el build de mutaciones es sólo dp; el código sólo-sp (`polar_rotation`) queda fuera.

## 6. Propuestas (no aplicadas)

| ID | Cambio | Esfuerzo |
|---|---|---|
| Q1 | **Hecho para 1–5 (sección 4.1).** Agregar los chequeos de los huecos 1–5 (y 6–11 según decida el usuario), con verdad independiente; volver a correr `mutation_check` de las áreas tocadas | 1–1.5 días |
| Q2 | **Hecho (2026-11-10, D-025)**, con los renombres espejo de `src` (`test_muscle_<x>[_<tema>]`). Con Q1 hecho, **retirar los 38 tests heredados** (criterio D-023) o conservarlos como documentación; **decisión del usuario** (la regla vigente dice «no se mueven ni se reescriben») | 0.25 |
| Q3 | Versionar las áreas nuevas (`types`, `linalg`, `materials`) y las opciones `--reference/--unit` | hecho en el árbol, falta commit |
| Q4 | `mutation_check` con build `sp;dp` opcional (o un área sp para `polar_rotation`) | 0.5 |
| Q5 | Borrar los `fd_check` no usados de `test_mat_return_mapping` (avisos, aspecto 4) | < 0.1 |
