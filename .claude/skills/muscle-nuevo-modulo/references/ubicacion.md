# Dónde vive cada cosa y qué no reescribir

Léela antes de escribir cualquier función que no sea la física propia del módulo. La lista completa y al día de lo
que existe la da `scripts/api_index.py <repo> [palabra ...]` (genéricos, tipos y constantes con su documentación;
`--all` agrega los específicos por kind). Busca por nombre **y por concepto** (`api_index.py <repo> eigen`,
`... deviatoric`, `... inverse`, `... rotation`, `... tolerance`): el nombre de lo que necesitas puede no ser el que
imaginas.

## Bloques que ya existen (no reescribirlos)

| Necesitas | Usa | Dónde |
|---|---|---|
| Suma, producto por escalar, contracción doble, simple, transformación, producto diádico | `+ - * /`, `.ddot.`, `.tdot.`, `.tdotsym.`, `.transform.` | `muscle_tensors` (`use muscle_tensors`) |
| Desviador, traza | `.dev.`; traza = `a .ddot. I` con `type(iden_2O) :: I` | `muscle_tensors` |
| Inversa, determinante, norma, cuadrado, transpuesta, conversiones 4O | `.inv.`/`inverse`, `det`, `norm`, `square`, `transpose`, `convert_2sym`/`convert_3sym` | `muscle_tensors` / `muscle_tensor_functions` (las funciones libres, para kernels) |
| Pasar entre tensor y matriz 3x3 (para `matmul` en un kernel) | `to_matrix(a)` (`2Osym`, `2O`), `from_upper(m)` (→ `ten_3D2Osym`) | `muscle_tensor_functions` |
| Comparar con tolerancia | `.approx.`, `is_approx(a, b, tol)` | tipos de tensores; tolerancias en `pages/tensors.md` |
| Identidades (2º y 4º orden, escaladas) | `iden_2O`, `iden_4O3T`, `iden_4O4T`, `iden_*S` | `muscle_tensors` |
| Autovalores, autovectores, autoproyecciones | `spectral_decomposition` (K, la única descomposición, precisa en toda brecha: D-032), `eigenprojections` (`P_a = v_a ⊗ v_a`), `eigenvals_refined`; `eigenvals` sólo valores, con umbral absoluto 1e-10 (no para desviadores chicos) | `muscle_math_eigen` |
| Derivadas de autovalores; gradiente y hessiano de una función de los autovalores | `dEigenvalues_dTensor`, `d2Eigenvalues_dTensor2`, `spectral_gradient`, `spectral_hessian`, `dsquare` | `muscle_math_spectral_derivs` |
| Sistema lineal pequeño en un kernel | `lu_solve`, `LU<n>SOLVE` | `muscle_math_lu_solve` (nunca LAPACK en kernels) |
| Inversa de matriz 3x3…9x9 con número de condición | `M33INV` … `M99INV`, `COND_INF` | `muscle_math_inverses` |
| Norma infinito de una matriz (máxima suma de fila), en un kernel | `NORM_INF` (bucles, sin temporales: D-030) | `muscle_math_inverses` |
| Medidas de deformación, conversiones de tensión, push-forward/pull-back de tangentes, rotaciones corrotacionales | `kin_*` | `muscle_kin_finite` (vía `muscle_kin_laws`) |
| Tolerancias y umbrales por kind | `TOL_APPROX_<k>`, `NEWTON_TOL_<k>`, `TOL_ELASTIC_<k>`, `SY_MIN_<k>`, `DGAMMA_MIN_<k>`, `EPS_ABS_<k>`, pasos de diferencias finitas | `muscle_kinds` |
| Códigos de error | `STATUS_OK`, `STATUS_INVALID_PARAMETER`, `STATUS_INVALID_TIME_STEP`, `STATUS_INVALID_DEFORMATION`, `STATUS_NOT_CONVERGED` | `muscle_status` |
| Ley de endurecimiento dentro de otra ley (compuesta) | recibir `sy`, `dsy_dep` ya evaluados (`hard_stress`, `hard_dstress_dep`) | `muscle_hard_laws`; listas `*_COMPOSED` |
| Derivadas por diferencias finitas | macros `FD_DERIV`, `FD_GRAD`, `FD_HESS`, `FD_JAC` | `tests/common/muscle_fd.fypp`: **sólo en tests**; la librería usa derivadas analíticas |
| NaN, `exponent()` o `scale()` en un kernel | `QUIET_NAN_<k>`, `fp_exponent`, `fp_scale` | `muscle_math_fp` (sin módulos IEEE, `exponent` ni `scale` en kernels: D-019.5b, D-031) |
| Rotación polar; (ln a − ln b)/(a − b) | `polar_rotation`, `log_divided_difference` | `muscle_math_matrix_functions` |

## Dónde va lo nuevo

| Lo que se escribe | Dónde |
|---|---|
| La física del modelo: su fórmula, sus derivadas, su inicialización | el módulo del modelo, `muscle_<prefijo>_<modelo>` |
| Una función auxiliar que usa los parámetros del modelo | el mismo módulo, privada, con nombre `<prefijo>_<modelo>_<cosa>_<k>` |
| Algo que comparten varios modelos **de una familia** | un módulo de soporte de la familia (precedentes: `muscle_yield_anisotropic` para CPB06/CPB06ex2/Yld2004, `muscle_visco_types`, `muscle_kin_finite`) |
| Álgebra general: autovalores, inversas, sistemas, funciones de matrices, derivadas espectrales | `src/math/`, en el módulo de su tema (`pages/math.md`: qué va en cada uno y las dos capas); módulo de soporte, `references/soporte.md` |
| Una operación de tensores (2D y 3D) | `src/tensors/` |
| Un umbral o tolerancia que depende de la precisión | `muscle_kinds`, por kind, documentado |
| Un umbral propio de un modelo (un piso de tasa de esa ley) | constante con nombre y documentación en el módulo del modelo, por kind |
| Un código de error nuevo | `muscle_status` (decisión: cambia la API) |

**Criterio rápido:** si la función no lee ningún parámetro del modelo (`prm%...`), no es del modelo. Si la usarían
dos módulos, no es de ninguno de ellos. Mover algo a `src/math/` o `src/tensors/` lo hace un módulo de soporte con
sus propios tests y su área de mutaciones; si cambia o agrega API pública, pregúntale antes al usuario.
