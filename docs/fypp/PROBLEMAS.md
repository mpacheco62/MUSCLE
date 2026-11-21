# Problemas encontrados durante la migración FYPP

Registro único de los bugs, anomalías y limitaciones encontrados al migrar MUSCLE (código original,
compiladores y herramientas). Cada entrada indica dónde está, su estado y la decisión que la cubre
en [`DECISIONES.md`](DECISIONES.md). Se actualiza al terminar cada tarea, junto con `PLAN.md`.

Estados:
- **Corregido**: el comportamiento cambió (cambio aprobado por el usuario; golden actualizado).
- **Documentado**: se conserva bit a bit en dp y se documenta en el FORD del módulo; no se corrige.
- **Pendiente**: hay que resolverlo en una fase posterior (se indica cuál).
- **Externo**: bug de un compilador o herramienta; se evita con un rodeo.

## Tensores (Fases 2–4)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Umbrales de `is_approx`/`.approx.` escritos como literales real32 (`1e-12` ≈ 9.9999999600e-13) | todas las familias de tensores | Corregido: literales del kind | D-010, D-011 |
| Inversa 3D con overflow en la diagonal para tensores grandes | `inverse` de 3D2Osym | Corregido | D-011 |
| Errata en la salida formateada (UDTIO) de 3D2O y mensaje de error de `iden_2O` | `muscle_tensor_io` | Corregido | D-011, D-014.5 |
| Diagonal de `.tdotsym.` 3D: `0.5*(a_I*b_I + b_I*a_I)` desborda si a_I·b_I > huge/2 | `contraction`/`dyadic` | Corregido: `a_I*b_I` como en 2D | D-015.6 |
| `.ddot.` 4O3sym:4O3sym con orden de sumas distinto en 2D (matmul) y 3D (orden de sympy) | `contraction_double` | Corregido: una sola regla (≤ 1 ulp de la escala) | D-015.8 |
| Componente 33 de `.transform.` 2D: `(F33*F33)*S33` en vez de la regla general | `transform` | Documentado (cambiarla movería ~1 ulp) | D-015.7 |
| Con FMA (`-mfma`, `-march=native`) el golden bit a bit de `.ddot.` deja de coincidir | `contraction_double` | Documentado (el proyecto no usa FMA) | D-015.9 |
| Inversa por adjunta en sp poco precisa para κ∞ ∈ [1e2, 1e3] | `muscle_math_inverses` | Corregido: umbral adjunta/LU sp = 1e2 | D-013 |

## Math (Fase 5)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Derivadas numéricas falsas (= 0) con argumento nulo: paso mínimo 1e-40 hace f(±h) = f(0) | `muscle_math_derivatives` | Corregido: h = √ε en x = 0 (1ª derivada), paso relativo en la 2ª | D-016.4 |
| Paso de `derivative_escalar_tensor3x3sym` relativo a \|f(x)\| (derivada 0 si f = 0) | `muscle_math_operations_host` | Corregido: relativo a \|x\| | D-016.4 |
| Umbral de `eigenvals` escrito como literal real32 (1.0000000134e-10) | `muscle_math_operations` | Corregido: `1.0e-10_dp`; sp 1e-9 | D-016.5 |
| Umbrales de las derivadas espectrales absolutos (unidades de tensión²): en sp, \|H\| hasta 32 con escala 1000 | `muscle_math_spectral_derivs` | Corregido: relativos a max λ² | D-016.6 |
| `spectral_decomposition` sp pierde la normalización con cortes ~1e-11 (\|c\|² subnormal) | `muscle_math_spectral_derivs` | Corregido: escalado exacto por potencia de 2 (desde D-032 se escala el desviador de K) | D-016.11, D-032 |
| Umbral absoluto 1e-10 de `eigenvals`: descarta cortantes chicos sobre una hidrostática | `eigenvals` (afectaba a CPB06, CPB06ex2, Yld2004) | Corregido para los criterios: su f toma los autovalores de `spectral_decomposition` (K, sin umbral; test `test_small_shear`); `eigenvals` lo conserva | D-017.9, D-032 |
| `spectral_decomposition` (rotación plana, antes de D-032): la tangente `sign(θ)/(|θ| + sqrt(θ² + 1))` desborda cuando el corte del par es casi `tiny` (la fórmula que D-028 corrigió en Jacobi); en sp, con un tensor isótropo rotado, 6 casos del banco levantan excepciones (H3 de `REANALISIS_MATH.md`) | `muscle_math_eigen` | Corregido: K usa la tangente sin desborde de D-028 y el desviador escalado; 0 excepciones en el banco | D-032 |
| `eigenvals` divide por cero (1/0 = Inf, resultado correcto) con la raíz mayor repetida; aborta con trampa de división por cero | `eigenvals` | Corregido (Fase 8): con `n2 <= n1/huge` va directo al límite, sin dividir; mismos bits. Test con trampas FPE | D-019.6 |
| Un argumento no nulo pero diminuto conserva el paso mínimo histórico de las diferencias finitas | `include/muscle_fd.fypp` | Documentado | D-016.4 |
| Macros FD vs API host: los últimos bits cambian si el compilador integra la función (≈4e-11 con -O3) | `include/muscle_fd.fypp` | Documentado (comparar con tolerancia en la Fase 6) | D-016.9 |
| Golden de `mklsolve` depende de la LAPACK enlazada | `tests/reference/math` | Documentado (la cabecera anota la librería) | — |
| `M33INV`/`M44INV` aceptaban la inversa adjunta con un NaN: `maxval(sum(abs(.), dim=2))` de gfortran salta las filas NaN, así que \(\kappa_\infty\) salía finito (p. ej. A = [2b b 0; b b 0; 0 0 1/b] con b² desbordado: el cofactor c33 = Inf − Inf es NaN y det(A) es finito) y no se recurría a LU | `muscle_math_inverses` | Corregido: `NORM_INF` devuelve NaN con una fila NaN y el adjunto se rechaza; test `test_adjugate_nan` (encontrado por el revisor de D-030) | D-030 |
| `M<n><n>INV` y `lu_solve` devolvían `OK_FLAG = .true.` con una entrada ±Inf: la LU sólo rechazaba un pivote nulo o NaN (`abs(p) > 0` es verdadero para Inf); la salida era NaN (Inf en (n,n) de 3x3 y 6x6) o finita sin sentido (Inf en (1,1)). Un Inf por desborde aguas arriba pasaba como iteración válida de un Newton | `muscle_math_lu_bodies`, `muscle_math_lu_solve` (n = 1) | Corregido: el pivote debe ser además finito; tests con ±Inf en cada posición, n = 1..9 y 3..9 (encontrado por el banco de nivel 0 de la campaña de verificación) | D-034 |
| LU con bucles e índice de pivote en tiempo de ejecución: en la GPU Intel la matriz baja a memoria privada (no se indexan registros con un índice variable) y el 6x6 sp tarda 144 ns por matriz; con `M66INV` el lote de `tm_jc_cpb06_ludwik` sp no tenía kernel en Release | `muscle_math_lu_bodies` | Mejorado (rendimiento): LU desenrollada con pivoteo por compare-and-swap hasta 7x7 (memoria privada 0, 23.6 ns el 6x6 sp) | D-035 |

## Endurecimiento (Fase 6)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Ludwik: las derivadas usan `ep + 1e-10` y la tensión usa `ep` (no son sus derivadas exactas) | `muscle_hard_ludwik` | Documentado (offset 1e-10 también en sp) | D-017.4 |
| Swift con e0 = 0 (Hollomon): derivadas ±Inf o NaN en ep = 0 | `muscle_hard_swift` | Documentado (hace falta e0 > 0; no se valida) | D-017.4 |
| Voce modificada: la doc describía mal `Sy` y omitía términos; texto en castellano en la doc | `muscle_hard_vocemod` | Corregido (sólo doc) | — |

## Viscoplásticas (Fase 6)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| `dstress_dep` omite ∂σ/∂ep | RK | Documentado | D-017.7 |
| `flow_stress = 0` con epd ≤ 0 (tensión nula en el primer incremento) | RK (JCfull: corregido, piso `epmin`, D-043) | Documentado | D-017.7, D-043 |
| Inf en epd = epdmax y NaN por encima (con tasas realistas) | MRK, RK | Documentado | D-017.7 |
| Salto en epd = 0⁺; `dstress_dep` = 0 con epd ≤ 1e-12 y finita donde σ es NaN/Inf | MRK | Documentado | D-017.7 |
| Con F = 1 − χ1·L ≤ 0 la flujo es plana pero las derivadas son grandes (χ2 > 1) | MRK | Documentado | D-017.7 |
| `dstress_dep` ignora los recortes de `flow_stress` (suelos de ep distintos, X ≤ 0.999) | NNL | Documentado | D-017.7 |
| Derivadas evaluadas en el suelo donde la flujo es plana | VA | Documentado | D-017.7 |
| División por cero si Y(ep) = 0 | Perzyna | Documentado | D-017.7 |
| MRK en sp: `1 - exp(-n·ep)` pierde precisión con ep chico (cancelación, 4e-5 en ep = 1e-4) | MRK | Documentado | — |
| El adaptador OO de JCfull evalúa `hard_law` con epd ≤ 0 (el original no): exige `hard_law` asignado | JCfull | Documentado | — |

## Criterios de fluencia (Fase 6)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Derivadas analíticas NaN (0/0) con tensión equivalente nula | von Mises, Hill48, CPB06, CPB06ex2, Yld2004 | Corregido: gradiente y hessiana = 0 en f = 0 | D-017.9 |
| No exactamente insensible a la presión: `linear_transform` aplica L = C·P_dev directamente sobre σ y, por redondeo, las filas del bloque normal de L suman ~1e-17 en vez de 0; con una hidrostática pura f ≈ 7e-17·p (AA2090) y el gradiente sale O(1) en dirección arbitraria. Irrelevante con cualquier parte desviadora. Corregirlo (aplicar L sobre diferencias como CPB06, o forzar filas de suma 0 en `init`) cambiaría bits en dp | Yld2004 | Documentado (decisión del usuario 2026-11-01: sólo documentar) | D-017.9 |
| Sin `init` (K = 0) las derivadas numéricas dan NaN/Inf | Drucker-Prager | Documentado | D-017.9 |
| En el ápice la derivada analítica (subgradiente) y la numérica difieren | Drucker-Prager | Documentado (no diferenciable) | D-017.9 |
| Hessiana con k = ±1 y un autovalor nulo a ~45% de la numérica (valor elegido en `curvature`) | CPB06 | Documentado | D-017.9 |
| En sp, Φ subnormal con \|σ\| ≲ 3e-6 (a = 8): derivadas con error o Inf | CPB06 | Corregido: λ reescalado por potencia de 2 de max\|λ\| | D-017.10c |
| Hessiana con seq³: Inf en sp con seq ≲ 1.5e-13 y término perdido con seq ≳ 7e12 | von Mises | Corregido: forma estable con n = dev/seq | D-017.10d |
| Parámetros derivados (B, L) públicos: el tipo de kernel se podía construir inconsistente | CPB06 (y CPB06ex2, Yld2004, Drucker-Prager) | Corregido: componentes privados + `yield_<m>_init` | D-017.10b |
| Con P no semidefinido, seq² = 0 con P:σ ≠ 0 daba Inf; ahora 0 (seq² < 0 da NaN) | Hill48 | Documentado | — |
| Hessiana con q⁴: en sp desbordaba con q ≳ 4.3e9 y salía finita pero incorrecta (10–50%) | Drucker-Prager | Corregido: forma estable con n = S/q | D-017.11d |
| `stress_eq` y gradiente con r³/q³: en sp valen sólo para 2e-13 ≲ q ≲ 7e12 | Drucker-Prager | Documentado (forma estable cambiaría más bits dp) | D-017.11d |
| Eje hidrostático en sp: con umbral 1e-9 un hidrostático no exacto en binario no se detectaba | Drucker-Prager | Corregido: umbral sp 1e-5 | D-017.11b |
| `pow` no es exactamente homogéneo: el reescalado exacto por potencia de 2 mueve f ≤ 12 ulp con exponentes no enteros | CPB06ex2, Yld2004 | Aceptado (hasta 1.02e-14 en una hessiana, excepción aprobada) | D-017.11a |
| Los ejecutables de test no enlazan con offload (`GOMP_offload_register_ver` indefinido): faltan flags de enlace en CMake | todos los tests | Pendiente (Fase 8) | — |
| CPB06: `prm%C1` (9 reales por material) sólo lo usa `yield_cpb06_init` para calcular B y L; desde D-032 el kernel ya no lo lee (la f usa `L`), pero viaja al dispositivo con el tipo | `muscle_yield_cpb06` | Documentado (quitarlo cambia el tipo: decisión aparte) | D-032 |

## Elasticidad (Fase 6)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Desviador por dos caminos (3D analítico; 2D y `_transient` por `.dev.`), distintos en redondeo | `Elasticity_linear` | Corregido: un solo camino analítico (~1e-16) | D-017.12a |
| ν = 0.5: división por cero en `set_parameters` (λ = +Inf, tensión NaN, tangente Inf); ν no se valida | `Elasticity_linear` | Documentado (exige ν < 0.5) | D-017.12b |
| Rutas `_transient` y `stress_and_update_3D` ignoran `dstrain`, `dtime` e `internal`; faltan `stress_and_update_2D` y `dstress_dstrain_dev` 2D | `Base_elasticity`/`Elasticity_linear` | Documentado | D-017.12b |

## Solvers (Fase 7)

Nota: las entradas de las Fases 5–6 que citan la API OO (`_transient`, `unpack_from_fea`,
`update_incremental`, `Elasticity_frozen`, `Base_*`) describen código retirado en la Fase 7 (D-018.18);
se conservan como historia.

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Tangente numérico con deformación nula: el paso debe ligarse al material (en sp √ε ≈ 3.5e-4 puede acercarse a la fluencia) | `closest_point`, `_time`, `_perzyna` | Resuelto: los kernels usan tangente analítico; DF sólo en tests (D-018.6) | D-016.4, D-018.6 |
| Al agotar iteraciones `iters` vale `iter_nw+1` (perzyna usa `min`) | `closest_point` | Corregido (`min(iters, iter_nw)`) | D-018.10a |
| `strain_p` sin inicializar en el tangente del caso elástico | `closest_point_time` | Documentado (se porta tal cual) | D-018.2 |
| Chequeo elástico con epd = 0 y sin guarda: σ_flujo(ep, 0) = 0 da 0/0 = NaN y deformación nula entra al Newton (jcfull: 200 iteraciones sin converger en pasos elásticos/descarga) | `closest_point_time` | Corregido (OO y kernel) | D-018.10b, D-018.12d |
| `strain_p` no recibe omega (a diferencia de `closest_point`); `init` no reinicia omega | `closest_point_time` | Documentado | D-018.2 |
| Sin convergencia devuelve tensión 0 (`data%stress` no se escribe) y tangente numérico 0 | `closest_point_time`, `closest_point_2` | Documentado (status ≠ OK) | D-017.14f |
| jcfull en sp: paso 1 sin converger, después ep difiere ~6e-4 de dp y el tangente por DF da NaN (flujo 0 a tasa nula) | material `time` + JC completo | Corregido: piso de tasa `epmin` (D-043); sp se comprueba | D-018.12b, D-017.7, D-043 |
| Perzyna en sp no converge con incrementos de ~40 deformaciones de fluencia (4e-2), aun con el piso de redondeo | `perzyna` sp (material `kin_ul_pz_*`) | Documentado; el test usa incrementos ≤ 3e-2 | D-018.12a |
| cp2 con endurecimiento blando converge (status OK) a otra solución (5 % de diferencia) en el chequeo F→I | `cp2` + cinemática | Documentado (sin dγ ≥ 0) | D-018.10c |
| Incrementos de ~1000 deformaciones de fluencia (cpb06) y exponentes extremos de Perzyna no convergen en sp | `closest_point`, `perzyna` sp | Documentado; las pruebas físicas corren sólo en dp | D-018.18g |
| dp: cinemáticas logarítmicas (kin_ml_*, kin_sl_*) con F = I + 1e-6 dan σ vs ddsdde:ε con error relativo 0.67 (6e-5 con 1e-5). Causa: con brechas < ~3e-6 la forma cerrada toma las raíces como repetidas (P = I/3) y ln C pierde la parte desviadora; σ está mal, no el tangente | log material/espacial dp | Corregido (Jacobi en la franja [1e-12, 1e-2)); desde D-032 sin franja: `eigenprojections` sobre `spectral_decomposition` (error relativo de H ≲ 10 ε/\|F − I\|, el límite de formar C) | D-018.19a, D-032 |
| `kin_tl_tm_vm_jc` no converge (200 iteraciones) con 1e-4 de deformación, σ ≈ 21 < y0 = 200, en ambos kinds | material TL + time + JC (también sin cinemática) | Corregido: era plástico (flujo a tasa mínima 15.8 MPa) y el Newton se estancaba porque la derivada de tasa es 0 en el piso; el solver usa la derivada por la derecha | D-018.19b |
| `tm_nnl_vm` con deformación ×40 no converge (200 iteraciones en dp, 177 en sp) | `time` + NNL | Pendiente (Fase 9); excluido del chequeo de status del smoke test | D-018.19 |
| MRK: NaN desde el primer paso plástico (probable epd > epdmax de la ley, ver viscoplásticas) | `closest_point_time` | Documentado | D-018.2 |
| Tangente analítico vs numérico: RK ~8 %, NNL 3.5e-4 (pasos convergidos) | `closest_point_time` | Documentado | D-018.10d |
| Sin restricción dγ ≥ 0: incremento grande bilineal + von Mises converge a ep < 0 (−1.84e-2; cp da +2.00e-2); no converge con cpb06/cpb06ex2 (incremento grande) ni yld2004-2 (NaN) | `closest_point_2` | Documentado | D-018.10c |
| Chaboche n=0.5 + Swift + DP en compresión: 4/6 pasos plásticos sin converger; tangentes con diferencias 1.4–1.6e-5 (Perzyna n=0.5 Hill48, n=5 DP) | `closest_point_perzyna` | Pendiente (Fase 9); con D-018.12a persiste en el caso 6 (tangente DF vs analítico 0.24) | Fase 7 golden |
| `closest_point_2`: la fila de consistencia del jacobiano 7x7 (`a(7, 1:6) = (-C:df)%vals`) no multiplica por 2 las componentes de corte (la derivada de f respecto de la componente Voigt j es −(C:n):e_j con peso 2 en los cortes): el punto fijo es correcto, pero la convergencia es lineal con cortes (paso Hill grande: 50 iteraciones frente a 7 de `closest_point`; retorno radial no exacto tras una iteración, 1e-10) | `src/solvers/muscle_closest_point_2_body.fypp` | **Marcado, sin corregir** (usuario, 2026-11-09: queda pendiente; corregirlo cambia resultados). `test_mat_return_mapping` lo excluye de la convergencia cuadrática y del conteo de iteraciones | ESTADO_OPTIMO R6 |
| `closest_point_time` + Johnson-Cook: con poca sobretensión (4·TOL_ELASTIC, solución en el piso de tasa Δγ/dt < epmin) no converge en 200 iteraciones: la derivada por la derecha σy·C/(epmin·dt) (D-018.19b) hace pasos mucho menores que la solución | `src/solvers/muscle_closest_point_time_body.fypp`, `muscle_vp_jc.fypp` | **Marcado, sin corregir** (usuario, 2026-11-09: queda pendiente; corregirlo cambia resultados). `test_mat_return_mapping` sólo comprueba «no elástico» en ese caso | ESTADO_OPTIMO R6 |


## Viscoelasticidad (Fase 6)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| γ = (1−e^{−x})/x cancela ~3.7e-12 justo encima del cambio a la serie (x = 1e-5) | Maxwell | Corregido: umbral 1e-4 en dp | D-017.13a |
| dt < 0 sin validar (Inf/NaN); Kelvin aborta con `error stop` (inválido en GPU) | Maxwell, Kelvin | Corregido: `status` en kernel; OO aborta en ambos | D-017.13b |
| No valida parámetros negativos ni ν = 0.5 | Maxwell | Documentado | D-017.13c |
| dt = NaN pasa las validaciones (`dt < 0` / `dt <= 0` son falsas con NaN) con `STATUS_OK`, igual que el original | Maxwell, Kelvin | Documentado (cubrirlo exigiría `.not. (dt >= 0)` y cambiaría el OO) | D-017.13 |
| `visco_step_coefficients` (antes `visco_coefficients`) chocaba con el tipo `Visco_coefficients` (Fortran no distingue mayúsculas) | familia visco | Corregido: renombrado | D-017.13d |
| Con dt muy chico la tensión es diferencia de términos ∝ ζ/dt (cancelación) | Kelvin | Documentado | D-017.13c |
| `stress_2D` ignora yz/xz de la tensión histórica; modelos y estado sólo 3D | `Elasticity_frozen` | Documentado | D-017.13c |
| `_transient` ignoran dstrain/dtime/internal; `unpack_from_fea` exige `init` previo | Base/estado | Documentado | D-017.13c |

## Cinemáticas (Fase 6)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Umbral absoluto 1e-12 sobre los autovalores de C: con F casi singular R no es ortogonal y ln λ se pierde, sin aviso | cinemáticas finitas | Corregido en kernel (`status`); OO igual | D-017.14a |
| det F ≤ 0 sin validar (J negativo, R impropia); `ok` de `inverse` no propagado (F singular → NaN) | cinemáticas finitas | Corregido en kernel (`status`); OO igual | D-017.14a |
| R pierde precisión cerca de C ≈ I (\|RᵀR − I\| ≈ 3.3e-11 con F ≈ I + 1e-3) | polar (vía `dEigenvalues_dTensor`) | Documentado | D-017.14b |
| Tangentes ignoran `stress_spatial`/`stress_material`; logarítmica material sin derivada de ln C en el tangente espacial; `update_incremental` aditivo e ignora dt | cinemáticas | Documentado | D-017.14b |

## GPU / compiladores (Fase 8)

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| ifx 2026.2: código de dispositivo roto (LLVM «Instruction does not dominate all uses») para `matmul` vector×matriz en `declare target` | `muscle_tensor_ops_contraction_double.fypp:108` | Corregido: sumas explícitas en todos los builds | D-019.4c |
| ifx rechaza `#line` entre las mitades de un literal partido por fypp | plantillas | Corregido: `--line-numbering-mode=nocontlines` | D-019.4d |
| ifx 2026.2: error interno en `call p%init(...)` si el tipo llega por `muscle_tensors` con dos nombres (`ten_3D2Osym` y `ten_3D2Osym_dp`) | 5 tests `tensors-types-*_kinds` | Corregido (2026-11-14): los alias históricos viven en `muscle_tensor_aliases`; la causa es importar un tipo con dos nombres en un mismo módulo (se aisló con reproductores mínimos de ~5 líneas contra los módulos de la biblioteca, que no están versionados; pendiente de reportar a Intel) | D-019.4f, D-036 |
| ifx 2026.2 (Windows, CPU): los tests con igualdad bit a bit fallan con el modelo de punto flotante por defecto (`math_inverses`, `math_fd_macros`, `tensors-functions`, `visco_kinds`, `vp_laws_kinds`, `math_kernels_target`) | build con ifx | Corregido: `/fp:precise` en `src/` y `tests/` con ifx en Windows, salvo `Profile` (D-037). Desde D-042 la biblioteca usa el modelo rápido y los tests siguen con `/fp:precise` (estos tests no cambiaron); otras comprobaciones pasan a 2–32 ulp (D-042) | D-037, D-042 |
| ifx 2026.2 (Windows, CPU): `error #8212` «Omitted component is not initialized [VISCO]» al llamar `mat_state_*()` sin argumentos. Causa: el componente de tamaño cero `h(6, 0)` de los materiales con 0 ramas (Kelvin); ifx no lo trata como inicializado, con ninguna forma de inicializarlo | `tests/materials` (`mat_closest_point*`, `mat_solver_checks`); gfortran lo acepta | Corregido: `h(6, max(1, visco_branches))` en `muscle_material.fypp` | D-037 |
| ifx 2026.2: `write(formatted)` no encuentra la rutina UDTIO si el tipo y la interfaz llegan por un mismo módulo que los reexporta (`muscle_tensors`): `iostat 61` (escribe como tipo intrínseco) o texto vacío. Reproducido con una fachada de 4 líneas. Funciona con el tipo importado de su módulo de familia, con la interfaz importada de `muscle_tensor_io`, o con `generic :: write(formatted)` ligado al tipo. Los 5 `*_kinds` sólo corren desde D-036 (antes no compilaban); en ellos además fallaba `0 .ddot. (-1) = -0`: ifx usa `/assume:nominus0` por defecto y `SIGN(1.0, -0.0)` da `+1.0` (la suma sí daba `-0`) | `tensors-types-{2o,2osym,4o,4o2sym,4o3sym,iden}_kinds` (6 tests); afectaba a quien hiciera `use muscle_tensors` y escribiera con `DT` usando ifx | Corregido: `write(formatted)` ligado a cada tipo (D-038) y `/assume:minus0` en la rama Windows de ifx de `src/` y `tests/` (D-037) | D-037, D-038 |
| Windows: la pila del hilo principal (1 MB) desborda en `mat_batch_openmp` (`STATUS_STACK_OVERFLOW`); con 8 MB pasa | `materials-mat_batch_openmp`, gfortran y ifx | Corregido: los tests reservan 16 MB (`/F16777216`, `-Wl,--stack,16777216`) | D-037 |
| ifx: `omp_target_alloc`/`memcpy`/`free` sólo funcionan en un programa compilado con offload; un build de host los enlaza con `omptarget.lib` y el programa se cae en `omptarget.dll` (también con `OMP_TARGET_OFFLOAD=DISABLED`) | 8 tests `mat_batch_dev_*` | Corregido: con ifx (Windows y Linux) se compilan con `-Qopenmp-targets:spir64` (`-fopenmp-targets=spir64` en Linux) y `ctest` los ejecuta con `OMP_TARGET_OFFLOAD=DISABLED` (el host es el único dispositivo), así que corren en CPU | D-037, D-039 |
| gfortran 16.2 (Windows, CPU): los tests `target` (región `target` con libgomp en el host) se cuelgan de forma intermitente: `gpu-family_target_hard_sp` (1 de 18 corridas seguidas agotó 60 s; las demás tardan ~0.3 s), y desde 2026-11-19 también `hard_dp`, `os_dp` y `math-math_kernels_target`, en Debug y `Release` (`hard_*`: 2 de 30 corridas del par en `HEAD` y 2 de 30 con D-043; `math_kernels_target`: pasó 2 veces y se colgó la tercera). Los ejecutables no cambian entre corridas: no depende del código; con ifx y en Linux no se observó. **Causa (investigada 2026-11-19)**: interbloqueo al salir del proceso, en la cadena de herramientas (libgomp 16.2 + winpthreads 12 de conda-forge), no en MUSCLE. Se reproduce con un programa de 12 líneas (10 regiones `!$omp target teams distribute parallel do` sobre n=64, `gfortran -O3 -fopenmp -foffload=disable`): ~30 de 100 corridas se cuelgan **después** de imprimir el resultado, y un `!$omp parallel do` con el mismo bucle no se cuelga (0/100). Con más hilos es más probable (32 núcleos; `OMP_NUM_THREADS` 4: 0/400, 8: 1/400, 16: 0/100, 20: 5/100, 24: 9/100, 31: 29/100, 64: 42/100); `OMP_NUM_TEAMS=1`, `OMP_TEAMS_THREAD_LIMIT=1` y `OMP_WAIT_POLICY` no cambian nada. En el cuelgue queda un solo hilo (`Wait/UserRequest`, sin CPU); su pila, tomada con `cdb -pv` (WinDbg, modo no invasivo: gdb no puede adjuntarse porque inyecta un hilo y el proceso ya está saliendo): `common_exit` → `ExitProcess` → `LdrShutdownProcess` → `libgomp!__DllMainCRTStartup` (`DLL_PROCESS_DETACH`) → `execute_onexit_table` → `libgomp!__do_global_dtors` → `libwinpthread!pthread_key_delete` → `pthread_mutex_lock` → espera para siempre. `ExitProcess` mata primero los demás hilos y luego corre los destructores de las DLL; el mutex de winpthreads queda tomado por un hilo muerto. Con una pausa de 1 s antes de salir no se cuelga (0/100): hipótesis (no vista, el hilo dueño ya no existe) = en el respaldo en el host, `target teams` deja hilos terminando (y corriendo sus destructores TLS bajo ese mutex) cuando el hilo principal llega a `exit`; los de un `parallel` quedan dormidos en el pool | `tests/gpu`, `tests/math` (tests `target`) | **Documentado, sin corregir** (decisión del usuario, 2026-11-19): es de la cadena de herramientas; si se cuelga, repetir el test. Mitigación posible, no aplicada: `OMP_NUM_THREADS=4` en los tests `target` de gfortran en Windows sin GPU. Candidato a reporte a mingw-w64 (winpthreads) o GCC (libgomp), a decidir por el usuario | — |
| gfortran 16.2 (Windows): error interno del compilador (`Segmentation fault`, pase GIMPLE `threadfull`) al compilar `test_muscle_material_smoke` en Debug (`-O3`), una vez, con otros builds corriendo en paralelo; al recompilar pasó | `tests/materials` | Vigilar; si se repite, reportar con `-freport-bug` | — |
| ifx 2026.2 (Linux, `Release`): aviso `#5462` «Global name too long, shortened» en `muscle_mat_pz_kelvin_cha_hill48_swift`. Causa: el símbolo de una variable local que genera el compilador en `mat_batchdev_<nombre>_dp`, `muscle_mat_<nombre>_mp_mat_batchdev_<nombre>_dp_$blk.var$1072`, lleva el nombre del material dos veces (módulo y específico privado) y con 26 caracteres mide 97, uno más que el máximo (ifx quita el primer carácter). No sale en Debug, en Windows ni con gfortran | `src/materials` (nombres de material de más de 25 caracteres) | Documentado e inocuo (símbolo local, no la API; decisión del usuario 2026-11-19: dejarlo); aviso en `pages/integration.md`. Si molestara: quitar el nombre del material de los específicos privados (cambia D-011.4) | — |
| ifx 2026.2 con la biblioteca en modo rápido (sin `/fp:precise`): dos caídas intermitentes del compilador, que no se repitieron al recompilar: «compilation aborted (code 3)» en `muscle_yield_yld2004` (Linux, 32 trabajos en paralelo) y error interno `C0000005` en `muscle_mat_cut_maxwell_vm_bilinear` (Windows, `Release`). Con `/fp:precise` no se observó ninguna | `src/` (ifx, banderas rápidas) | **Vigilar**: desde D-042 la biblioteca usa el modelo rápido; en las corridas de D-042 (8 builds) no se repitieron. Si reaparece, repetir la compilación | D-041, D-042 |
| Los bancos no tienen la pila de 16 MB de los tests: `bench_kin_kernel` desborda la pila en Windows (`forrtl: severe (170)`) tras los casos sp | `benchmarks/CMakeLists.txt` (`MUSCLE_LINK_STACK` sólo en los tests) | **Pendiente**: añadir `MUSCLE_LINK_STACK` a los bancos | D-041 |
| El caso JC completo + Voce modificado de `closest_point_time` tiene `check_sp=.false.`: en sp sólo se acotan las iteraciones de los pasos que convergen (D-041); nada comprueba el estado ni los valores de sp, ni que sp converja (con el modelo rápido de ifx, sp deja de converger en un paso sin que nada lo detecte) | `tests/materials/test_muscle_material_closest_point_time.fypp` | Corregido: con el piso de tasa (D-043) sp converge y se comprueba | D-041, D-043 |
| Caso sensible JC completo + Voce modificado (`tm_jcfull_vm_vocemod`): su camino del Newton cambia con los últimos bits de las operaciones (con el modelo preciso, una perturbación relativa de 1e-11 de la deformación hace que dejen de converger 2 pasos en dp; 2 pasos no convergen nunca). Con el modelo rápido de ifx su referencia difiere en 1.7e-9 (categoría 1e-9) y su tangente numérica en 5.6e2 | `tests/materials` (`solver_case%sensitive`); Newton de `closest_point_time` con JC completo | Corregido: la causa era el modelo (flujo negativo bajo 1e-31, 0 en tasa nula); con el piso `epmin` converge en 6–8 iteraciones y se retiró `sensitive` | D-042, D-043 |
| ifx con fallback en host: NaN en las regiones target de `kin_kinds`. Causa: en el código host de un procedimiento `declare target`, un constructor `[x, y]` con x de tipo array conserva sólo el último elemento de x. El código de dispositivo está bien | tests (`pack_state`); `src/` no usa el patrón | Rodeado en el test; lint `arrayctor`; reproductor en el scratchpad (`r6.f90`) | D-019.4e |
| dp emulado en la UHD 770 deja `__isinfd`/`__isnand` sin resolver. Causa: `exponent()` → `for_exponent8_v` | spectral, cpb06/cpb06ex2/yld2004 | Corregido: `fp_exponent` con los mismos bits; lint `exponent` | D-019.5b |
| IGC (JIT de la UHD 770) compila mal closest_point y closest_point_2 con más de un work item: tensiones erróneas, status 104 o DEVICE_LOST. Bien con 1 work item o `-ze-opt-disable` (salvo la rampa); cutting_plane bien; valgrind limpio | `mat_batch*` en SPIR64 | Causa principal: el temporal de `COND_INF` (D-030); sus kernels de material sp vuelven a pasar en GPU. `mat_batch`/`mat_batch_openmp` siguen excluidos (sin volver a probar) | D-019.5c, D-030 |
| IGC/UHD 770 con kernels de material: closest_point, cp2 y time en sp pierden el dispositivo o dan otros iterados; en dp (emulado) hay fallos no reproducibles entre corridas. Sólo pasan en GPU 4 materiales sp (pz_cha_yld2004, pz_cpb06ex2, cut_yld2004, tm_jc_cpb06) | `materials-mat_kernels_target_*` en SPIR64 | sp: tras D-030/D-031 pasan todos salvo `tm_rk_vm` y `tm_jcfull_vm_vocemod` (otros iterados también sin fusión de EUs), que con el código de D-048 pasan en Linux y en Windows (2026-11-21; sin medir si ya pasaban sin su rodeo, tras D-043) y vuelven a correr; dp: siguen excluidos (`tools/gpu_skip_tests.sh`) | D-019.5f, D-030 |
| `scale()` en ifx/spir64 llama a `for_scale4_v/8_v` (dp emulado también en kernels sp). Reemplazarlo hacía perder el dispositivo en 3 materiales sp: la causa era el temporal de `COND_INF` (fila siguiente), no el reemplazo | cpb06, cpb06ex2, yld2004, `spectral_decomposition` | Corregido: `fp_scale` (`muscle_math_fp`); lint `scale`; la auditoría SPIR64 ya no permite `for_scale*_v` | D-019.6, D-031 |
| IGC (UHD 770): un temporal de tamaño variable en código de dispositivo no se puede acotar. `maxval(sum(abs(A), dim=2))` de `COND_INF` sobre un `A(:,:)` de forma asumida genera un `alloca` dinámico (el único de la imagen). Si el kernel no tiene stack calls, IGC reserva 4 KB de memoria privada por work item (`private_size 32768` por hilo SIMD8) y el kernel la desborda: DEVICE_LOST, también con 1 work item, `-ze-opt-disable` y la fusión de EUs apagada (`IGC_StackOverflowDetection=1`: «Stack overflow detected!»; con `IGC_ForcePerThreadPrivateMemorySize=65536` pasa). Con stack calls (las de la emulación FP64 de `scale()`) reserva 1,4 MB y alcanza: es la causa del device lost al reemplazar `scale()` (D-019.6, comprobado en sus 3 kernels). Con D-030 y D-031 juntas vuelven a pasar además 8 kernels sp excluidos (closest point, closest point 2, closest point time y Perzyna; D-019.5c/5f); no se midió cuántos se deben sólo a D-030 | `COND_INF`, `M33INV`/`M44INV`, LU con `KAPPA` (inversos y return mapping) | Corregido: `NORM_INF` (bucles explícitos, sin temporales); regla en `gpu-static-dispatch.md` | D-030 |
| IGC (UHD 770, Gen12LP): con emulación FP64, llamadas a las funciones de emulación dentro de una rama divergente y arreglos locales indexados con tablas `PARAMETER` (-O1 o más), si sólo uno de dos hilos SIMD8 de una pareja de EU fusionadas ejecuta la rama, los 8 carriles del otro devuelven un valor erróneo. Mismo código máquina correcto con `IGC_DisableEuFusion=1`; también con IGC 2.42.6 / NEO 26.40. Rompía los candidatos de Jacobi de `PLAN_ESPECTRAL.md` en GPU | kernels sp con `scale()` y kernels dp (emulados) | Rodeado en sp: `fp_scale` saca las llamadas a la emulación; rodeos `IGC_EmulationFunctionControl=1` / `IGC_DisableEuFusion=1`; reproductores e informe en `docs/fypp/igc_mre/` | D-031 |
| ifx/IGC Release (-O2) en la UHD 770: las imágenes de `pz_cpb06ex2_bilinear` y `tm_jc_cpb06_ludwik` cargan sin los kernels de lote (segfault); sólo dan resultados correctos en Debug (-O0). Sin dp emulado: «Double type is not supported» (sp y dp comparten imagen); con `-ze-opt-disable`, device lost | benchmark GPU | Parcial (D-035): con la LU desenrollada aparecen los kernels sp de `tm_jc_cpb06_ludwik` (corre en Release); siguen faltando sus kernels dp y los de `pz_cpb06ex2_bilinear` y `pz_cha_yld2004_ludwik` (otra causa, sin investigar). Validar en NVIDIA/AMD | D-019.8, D-035 |
| `pz_cha_yld2004_ludwik` no corre en la GPU ni en Release ni en Debug fuera del test target | benchmark GPU | Documentado | D-019.8 |
| Un build de CPU con ifx del benchmark no enlaza (`omp_target_alloc` sin definir con `-fiopenmp` y sin target de offload) | benchmark | Documentado; las filas de CPU usan gfortran | D-019.8 |
| La regla de lint `quote` nunca corría (error de sintaxis en perl) | `tools/check.sh` | Corregido | D-019.4 |
| IGC (UHD 770): SIGSEGV intermitente dentro de libigc al compilar en el JIT `mat_kernels_target_pz_cha_yld2004_ludwik_sp` durante `check.sh` (con carga de CPU); pasa 3/3 al repetirlo solo. Código generado idéntico (Fase 9, E3) | tests target SPIR64 | Documentado: repetir el test; no es un cambio de código | D-020 |
| Driver Intel de Windows (UHD 770, 32.0.101.5972 y 32.0.101.7092): sin emulación de FP64. Una imagen con cualquier operación dp no carga (`Double type is not supported on this platform`, al crear el módulo de Level Zero), también la de un test sp de materiales: cada módulo de material trae sus kernels por lotes sp y dp en la misma imagen. Probado sin efecto (2026-11-21): las variables de Linux (`OverrideDefaultFP64Settings=1`, `IGC_EnableDPEmulation=1`), `NEO_FP64_EMULATION=1` y `NEOReadDebugKeys=1`, las claves de registro `HKLM\SOFTWARE\Intel\IGFX\{OCL,IGC,L0}` y la opción de
compilación del módulo `-ze-fp64-gen-emu` (que el IGC de Windows sí conoce, pasada con `LIBOMPTARGET_LEVEL_ZERO_COMPILATION_OPTIONS`). El IGC publicado para Windows no contiene la opción `EnableDPEmulation` (la de Linux); el runtime de Level Zero lee `...\IGFX\L0` y conoce `NEO_FP64_EMULATION`, pero no la activa en esta GPU (¿sólo Xe-HPG?, sin confirmar) | GPU Intel en Windows | Externo: build sin dp (`-DMUSCLE_REAL_KINDS=sp`) para la GPU; dp en GPU sólo en Linux | D-047 |
| Windows: el `libiomp5md.dll` de conda (paquete `llvm-openmp`, llega con MKL al entorno `gfortran`) tapa al de oneAPI si el entorno se activa después de `setvars.bat`: sirve para OpenMP en CPU, pero un programa de ifx con regiones `target` no arranca (`0xC0000139`, punto de entrada no encontrado; también en el build de CPU y con `OMP_TARGET_OFFLOAD=DISABLED`) y bajo ctest parece colgado hasta el límite de tiempo: con ifx en CPU fallaron los 58 tests `target` que no son `batch_dev` (2026-11-21). Es del entorno `gfortran` de esta máquina (`llvm-openmp` llegó con MKL el 2026-11-14); el entorno `muscle` de `environment.yml` (OpenBLAS) no lo trae (deducido de sus dependencias, no instalado). Sin explicar: con el mismo orden ifx pasó 133/133 el 2026-11-20 | tests target con ifx en Windows (CPU y GPU) | Externo: `compiler\bin` de oneAPI primero en el `PATH` (`tools/limit_run.ps1` no lo arregla: lo hace quien llama) | D-047 |
| GPU de Windows (UHD 770, driver 32.0.101.7092, ifx 2026.2, build sin dp): 7 de 31 tests `target` difieren de la CPU, y los 7 pasan en Linux (WSL, misma GPU, NEO 26.05, mismo código, 48/48, 2026-11-21). Tres patrones: `cp_kelvin_hill48_swift` y `cp2_maxwell_vm_bilinear`, la historia viscoelástica (`strain_e`, `h(i)`) sale del orden de la tensión en todos los puntos (tensión y tangente ≤ 14 eps); `cut_yld2004_ludwik`, NaN y otros status; `tm_mrk_vm`, `tm_nnl_vm`, `tm_va_vm`, `tm_jc_cpb06_ludwik`, la tangente de los puntos plásticos difiere 1e-5–1e-3 relativo (tensión ≤ 4 eps). Investigado (2026-11-21): no son las banderas (en un build con `MUSCLE_GPU` ifx compila igual en los dos sistemas: biblioteca
`/Od`/`-O0`, tests `-O3`, sin modelo de punto flotante) ni el código (el mismo pasa en Linux): apunta al IGC del driver de
Windows. Determinista (el mismo error, al dígito, entre corridas y sin la caché de kernels, `NEO_CACHE_PERSISTENT=0`). No sirve
para acotarlo: con cualquier opción real de IGC (`LIBOMPTARGET_LEVEL_ZERO_COMPILATION_OPTIONS` = `-ze-opt-disable`,
`-ze-opt-disableEUFusion`, `-ze-opt-disable-compaction`) estos kernels, y también uno que pasa (`cp_dp_ludwik`), pierden el
dispositivo (`DEVICE_LOST`); la opción neutra `-cl-std=CL2.0` no cambia nada (la variable suma opciones). Más memoria privada
(`IGC_ForcePerThreadPrivateMemorySize=16384`) y `IGC_StackOverflowDetection=1` no cambian el resultado (sin confirmar que el IGC
de Windows lea esas variables). El IGC de Windows no tiene las palancas de Linux (`IGC_DisableEuFusion`,
`IGC_ShaderDumpEnable`). **Causa (2026-11-21, imprimiendo la salida de un punto de `cp_kelvin_hill48_swift`):** en el
dispositivo, el argumento `strain` (`intent(in)`) de `mat_update` sale del return mapping con la tensión que éste escribió en
`stress`: `strain_e = strain - strain_p` daba la tensión (también restando componente a componente; un valor centinela
comprobó que se ejecutaba la copia editada). Interpretación (sin reproductor mínimo): el IGC de Windows ubica el `strain` y el
`stress` de quien llama en la misma memoria privada. Lo que el kernel lee de `strain` después de escribir la tensión sale mal:
la historia viscoelástica, la tangente de los `tm_*` y el plano de corte (NaN, status) | `mat_kernels_target_*_sp` en Windows
| Corregido con un rodeo: `mat_update` trabaja sobre una copia local de `strain` (`muscle_material_kernel.fypp`); 31/31 en la
GPU de Windows | D-048 |

## Registro de soluciones provisionales de compiladores (vigentes)

Cada fila es código o configuración que existe sólo por un defecto de un compilador o driver. «Verifica»: qué
detecta si se rompe o si se puede quitar. Revisar esta tabla al cambiar de versión de compilador.

| Solución provisional | Compilador / versión | Dónde | Verifica | Retirar cuando |
|---|---|---|---|---|
| Sumas explícitas en lugar de `matmul` vector×matriz | ifx 2026.2 (código de dispositivo) | `muscle_tensor_ops_contraction_double.fypp` | build SPIR64 y tests target de `check.sh` | ifx compile bien `matmul` en `declare target`; el usuario decidió mantener las sumas en todos los builds (D-019.4c), así que retirarla es opcional |
| Copia local `strain_in` de `strain` en `mat_update` (el return mapping no recibe el argumento) | IGC del driver de Windows 32.0.101.7092 (UHD 770): `strain` sale con la tensión | `muscle_material_kernel.fypp` (variante de deformación pequeña) | tests target de materiales en la GPU de Windows | el IGC de Windows no comparta la memoria de `strain` y `stress` (correr los target sin la copia) |
| `--line-numbering-mode=nocontlines` de fypp | ifx (rechaza `#line` dentro de un literal continuado) | `cmake/MuscleFypp.cmake` | build SPIR64 | ifx acepte marcadores entre líneas continuadas |
| `fp_exponent` en lugar de `exponent()`; regla de lint `exponent` | ifx con dp emulado en GPU Intel (`for_exponent8_v` no existe en el dispositivo) | `muscle_math_spectral_derivs`, cpb06, cpb06ex2, yld2004 | lint y auditoría SPIR64 | el runtime de dispositivo de ifx tenga la rutina, o GPU con FP64 nativo |
| Regla de lint `arrayctor` (sin constructor `[x, y]` encabezado por un array en código host de `declare target`) | ifx 2026.2 (conserva sólo el último elemento) | todo el código | lint de `check.sh` | ifx corrija el constructor (reproductor `r6.f90`) |
| Módulo `muscle_tensor_minus` (fusiona sólo `operator(-)`) | ifx 2026.2 (error #6549: se pierde el `-` unario) | `src/tensors/muscle_tensor_minus.fypp`, fachada | build con ifx | ifx resuelva la fusión de genéricos |
| Componente `h(6, max(1, visco_branches))` del estado visco, con una columna sin uso cuando hay 0 ramas | ifx 2026.2 (`#8212`: no trata como inicializado un componente de tamaño cero) | `src/materials/muscle_material.fypp` | build con ifx de `mat_closest_point*` y `mat_solver_checks` | ifx inicialice por defecto los componentes de tamaño cero |
| `/fp:precise` (`-fp-model=precise` en Linux) de ifx en `tests/` (salvo `Profile`); la biblioteca, modelo rápido y, en Debug, especulación segura (D-042; sin `/Qftz-` desde D-043) | ifx (modelo de punto flotante rápido por defecto) | `cmake/MuscleCompiler.cmake` | build con ifx en Windows y en Linux | nunca: es la semántica que los tests suponen |
| Los tests `mat_batch_dev_*` se compilan con `-Qopenmp-targets:spir64` y se ejecutan con `OMP_TARGET_OFFLOAD=DISABLED` con ifx, en Windows y en Linux, sin `MUSCLE_GPU` (opción `HOST_OFFLOAD` de `muscle_target_test_options`) | ifx (`omp_target_*` sólo funciona en un programa compilado con offload) | `tests/CMakeLists.txt`, `cmake/MuscleCompiler.cmake` | `mat_batch_dev_*` en CPU con ifx | ifx sirva `omp_target_*` en un programa de host |
| Pila de 16 MB al enlazar los tests en Windows | Windows (hilo principal de 1 MB) | `cmake/MuscleCompiler.cmake` (`MUSCLE_LINK_STACK`) | `mat_batch_openmp` | permanente (es del sistema operativo) |
| `/Ob0` en lugar de `-fno-inline` para los tests `target` con ifx en Windows: `-fno-inline` se ignoraba (aviso `#10006`) y los tests suponían una sola copia de cada procedimiento | ifx en Windows (la opción de no inlining es `/Ob0`) | `cmake/MuscleCompiler.cmake` (`MUSCLE_FLAG_NO_INLINE`) | los 66 tests `target` con ifx | nunca: es el nombre de la opción en ese compilador |
| Trampas de punto flotante (`/fpe:0`, `-fpe0` en Linux; desde D-042 sólo en Debug, donde la biblioteca lleva especulación segura) en los 5 tests de `math` que en gfortran llevan `-ffpe-trap=zero,overflow,invalid`. Con ifx `-fpe0` también aplana los denormales (según la documentación de Intel; no medido): por eso sólo en esos 5 tests, no en todos | ifx (sin las trampas esos tests comprobaban menos que con gfortran) | `cmake/MuscleCompiler.cmake` (`muscle_fpe_traps`) | `math-math_spectral_kinds`, `eigen_eigenvals`, `spectral_cases`, `eigen_projections`, `matrix_functions` | nunca |
| `generic :: write(formatted)` ligado a cada tipo, con `muscle_tensor_io` sin dependencia de los tipos (la fachada no reexporta la interfaz) | ifx 2026.2 (no encuentra una interfaz `write(formatted)` reexportada junto con su tipo: `iostat 61` o salida vacía) | `src/tensors/types/*.fypp`, macro `TENSOR_UDTIO_TBP`, `muscle_tensor_io.fypp` | build con ifx de `tensors-types-*_kinds` e `iden_kinds` | ifx resuelva la interfaz reexportada; el diseño enlazado es válido en ambos compiladores y puede quedarse |
| `/assume:minus0` (`-assume minus0` en Linux) de ifx en `src/` y `tests/` | ifx (`SIGN(1.0, -0.0)` da `+1.0` por defecto; estándar y gfortran dan `-1.0`) | `cmake/MuscleCompiler.cmake` | `tensors-types-2o_kinds` (`0 .ddot. (-1) = -0`) | nunca: es la semántica estándar |
| Módulo `muscle_tensor_aliases` (renombres históricos `ten_3D2O => ten_3D2O_dp`, ...; la fachada no renombra) | ifx 2026.2 (error interno #5623 en `call p%init(...)` si un módulo importa un tipo con dos nombres, D-036) | `src/tensors/muscle_tensor_aliases.fypp`, fachada | build con ifx de los 5 tests `tensors-types-*_kinds` (compilan sin error interno) | ifx acepte dos nombres de un tipo en un mismo módulo |
| Tipos derivados como locales de un procedimiento `declare target`, nunca en `private(...)` | ifx (error interno) | `BATCH_LOOP` y kernels de lote | build SPIR64 | ifx compile tipos derivados en `private` |
| Asignar componentes (`t%vals = x(1:6)`) en lugar de un constructor con una sección de argumento de forma asumida | IGC (recursión infinita del JIT) | kernels | tests target en la GPU Intel (con límite de memoria) | IGC corregido |
| Exclusiones de `tools/gpu_skip_tests.sh` (4 patrones) | IGC del driver 26.05 de la UHD 770 | tests target de materiales | `check.sh` lista lo excluido | driver nuevo de Intel, o validación en NVIDIA/AMD |
| NaN con `transfer` del patrón de bits (sin módulos IEEE en kernels) | gfortran nvptx/amdgcn (llamadas IEEE sin implementación en el dispositivo) | `QUIET_NAN` | builds NVPTX/AMDGCN y auditoría | probablemente permanente (regla GPU) |
| `MAT33_TO_VEC9` en lugar de `reshape` (la otra dirección es `to_matrix`, D-028) | gfortran (llama a `_gfortran_reshape_*` en el dispositivo) | `include/muscle_common.fypp` y usuarios | auditoría de imágenes | probablemente permanente |
| `MUSCLE_UDTIO=AUTO` desactiva la salida formateada | gfortran < 16.3 (locales estáticos con UDTIO y OpenMP) | `src/CMakeLists.txt` | `check.sh` compila el camino `MUSCLE_UDTIO=ON` | se retira solo: con gfortran ≥ 16.3 `AUTO` la activa |
| Temporales tipados en lugar de `associate` con resultados de función de tipo derivado | gfortran | código de tensores | compilación | gfortran corregido |
| Ejecutar programas GPU con `systemd-run` (límite de memoria) y `timeout`; nunca `ulimit -s unlimited` | IGC (un bug llegó a usar 26 GB) | `check.sh`, reglas de tests | — | práctica permanente |

Corregido: el error interno de ifx en `call p%init(...)` con dos nombres del mismo tipo (D-019.4f) se evita con
`muscle_tensor_aliases` (D-036).

## Tests

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Rotación de Jacobi de las autoproyecciones: `theta = (m(q,q) - m(p,p))/(2 m(p,q))` desborda en sp cuando `m(p,q)` queda cerca de `tiny` con autovalores ~1e3 (excepción de punto flotante; Inf sin trampa); afectaba a la cinemática finita en sp con estiramientos grandes | `src/math/muscle_math_spectral_derivs.fypp` (`jacobi_eigen`, antes en `muscle_kin_finite`) | Corregido (2026-11-11): para theta > 1e15, `t = m(p,q)/(m(q,q) - m(p,p))` sin formar theta; test `math-math_spectral_projections` (escala 1e3 en sp) | D-028 (`jacobi_eigen` retirado en D-032; la misma tangente sigue en `spectral_decomposition`) |
| `-Wcompare-reals` (de `-Wextra`) avisa en 26 comparaciones exactas deliberadas de `src` (cero exacto antes de dividir en yield, interruptores `prm%m /= 0` en JC completo y VA, pivote nulo en `gauss_solve`/`gauss_seidel`, `.approx.` con valores no finitos dentro de `TENSOR_COMMON_PROCEDURES`) | Flags Debug de `cmake/MuscleCompiler.cmake` | **Pendiente** (usuario, 2026-11-10): se silencia con `-Wno-compare-reals` para que el build quede sin avisos, pero no debe quedar siempre; retirar el flag cuando esas comparaciones se reescriban o se marquen una a una | REANALISIS_MANTENIBILIDAD 2.1 |
| Goldens bit a bit fallaban con otros flags (`Profile_unittests`: `kin_golden`, `math_golden_derivatives`) | `tests/common/muscle_test_golden.f90` | Corregido: tolerancia por categoría y modo estricto opcional | D-021 |
| El target `coverage_report` fallaba con lcov 2.x (marcadores `#line`: «mismatched end line») y medía también los tests | `CMakeLists.txt` | Corregido: `--ignore-errors inconsistent` con lcov ≥ 2 y extracción de `src/` e `include/` | ESTADO_OPTIMO 1.4 |
| Con FMA (`-march=native -O3`): registros cuyo valor exacto es 0 dan un error relativo enorme (escala ~0): `elas_golden` (`volumetric \| stress_dev`), `yield_golden` (yld2004 `s31`) | goldens de elasticidad y yield | Resuelto (plan R): el test se reemplazó por tests unitarios con tolerancia por kind | D-021 |
| Con FMA: `math_golden_spectral_derivs` compara autovectores de autovalores repetidos (no únicos; error 1.57 con una base válida distinta) | `tests/math/test_math_golden_spectral_derivs.f90` | Resuelto (plan R): el test se reemplazó por tests unitarios de invariantes (V Λ Vᵀ, ortogonalidad); `math_reference` guarda 2 empates exactos que con FMA podrían elegir otra base válida (pendiente de medir con `-march=native`) | D-021 |
| Con FMA: `kin_golden` `nearly singular general \| pull_back_tangent` 2.7e-8 (mal condicionado) | `tests/kinematics/test_kin_golden.fypp` | Resuelto (plan R): el test se reemplazó por tests unitarios con tolerancia por kind | D-021 |
| Con FMA: tangente numérica de `closest_point_time` (jcfull + vocemod) cambia 1e3 relativo y aparece NaN: diferencias finitas en un punto no suave (piso de JC) | `tests/materials/test_mat_closest_point_time.fypp` | Corregido: con el piso `epmin` la tangente numérica es finita y ≤ 4e-9 de la referencia (D-043) | D-021, D-043 |
| Con FMA: diferencias finitas de math 1.4e-9 (> 1e-9) y `math_golden_operations` (derivadas por diferencias en la categoría álgebra) 2.6e-10 | goldens de math | Resuelto (plan R): el test se reemplazó por tests unitarios con tolerancia por kind | D-021 |
| Los tests de compilación fallida (parámetros privados) pasaban con cualquier error porque la ruta del archivo contiene "private" | `tests/plasticity/yield_criteria`, `tests/elasticity` | Corregido: se exige "PRIVATE component" y fallan con "Cannot open module/Fatal Error" (sólo gfortran) | — |
| `check_approx` (escalares y arrays reales) aceptaba cualquier valor si un operando era Inf: la escala `max(1, |got|, |expected|)` salía infinita (encontrado por el chequeo de mutaciones del plan R: `C*dx` → `C/dx` en JC devolvía Inf y pasaba) | `tests/common/muscle_test_utils.fypp` | Corregido (R3): un componente no finito sólo pasa si es igual al esperado; la escala de un array usa sólo los componentes finitos. Ningún test escondía un Inf | ESTADO_OPTIMO R3 |
| `eigenvals`: la rama `abs(1 - d) <= EPS` (raíz media ≈ media de las otras) da autovalores ~1e-12 relativo lejos de la forma cerrada en `diag(5, 3 + 1e-11, 1)`; sin la rama la forma cerrada acierta (y sigue exacta en los casos exactos) | `src/math/muscle_math_eigen.fypp:101` (antes `muscle_math_operations`) | **Marcado, sin corregir** (desde D-032 no lo usan los criterios anisótropos; tras F2, ningún kernel de la librería) (usuario, 2026-11-09: queda pendiente; corregirlo cambia resultados). Mientras tanto la referencia mínima `math_reference` la protege | ESTADO_OPTIMO R2 |
| El mismo hueco en `is_approx`/`.approx.` de los tensores (escala con un Inf ⇒ todo pasa) y en `.approx.` de las identidades escaladas | `include/muscle_tensor_common.fypp`, `src/tensors/types/muscle_tensor_iden.fypp` | Corregido (2026-11-09, decisión del usuario): un componente Inf/NaN sólo coincide con el mismo Inf y la escala usa sólo los componentes finitos (`abs(x) <= huge(x)`, apto para kernels); tests en `test_approx` (todas las familias) y `test_approx_2O` (identidades), que fallan con el código anterior | ESTADO_OPTIMO R3 |

## Compiladores y herramientas (Externo)

| Problema | Dónde | Rodeo | Ref. |
|---|---|---|---|
| gfortran: `associate` con resultados de función de tipo derivado falla | código de tensores | Temporales tipados | — |
| gfortran 15: reasignar un `class(...), allocatable` ya asignado con otro tipo dinámico no realoca y desborda el heap | p.ej. `law%hard_law = otra_ley` | Desasignar antes de reasignar | D-017.7 |
| gfortran < 16.3: salida formateada (UDTIO) y tests OpenMP fallan juntos. Mecanismo (2026-11-02): si un ámbito ve `write(formatted)` de un tipo, gfortran da almacenamiento estático a las variables locales de ese tipo, incluso en `pure`. Como `muscle_tensors` lo reexporta, los locales de los kernels se comparten entre hilos (carrera). Con `MUSCLE_UDTIO=ON` forzado fallan `kin_kinds`, `math_kernels_target` y `mat_batch`; viene de antes de la Fase 7 | `MUSCLE_UDTIO=ON` | Mitigado: `AUTO` lo desactiva con GNU < 16.3 (`src/CMakeLists.txt`) y `check.sh` sólo corre ahí los tests de tensores. La versión 16.3 como corrección no está verificada | D-011 |
| ifx 2026.2: se pierde el `-` unario al fusionar `operator(-)` de familias y de operadores en un módulo | `muscle_tensors` | Módulo `muscle_tensor_minus` que fusiona sólo `operator(-)` | gpu-static-dispatch.md |
| ifx: tipos derivados en `private(...)` de una región `target` dan error interno | kernels GPU | Declararlos como locales de un procedimiento `declare target` | gpu-static-dispatch.md |
| IGC (Intel GPU): constructor de estructura con sección de un argumento de forma asumida → recursión infinita del JIT (26 GB) | kernels GPU | Asignar el componente (`t%vals = x(1:6)`); ejecutar siempre con límite de memoria | gpu-static-dispatch.md, testing.md |
| Intel UHD 770: sin FP64 nativo; la emulación descarta kernels grandes sin error | GPU Intel | Probar dp en GPU con FP64 nativo (Fase 8) | PLAN Fase 8 |
| gfortran nvptx/amdgcn: los módulos IEEE generan llamadas sin implementación en el dispositivo | `QUIET_NAN` | NaN con `transfer` del patrón de bits | — |
| FORD no preprocesaba los `.fypp` (los módulos migrados desaparecían de la doc) | `my_project.md` | `preprocessor: fypp`, target `doc` | PLAN Fase 9 |
| fypp: los `#line` de una macro apuntan a la línea del `#:call` | macros de modelos | Patrón invertido: el autor escribe los cuerpos | D-017.6 |
