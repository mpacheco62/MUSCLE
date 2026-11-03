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
| `spectral_decomposition` sp pierde la normalización con cortes ~1e-11 (\|c\|² subnormal) | `muscle_math_spectral_derivs` | Corregido: escalado exacto por potencia de 2 | D-016.11 |
| Umbral absoluto 1e-10 de `eigenvals`: descarta cortantes chicos sobre una hidrostática | `eigenvals` (afecta a CPB06, CPB06ex2, Yld2004) | Documentado | D-017.9 |
| `eigenvals` divide por cero (1/0 = Inf, resultado correcto) con la raíz mayor repetida; aborta con trampa de división por cero | `eigenvals` | Pendiente (Fase 8, GPU) | PLAN Fase 5 |
| Un argumento no nulo pero diminuto conserva el paso mínimo histórico de las diferencias finitas | `include/muscle_fd.fypp` | Documentado | D-016.4 |
| Macros FD vs API host: los últimos bits cambian si el compilador integra la función (≈4e-11 con -O3) | `include/muscle_fd.fypp` | Documentado (comparar con tolerancia en la Fase 6) | D-016.9 |
| Golden de `mklsolve` depende de la LAPACK enlazada | `tests/reference/math` | Documentado (la cabecera anota la librería) | — |

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
| `flow_stress = 0` con epd ≤ 0 (tensión nula en el primer incremento) | RK, JCfull | Documentado | D-017.7 |
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
| jcfull en sp: paso 1 sin converger, después ep difiere ~6e-4 de dp y el tangente por DF da NaN (flujo 0 a tasa nula) | material `time` + JC completo | Documentado (sólo dp se comprueba) | D-018.12b, D-017.7 |
| Perzyna en sp no converge con incrementos de ~40 deformaciones de fluencia (4e-2), aun con el piso de redondeo | `perzyna` sp (material `kin_ul_pz_*`) | Documentado; el test usa incrementos ≤ 3e-2 | D-018.12a |
| cp2 con endurecimiento blando converge (status OK) a otra solución (5 % de diferencia) en el chequeo F→I | `cp2` + cinemática | Documentado (sin dγ ≥ 0) | D-018.10c |
| Incrementos de ~1000 deformaciones de fluencia (cpb06) y exponentes extremos de Perzyna no convergen en sp | `closest_point`, `perzyna` sp | Documentado; las pruebas físicas corren sólo en dp | D-018.18g |
| dp: cinemáticas logarítmicas (kin_ml_*, kin_sl_*) con F = I + 1e-6 dan σ vs ddsdde:ε con error relativo 0.67 (6e-5 con 1e-5). Causa: con brechas < ~3e-6 la forma cerrada toma las raíces como repetidas (P = I/3) y ln C pierde la parte desviadora; σ está mal, no el tangente | log material/espacial dp | Corregido (Jacobi en la franja [1e-12, 1e-2)) | D-018.19a |
| `kin_tl_tm_vm_jc` no converge (200 iteraciones) con 1e-4 de deformación, σ ≈ 21 < y0 = 200, en ambos kinds | material TL + time + JC (también sin cinemática) | Corregido: era plástico (flujo a tasa mínima 15.8 MPa) y el Newton se estancaba porque la derivada de tasa es 0 en el piso; el solver usa la derivada por la derecha | D-018.19b |
| `tm_nnl_vm` con deformación ×40 no converge (200 iteraciones en dp, 177 en sp) | `time` + NNL | Pendiente (Fase 9); excluido del chequeo de status del smoke test | D-018.19 |
| MRK: NaN desde el primer paso plástico (probable epd > epdmax de la ley, ver viscoplásticas) | `closest_point_time` | Documentado | D-018.2 |
| Tangente analítico vs numérico: RK ~8 %, NNL 3.5e-4 (pasos convergidos) | `closest_point_time` | Documentado | D-018.10d |
| Sin restricción dγ ≥ 0: incremento grande bilineal + von Mises converge a ep < 0 (−1.84e-2; cp da +2.00e-2); no converge con cpb06/cpb06ex2 (incremento grande) ni yld2004-2 (NaN) | `closest_point_2` | Documentado | D-018.10c |
| Chaboche n=0.5 + Swift + DP en compresión: 4/6 pasos plásticos sin converger; tangentes con diferencias 1.4–1.6e-5 (Perzyna n=0.5 Hill48, n=5 DP) | `closest_point_perzyna` | Pendiente (Fase 9); con D-018.12a persiste en el caso 6 (tangente DF vs analítico 0.24) | Fase 7 golden |

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

## Tests

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| Los tests de compilación fallida (parámetros privados) pasaban con cualquier error porque la ruta del archivo contiene "private" | `tests/plasticity/yield_criteria`, `tests/elasticity` | Corregido: se exige "PRIVATE component" y fallan con "Cannot open module/Fatal Error" (sólo gfortran) | — |

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
