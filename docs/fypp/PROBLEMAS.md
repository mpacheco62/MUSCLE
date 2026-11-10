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
| `eigenvals` divide por cero (1/0 = Inf, resultado correcto) con la raíz mayor repetida; aborta con trampa de división por cero | `eigenvals` | Corregido (Fase 8): con `n2 <= n1/huge` va directo al límite, sin dividir; mismos bits. Test con trampas FPE | D-019.6 |
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
| ifx 2026.2: error interno en `call p%init(...)` si el tipo llega por `muscle_tensors` con dos nombres (`ten_3D2Osym` y `ten_3D2Osym_dp`) | 5 tests `tensors-types-*_kinds` | Pendiente (reproductor para Intel) | D-019.4f |
| ifx con fallback en host: NaN en las regiones target de `kin_kinds`. Causa: en el código host de un procedimiento `declare target`, un constructor `[x, y]` con x de tipo array conserva sólo el último elemento de x. El código de dispositivo está bien | tests (`pack_state`); `src/` no usa el patrón | Rodeado en el test; lint `arrayctor`; reproductor en el scratchpad (`r6.f90`) | D-019.4e |
| dp emulado en la UHD 770 deja `__isinfd`/`__isnand` sin resolver. Causa: `exponent()` → `for_exponent8_v` | spectral, cpb06/cpb06ex2/yld2004 | Corregido: `fp_exponent` con los mismos bits; lint `exponent` | D-019.5b |
| IGC (JIT de la UHD 770) compila mal closest_point y closest_point_2 con más de un work item: tensiones erróneas, status 104 o DEVICE_LOST. Bien con 1 work item o `-ze-opt-disable` (salvo la rampa); cutting_plane bien; valgrind limpio | `mat_batch*` en SPIR64 | Excluido y documentado; reproductor `c3b.f90` en el scratchpad | D-019.5c |
| IGC/UHD 770 con kernels de material: closest_point, cp2 y time en sp pierden el dispositivo o dan otros iterados; en dp (emulado) hay fallos no reproducibles entre corridas. Sólo pasan en GPU 4 materiales sp (pz_cha_yld2004, pz_cpb06ex2, cut_yld2004, tm_jc_cpb06) | `materials-mat_kernels_target_*` en SPIR64 | Excluidos en `tools/gpu_skip_tests.sh`; validar en NVIDIA/AMD | D-019.5f |
| `scale()` en ifx/spir64 llama a `for_scale4_v/8_v`. Reemplazarlo por `fp_scale` (o por `x*2**n`) hace que IGC pierda el dispositivo en 3 materiales sp que hoy pasan | cpb06, cpb06ex2, yld2004, `spectral_decomposition` | Se mantiene `scale()`; permitido por nombre en la auditoría SPIR64 | D-019.6 |
| ifx/IGC Release (-O2) en la UHD 770: las imágenes de `pz_cpb06ex2_bilinear` y `tm_jc_cpb06_ludwik` cargan sin los kernels de lote (segfault); sólo dan resultados correctos en Debug (-O0). Sin dp emulado: «Double type is not supported» (sp y dp comparten imagen); con `-ze-opt-disable`, device lost | benchmark GPU | Documentado (validar en NVIDIA/AMD) | D-019.8 |
| `pz_cha_yld2004_ludwik` no corre en la GPU ni en Release ni en Debug fuera del test target | benchmark GPU | Documentado | D-019.8 |
| Un build de CPU con ifx del benchmark no enlaza (`omp_target_alloc` sin definir con `-fiopenmp` y sin target de offload) | benchmark | Documentado; las filas de CPU usan gfortran | D-019.8 |
| La regla de lint `quote` nunca corría (error de sintaxis en perl) | `tools/check.sh` | Corregido | D-019.4 |
| IGC (UHD 770): SIGSEGV intermitente dentro de libigc al compilar en el JIT `mat_kernels_target_pz_cha_yld2004_ludwik_sp` durante `check.sh` (con carga de CPU); pasa 3/3 al repetirlo solo. Código generado idéntico (Fase 9, E3) | tests target SPIR64 | Documentado: repetir el test; no es un cambio de código | D-020 |

## Registro de soluciones provisionales de compiladores (vigentes)

Cada fila es código o configuración que existe sólo por un defecto de un compilador o driver. «Verifica»: qué
detecta si se rompe o si se puede quitar. Revisar esta tabla al cambiar de versión de compilador.

| Solución provisional | Compilador / versión | Dónde | Verifica | Retirar cuando |
|---|---|---|---|---|
| Sumas explícitas en lugar de `matmul` vector×matriz | ifx 2026.2 (código de dispositivo) | `muscle_tensor_ops_contraction_double.fypp` | build SPIR64 y tests target de `check.sh` | ifx compile bien `matmul` en `declare target`; el usuario decidió mantener las sumas en todos los builds (D-019.4c), así que retirarla es opcional |
| `--line-numbering-mode=nocontlines` de fypp | ifx (rechaza `#line` dentro de un literal continuado) | `cmake/MuscleFypp.cmake` | build SPIR64 | ifx acepte marcadores entre líneas continuadas |
| `fp_exponent` en lugar de `exponent()`; regla de lint `exponent` | ifx con dp emulado en GPU Intel (`for_exponent8_v` no existe en el dispositivo) | `muscle_math_spectral_derivs`, cpb06, cpb06ex2, yld2004 | lint y auditoría SPIR64 | el runtime de dispositivo de ifx tenga la rutina, o GPU con FP64 nativo |
| Se mantiene `scale()` (permitido por nombre en la auditoría SPIR64) | ifx/IGC (reemplazarlo hace perder el dispositivo) | cpb06, cpb06ex2, yld2004, `spectral_decomposition` | auditoría SPIR64 y tests target | IGC no pierda el dispositivo con el reemplazo, o `for_scale*_v` exista en el dispositivo |
| Regla de lint `arrayctor` (sin constructor `[x, y]` encabezado por un array en código host de `declare target`) | ifx 2026.2 (conserva sólo el último elemento) | todo el código | lint de `check.sh` | ifx corrija el constructor (reproductor `r6.f90`) |
| Módulo `muscle_tensor_minus` (fusiona sólo `operator(-)`) | ifx 2026.2 (error #6549: se pierde el `-` unario) | `src/tensors/muscle_tensor_minus.fypp`, fachada | build con ifx | ifx resuelva la fusión de genéricos |
| Tipos derivados como locales de un procedimiento `declare target`, nunca en `private(...)` | ifx (error interno) | `BATCH_LOOP` y kernels de lote | build SPIR64 | ifx compile tipos derivados en `private` |
| Asignar componentes (`t%vals = x(1:6)`) en lugar de un constructor con una sección de argumento de forma asumida | IGC (recursión infinita del JIT) | kernels | tests target en la GPU Intel (con límite de memoria) | IGC corregido |
| Exclusiones de `tools/gpu_skip_tests.sh` (4 patrones) | IGC del driver 26.05 de la UHD 770 | tests target de materiales | `check.sh` lista lo excluido | driver nuevo de Intel, o validación en NVIDIA/AMD |
| NaN con `transfer` del patrón de bits (sin módulos IEEE en kernels) | gfortran nvptx/amdgcn (llamadas IEEE sin implementación en el dispositivo) | `QUIET_NAN` | builds NVPTX/AMDGCN y auditoría | probablemente permanente (regla GPU) |
| `VEC9_TO_MAT33` / `MAT33_TO_VEC9` en lugar de `reshape` | gfortran (llama a `_gfortran_reshape_*` en el dispositivo) | `include/muscle_common.fypp` y usuarios | auditoría de imágenes | probablemente permanente |
| `MUSCLE_UDTIO=AUTO` desactiva la salida formateada | gfortran < 16.3 (locales estáticos con UDTIO y OpenMP) | `src/CMakeLists.txt` | `check.sh` compila el camino `MUSCLE_UDTIO=ON` | se retira solo: con gfortran ≥ 16.3 `AUTO` la activa |
| Temporales tipados en lugar de `associate` con resultados de función de tipo derivado | gfortran | código de tensores | compilación | gfortran corregido |
| Ejecutar programas GPU con `systemd-run` (límite de memoria) y `timeout`; nunca `ulimit -s unlimited` | IGC (un bug llegó a usar 26 GB) | `check.sh`, reglas de tests | — | práctica permanente |

Abierto, sin solución provisional: error interno de ifx en `call p%init(...)` con dos nombres del mismo tipo
(D-019.4f).

## Tests

| Problema | Dónde | Estado | Ref. |
|---|---|---|---|
| `-Wcompare-reals` (de `-Wextra`) avisa en 26 comparaciones exactas deliberadas de `src` (cero exacto antes de dividir en yield, interruptores `prm%m /= 0` en JC completo y VA, pivote nulo en `gauss_solve`/`gauss_seidel`, `.approx.` con valores no finitos dentro de `TENSOR_COMMON_PROCEDURES`) | Flags Debug de `src/CMakeLists.txt` | **Pendiente** (usuario, 2026-11-10): se silencia con `-Wno-compare-reals` para que el build quede sin avisos, pero no debe quedar siempre; retirar el flag cuando esas comparaciones se reescriban o se marquen una a una | REANALISIS_MANTENIBILIDAD 2.1 |
| Goldens bit a bit fallaban con otros flags (`Profile_unittests`: `kin_golden`, `math_golden_derivatives`) | `tests/common/muscle_test_golden.f90` | Corregido: tolerancia por categoría y modo estricto opcional | D-021 |
| El target `coverage_report` fallaba con lcov 2.x (marcadores `#line`: «mismatched end line») y medía también los tests | `CMakeLists.txt` | Corregido: `--ignore-errors inconsistent` con lcov ≥ 2 y extracción de `src/` e `include/` | ESTADO_OPTIMO 1.4 |
| Con FMA (`-march=native -O3`): registros cuyo valor exacto es 0 dan un error relativo enorme (escala ~0): `elas_golden` (`volumetric \| stress_dev`), `yield_golden` (yld2004 `s31`) | goldens de elasticidad y yield | Resuelto (plan R): el test se reemplazó por tests unitarios con tolerancia por kind | D-021 |
| Con FMA: `math_golden_spectral_derivs` compara autovectores de autovalores repetidos (no únicos; error 1.57 con una base válida distinta) | `tests/math/test_math_golden_spectral_derivs.f90` | Resuelto (plan R): el test se reemplazó por tests unitarios de invariantes (V Λ Vᵀ, ortogonalidad); `math_reference` guarda 2 empates exactos que con FMA podrían elegir otra base válida (pendiente de medir con `-march=native`) | D-021 |
| Con FMA: `kin_golden` `nearly singular general \| pull_back_tangent` 2.7e-8 (mal condicionado) | `tests/kinematics/test_kin_golden.fypp` | Resuelto (plan R): el test se reemplazó por tests unitarios con tolerancia por kind | D-021 |
| Con FMA: tangente numérica de `closest_point_time` (jcfull + vocemod) cambia 1e3 relativo y aparece NaN: diferencias finitas en un punto no suave (piso de JC) | `tests/materials/test_mat_closest_point_time.fypp` | Pendiente: revisar el caso (paso de las diferencias o punto de evaluación) | D-021 |
| Con FMA: diferencias finitas de math 1.4e-9 (> 1e-9) y `math_golden_operations` (derivadas por diferencias en la categoría álgebra) 2.6e-10 | goldens de math | Resuelto (plan R): el test se reemplazó por tests unitarios con tolerancia por kind | D-021 |
| Los tests de compilación fallida (parámetros privados) pasaban con cualquier error porque la ruta del archivo contiene "private" | `tests/plasticity/yield_criteria`, `tests/elasticity` | Corregido: se exige "PRIVATE component" y fallan con "Cannot open module/Fatal Error" (sólo gfortran) | — |
| `check_approx` (escalares y arrays reales) aceptaba cualquier valor si un operando era Inf: la escala `max(1, |got|, |expected|)` salía infinita (encontrado por el chequeo de mutaciones del plan R: `C*dx` → `C/dx` en JC devolvía Inf y pasaba) | `tests/common/muscle_test_utils.fypp` | Corregido (R3): un componente no finito sólo pasa si es igual al esperado; la escala de un array usa sólo los componentes finitos. Ningún test escondía un Inf | ESTADO_OPTIMO R3 |
| `eigenvals`: la rama `abs(1 - d) <= EPS` (raíz media ≈ media de las otras) da autovalores ~1e-12 relativo lejos de la forma cerrada en `diag(5, 3 + 1e-11, 1)`; sin la rama la forma cerrada acierta (y sigue exacta en los casos exactos) | `src/math/muscle_math_operations.fypp:77` | **Marcado, sin corregir** (usuario, 2026-11-09: queda pendiente; corregirlo cambia resultados). Mientras tanto la referencia mínima `math_reference` la protege | ESTADO_OPTIMO R2 |
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
