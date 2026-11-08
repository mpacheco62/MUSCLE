# Análisis: código repetido (2026-11-04, rama `fypp/f9-readability`)

Aspecto 2 de 4 (1: `ANALISIS_COMPRENSION.md`; quedan sobreingeniería y mantenibilidad). Método:
`tools/dup_report.py` (bloques de ≥ 6 líneas normalizadas que aparecen más de una vez) más lectura de los
focos. **Advertencia de método:** el detector normaliza números y los nombres de kind, así que el código
desenrollado que sólo difiere en índices cuenta como repetido; las cifras son una **cota superior** de las
copias textuales. Mide texto repetido, no repetición de información (sección 3).

## 1. Cifras (VERIFICADO)

| Zona | Líneas de código | Redundantes | Dentro del mismo archivo | Entre archivos |
|---|---|---|---|---|
| `src/` + `include/` | 9 028 | 1 122 (12 %) | 506 | 616 |
| `tests/` | 29 395 | 6 454 (22 %) | 4 468 | 1 986 |
| `tools/` | 2 262 | 0 | – | – |

Expansión por plantillas (no es repetición de fuente): 15 479 líneas de plantilla en `src`+`include` generan
188 057 de Fortran. Los 80 módulos de material salen de 4 plantillas (1 147 líneas) y suman 148 886 líneas
generadas. La repetición sp/dp y 2D/3D está resuelta en el fuente con `#:for`.

**Conclusión de nivel:** el código de la librería está poco repetido (12 % como cota); la repetición pesa
en los **tests** y es sobre todo interna a cada archivo.

## 2. Dónde está, clasificado

### 2.1 Justificada o por diseño (no tocar)
- **Cuerpos de tensores desenrollados a mano** (`tensors/operators/*`, `tensors/types/*`): decisión de
  rendimiento registrada (reglas de estilo). Es la mayor parte de los bloques de 6 líneas ×17.
- **Un `pure function` por ley y por función de familia** (hard, vp, visco): la forma repetida es el
  contrato de la familia; lo que cambia es la fórmula.

### 2.2 Evitable, de rendimiento moderado (candidatos)
| # | Foco | Evidencia | Posible corrección | Ahorro estimado |
|---|---|---|---|---|
| R1 | **Evaluación conjunta del criterio de fluencia** (`yield_joint_evaluation_*` + `yield_stress_eq_and_derivatives_*`) | 5 copias de 36 líneas idénticas salvo el nombre de modelo (`vonmises`, `hill48`, `druckerprager`, `cpb06`, `cpb06ex2`; `yld2004` tiene la suya) | Macro de **texto fijo con sustitución de nombres** (permitida por D-020), en un include de un nivel | ~145 líneas |
| R2 | **Preámbulo y cierre de los cuerpos de solver** | 7 líneas ×5 solvers, 21 + 15 líneas entre `closest_point` y `closest_point_2`, 30 entre `perzyna` y `time` | Compartir fragmentos de texto fijo (un nivel de include) | ~120 líneas |
| R3 | **Cuerpos 3D y 2D escritos a mano** en `elasticity_linear` y `visco_frozen` (la Fase 9 los desenrolló a mano: 123 → 314 y 110 → 265 líneas) | 114/198 y 105/185 líneas en bloques repetidos | Volver a desenrollar desde `layout()` (lógica en plantilla, que D-020 desaconseja) o mantener | Decisión de estilo: ~180 líneas a cambio de lógica en Python |
| R4 | **Tipos de tensores** (`2o`/`2osym`, `4o`/`4o2sym`/`4o3sym`) | bloques de 15–19 líneas repetidos entre familias | Mover a `TENSOR_COMMON_*` lo que sea idéntico (ya es el criterio) | ~60 líneas |
| R5 | **Andamio de firmas de operadores** | bloque de 6 líneas ×17 en los 6 módulos de operadores | Ya se escribió a mano por D-020; dejarlo | 0 (decisión tomada) |

### 2.3 Tests: la mayor bolsa (22 %)
- Los tests `.f90` **heredados** (49 archivos, 14 277 líneas, sólo dp): `test_muscle_math_spectral_derivs.f90`
  tiene 80 % de sus líneas en bloques repetidos (casos casi idénticos), `test_muscle_yield_vonmises.f90` 83 %,
  `test_muscle_hard_swift.f90` y `test_muscle_hard_ludwik.f90` ~91 líneas redundantes cada uno. Son la red de
  seguridad de la migración (la regla pedía no modificarlos). Pasarlos a tablas de casos con un bucle
  reduciría miles de líneas, pero se pierde el detalle de qué caso falla salvo que la tabla lo etiquete.
- Los tests por plantilla (`*_kinds.fypp`, `ops_golden_cases.fypp`) repiten estructura por familia
  (345/775 líneas repetidas en `4o2sym`, 245/491 en `ops_golden_cases`); parte ya la absorben los helpers de
  `tests/common/`.
- `tools/`: 0 % de repetición.

## 3. Repetición de información (más importante que la de texto)

El detector no la ve. Un modelo nuevo se declara en **6 lugares** (lo mostró la prueba del recién
llegado): la plantilla del modelo, el módulo de familia (`use`, `public`, una `interface` por función),
`MODEL_PARAMS` en `muscle_meta.py`, `CMakeLists`, los tests y el registro de materiales con su golden. Tres
de esos lugares repiten la *lista de modelos de la familia*:
- `muscle_hard_laws.fypp` (61 líneas), `muscle_yield_criteria.fypp` (136), `muscle_vp_laws.fypp` (110):
  cada modelo aparece 6–12 veces por módulo de familia.
- `MODEL_PARAMS` repite los nombres de los componentes del tipo (cubierto por
  `tests/tools/test_muscle_model_params.py`, que compara nombre, orden, tipo y valor por defecto).
- Alternativa que D-020 permite («listas de nombres por patrón»): generar los `use`/`public`/`interface` del
  módulo de familia con `#:for m in <modelos de la familia>`. Hoy se prefirió escribirlos a mano para que el
  fuente se lea completo; es un compromiso explícito entre **explicitud** y **un solo lugar por dato**.

## 4. Recomendación

1. **R1** y **R2**: son texto fijo con sustitución de nombres y encajan en la regla; juntos ~265 líneas y
   una sola copia de cada fragmento. Riesgo bajo (el snapshot del generado debe salir idéntico).
2. **Mantener** R3 y R5 como están (decisión de estilo de la Fase 9) salvo que se prefiera la lógica en
   plantilla; **R4** sólo donde sea texto idéntico.
3. **Tests:** no tocar los `.f90` heredados ahora. Si se quiere reducir, empezar por
   `test_muscle_math_spectral_derivs.f90` y `test_muscle_yield_vonmises.f90` con tablas de casos
   etiquetadas, verificando que los mismos casos siguen cubiertos.
4. **Información:** decidir entre explicitud y lista única en los módulos de familia (sección 3); si se
   elige la lista única, el cambio es pequeño y el compilador sigue verificando los nombres.
5. Repetir `tools/dup_report.py src,include 20` después de cada cambio para medir.

## 5. Pendiente de decidir por el usuario
- Si se hacen R1 y R2 (ahorro ~265 líneas, riesgo bajo).
- Explicitud frente a lista única en los módulos de familia.
- Si se pasan los tests heredados a tablas de casos, y cuáles primero.
