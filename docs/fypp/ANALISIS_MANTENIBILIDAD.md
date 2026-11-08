# Análisis: mantenibilidad (2026-11-08, rama `fypp/f9-readability`)

Aspecto 4 de 4 (1: `ANALISIS_COMPRENSION.md`, 2: `ANALISIS_REPETICION.md`, 3: `ANALISIS_SOBREINGENIERIA.md`).
Pregunta: ¿qué tan caro y seguro es cambiar el proyecto, y qué tan dependiente es de una persona o de una
máquina? Método: historial de git, estado de ramas, red de seguridad (incluida una corrida de cobertura),
dependencias del entorno y deuda registrada. VERIFICADO = medido; ESTIMADO = mío.

## 1. Fortalezas (VERIFICADO)
- **Puerta de verificación reproducible** (`tools/check.sh`): lint de reglas mecánicas, build `sp;dp`,
  181 tests, tres builds GPU, tests target en la GPU Intel y auditoría de imágenes de dispositivo.
- **Goldens dp bit a bit**: detectan cualquier cambio numérico. Más `tools/gen_snapshot.sh` (generado y
  tabla de símbolos), `tools/fypp_metrics.sh` y `tools/dup_report.py` para medir cambios estructurales.
- **Decisiones y problemas registrados** (`DECISIONES.md`, `PROBLEMAS.md`: 34 corregidos, 47 documentados,
  4 pendientes) y documentación de usuario buena.
- **Código de modelos legible y local** tras la Fase 9 (un modelo = un archivo).

## 2. Riesgos, ordenados por gravedad

| # | Riesgo | Evidencia (VERIFICADO) | Efecto |
|---|---|---|---|
| M1 | **Las reglas del proyecto no están versionadas** | `.claude/` y `CLAUDE.md` están en `.gitignore`: 7 reglas, 476 líneas (estilo Fortran, kinds, 2D/3D, GPU, tests, git, «FYPP mínimo»). `pages/` y `docs/fypp/` las citan 8 veces | Un clon nuevo, otra máquina o un colaborador no tiene las reglas que el lint y la revisión aplican. Es el riesgo de pérdida más grande |
| M2 | **Rama de trabajo muy alejada de las integradas** | `fypp/f9-readability` va 282 commits por delante de `main` y 146 de `dev`; diez ramas de fase encadenadas sin integrar; `main` es el código OO anterior | Lo que ve un externo (y un worktree nuevo, como se comprobó) es código viejo; el costo de integrar crece |
| M3 | **Sin integración continua** | No hay `.github/workflows` ni equivalente; la verificación depende de `check.sh` en esta máquina (GPU Intel, gfortran con offload, ifx 2026.2) | Nada verifica los cambios fuera de esta máquina; un colaborador no puede reproducir la puerta |
| M4 | **Los goldens bit a bit dependen de los flags de compilación** | Con `Profile_unittests` (también `-O0`, pero con flags de cobertura en lugar de los de Debug) fallan 2 de 181: `kin_golden` (diferencias de 1e-17 a 4e-16 relativo) y `math_golden_derivatives` (2.4e-13, diferencias finitas) | La red de seguridad vale sólo para la configuración Debug de gfortran; otra configuración, versión de compilador o máquina dará falsos fallos. No probado en Release (ESTIMADO: también fallaría) |
| M5 | **La medición de cobertura está rota** | El target `coverage_report` falla con lcov 2.x («mismatched end line» en una plantilla, por los marcadores `#line`) | No se sabe qué parte del código cubren los tests |
| M6 | **Entorno no fijado** | `environment.yml` fija sólo python, fypp ≥ 3.1 y FORD; ningún compilador. El README dice «GFortran 9+», pero el código usa OpenMP 5.1 (`has_device_addr`) y sólo se probó con 15.2; no explica cómo instalar el toolchain de offload | Una versión distinta puede romper el build o los goldens sin aviso |
| M7 | **Soluciones provisionales de compiladores sin criterio de retiro** | Reglas de lint `exponent` y `arrayctor`, `fp_exponent`/`fp_scale`, sumas explícitas en lugar de `matmul`, `muscle_tensor_minus`, `gpu_skip_tests.sh` (4 patrones), UDTIO `AUTO`, `nocontlines`. Ninguna dice qué versión la corrige ni qué test lo comprobaría | Se acumulan aunque el compilador se corrija; nadie sabrá cuándo quitarlas |
| M8 | **Costo de cambio concentrado** | Commits de la migración: mediana 8 archivos, p90 20, máximo 57. `muscle_meta.py` cambió en 42 de 103 commits; `check.sh` en 17. `hard_stress` lo usan 17 archivos de `src` y 15 de tests. Un modelo nuevo toca 6 lugares; un material nuevo, 3 pasos extra por los goldens | Cambios de firma de familia son caros; `muscle_meta.py` es el punto caliente |
| M9 | **Ruido del build** | 52 avisos «Unused variable» en el build Debug | Tapan avisos nuevos y errores reales |
| M10 | **Conocimiento en español y en documentos largos** | `DECISIONES.md` 1 204 líneas y `PLAN.md` 709, en español, con 384 referencias `D-0xx` en el código en inglés | Barrera para mantenedores que no lean español o no tengan ese contexto |

## 3. Propuestas, ordenadas por valor/esfuerzo

| Prioridad | Cambio | Esfuerzo | Ataca |
|---|---|---|---|
| 1 | **Versionar las reglas**: sacar `.claude/rules/` (o una copia en `docs/dev/`) del `.gitignore`, dejando fuera sólo lo personal (settings locales) | Muy bajo | M1 |
| 2 | **Integrar la cadena de ramas en `dev`** (el push lo hace el usuario) y decir en el README qué rama es la vigente | Bajo (decisión del usuario) | M2 |
| 3 | **CI mínima** (lint + build CPU `sp;dp` + ctest con gfortran fijado) que corra sin GPU | Bajo–medio | M3, M6 |
| 4 | **Goldens por configuración**: bit a bit sólo en la configuración de referencia (gfortran X, Debug) y tolerancia documentada en las demás (≤ 1e-14 relativo, 1e-12 para diferencias finitas); dejarlo escrito en `testing` | Bajo | M4 |
| 5 | **Arreglar la cobertura** (`--ignore-errors inconsistent` en lcov, o gcovr) y medirla una vez | Bajo | M5 |
| 6 | **Declarar versiones probadas** (gfortran 15.2, ifx 2026.2, fypp, cmake 3.20) y corregir «GFortran 9+»; guía de instalación del toolchain de offload | Bajo | M6 |
| 7 | **Registro de soluciones provisionales**: tabla con compilador, versión, archivo, test que la verifica y condición de retiro | Bajo | M7 |
| 8 | Limpiar los 52 avisos y agregar `-Werror` a los avisos nuevos en CI | Bajo | M9 |
| 9 | `check.sh --quick` y reducir el registro de familias en Python (análisis 3) | Medio | M8 |
| 10 | Resumen en inglés de las decisiones vigentes (no de la historia) para mantenedores | Medio | M10 |

## 4. Plan integrado de los cuatro análisis

Ordenado por valor/esfuerzo; las cifras de ahorro son ESTIMADAS.
1. **Inmediato (bajo esfuerzo, sin riesgo numérico):** versionar las reglas (M1); README con ramas,
   versiones probadas y build (`build_all_tests`, Ninja); `contributing.md` completo (golden, un material por
   solver, `golden_params.py`); mensajes de error que nombren el archivo a editar; glosario; limpiar avisos;
   `check.sh --quick`; arreglar la cobertura; registro de soluciones provisionales.
2. **Corto plazo (requiere decisión del usuario):** integrar las ramas en `dev`; CI mínima; política de
   goldens por configuración; vía de escape (o script) para el golden de un material nuevo.
3. **Refactors acotados (riesgo bajo, verificables con el snapshot):** evaluación conjunta de yield y
   preámbulo de solvers como texto fijo (~265 líneas, análisis 2); reducir el registro de familias en Python
   (~250–300 líneas, análisis 3); página de arquitectura con un `.f90` generado de ejemplo.
4. **Estructural (sólo si aparece la necesidad):** alternativa D (parche guardado), B o C para los
   materiales; tests heredados a tablas de casos.

## 5. Pendiente de decidir por el usuario
- Versionar las reglas de `.claude/` (y cuáles).
- Integrar las ramas y cuándo (el push es del usuario).
- CI: dónde (GitHub Actions u otro) y con qué compilador.
- Goldens por configuración: tolerancias fuera de la configuración de referencia.
- Qué bloques del plan integrado se ejecutan y en qué orden.
