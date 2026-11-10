# Re-análisis 4: mantenibilidad (2026-11-10, rama `fypp/f9-readability`, HEAD `82848cd`)

Aspecto 4 de 15 del re-análisis posterior al plan R. Compara con `ANALISIS_MANTENIBILIDAD.md` (2026-11-08).
Pregunta: ¿qué tan caro y seguro es cambiar el proyecto, y qué tan dependiente es de una persona o de una
máquina? Método: builds limpios en Debug, Release y `Profile_unittests` con `ctest` completo y cobertura
(`build-maint-*`), historial de git, estado de ramas y `PROBLEMAS.md`. Todo VERIFICADO salvo lo marcado.

## 1. Los riesgos del análisis anterior

| # | Riesgo (2026-11-08) | Hoy | Estado |
|---|---|---|---|
| M1 | Reglas del proyecto sin versionar | Igual (`.claude/`, `CLAUDE.md` ignorados) | **Aceptado** por el usuario (0.1) |
| M2 | Rama alejada de las integradas | 325 commits por delante de `main`, 189 de `dev` (antes 282/146); 12 ramas | **Aceptado** por ahora (0.2); crece ~40 commits por semana |
| M3 | Sin CI | Igual | Pospuesto (1.5) |
| M4 | Referencias dependientes de flags | **169/169 en Debug, Release y `Profile_unittests`** (builds limpios) | Resuelto (D-021 + plan R) |
| M5 | Cobertura rota | Funciona: **81.8 % de funciones** (3 749/4 582), 97.3 % líneas, 90.8 % ramas (antes 81.7 / 98.2 / 91.9) | Resuelto. Líneas y ramas bajaron ~1 punto al retirar referencias (no se midió qué archivos) |
| M6 | Entorno sin fijar | Versiones probadas en el README | Resuelto (0.3); el entorno sigue sin fijar compiladores |
| M7 | Soluciones provisionales sin criterio de retiro | Tabla con 14 filas en `PROBLEMAS.md` | Resuelto (2.6); ninguna se pudo retirar (mismos compiladores) |
| M8 | Costo de cambio concentrado | Desde 2026-11-05: 43 commits, mediana 3 archivos, p90 18, máx. 71 (renombre del plan R). `muscle_meta.py` cambió en 4 (antes 42 de 103). Lo más tocado es la documentación (`ESTADO_OPTIMO.md` 21 veces) | **Mejor**: el punto caliente pasó del código a los documentos |
| M9 | Ruido del build | **No resuelto: ver 2.1** | Error de medición en 1.6 |
| M10 | Conocimiento en español y en documentos largos | `docs/fypp/` 2 995 líneas (`DECISIONES.md` 1 294) | Igual; crece |

## 2. Hallazgos nuevos

### 2.1 El build Debug tiene 514 avisos, no 0 (corrige 1.6)
Build limpio Debug (`-Wall -Wextra`, `src/CMakeLists.txt:91`): **514 avisos, 504 en `src/`** y 78 en `tests/`
(algunas líneas cuentan en ambos). La afirmación de 1.6 («la librería no tiene avisos», criterio «Avisos del
build Debug: 0» en `ESTADO_OPTIMO.md`) fue un error de medición: casi con seguridad se contó un build
incremental, donde la librería ya estaba compilada. Clasificación:

| Aviso | Cantidad | Origen | Lectura |
|---|---|---|---|
| `-Wunused-dummy-argument` | 393 | Operadores con identidades (`ops_addition_subtraction` 120, `iden` 36, `contraction_*` 40), cinemáticas (`kin_*` ~80: argumentos del contrato de familia que una medida no usa), `visco_kelvin` | Intencional (contrato uniforme); ruido |
| `-Wcompare-reals` | 62 | Comparaciones exactas deliberadas (`== 0`) | Intencional; ruido |
| `-Wunused-variable` | 43 | Sobre todo tests | Limpieza |
| `-Wunused-function` | 10 | `tests/materials/test_mat_return_mapping.fypp`: `fd_check_<m>` de los 5 materiales Hill48 se generan pero **no se llaman** (sólo el caso radial, de von Mises, los usa); `test_kin_kinds` (`input_P`, `vals6`) | **Posible hueco de test**: la tangente por diferencias finitas no se comprueba ahí para Hill48 (verificar en el aspecto 6 si otro test la cubre) |
| `-Wuninitialized` | 5 | Temporal interno de gfortran en `test_kin_finite.fypp:123` | Ruido del compilador |
| `-Wreturn-type` | 1 | `mul_i2o_i2o` (tipo vacío `iden_2O`, no hay nada que asignar) | Falso positivo |

Efecto: un aviso nuevo y real (p. ej. una variable sin inicializar en un kernel) queda enterrado entre 500.

### 2.2 Deuda registrada
`PROBLEMAS.md`: ~36 corregidos o resueltos, 46 documentados, 5 pendientes, 3 marcados sin corregir (decisión del
usuario). De los pendientes, uno está **obsoleto**: «los ejecutables de test no enlazan con offload
(`GOMP_offload_register_ver`), pendiente (Fase 8)» (l. 81) lo resuelve `cmake/MuscleGpu.cmake:72` y los
tests target corren. Otros dos llevan desde la Fase 9 sin fecha: `tm_nnl_vm` con deformación ×40 no converge
(l. 111) y Chaboche n=0.5 + DP en compresión (l. 115); y uno de D-021 (tangente numérica de
`closest_point_time` con FMA, l. 194).

### 2.3 Dependencia de personas y máquinas
- Autoría en `main..HEAD`: Matias Pacheco (205 + 48 + 8 con otros nombres) y Jorge Astete (64). Toda la
  migración y el re-análisis los hizo una persona: el riesgo de dependencia de una persona no cambió.
- La puerta completa (`check.sh`) sólo corre en esta máquina (GPU Intel, gfortran con offload, ifx 2026.2).
  Sin CI no hay segunda máquina (M3).
- Builds limpios (gfortran 15.2, `-j8`): Debug 71 s, Release 68 s (incluye fypp y los 81 materiales); detalle
  en el aspecto 5.

## 3. Conclusión

La red de seguridad numérica es ahora **portable entre configuraciones** (Release y cobertura en verde, antes
fallaban) y la cobertura es medible. Los riesgos que siguen abiertos son los que el usuario aceptó o pospuso
(reglas sin versionar, ramas sin integrar, sin CI), más uno nuevo que estaba mal medido: **el build tiene
~500 avisos**, casi todos intencionales, que esconden los reales. El punto caliente de cambio pasó de
`muscle_meta.py` a los documentos de proceso.

## 4. Propuestas (no aplicadas)

| Prioridad | Cambio | Efecto | Riesgo |
|---|---|---|---|
| M-a | **Hecho (2026-11-10).** **Corregir 1.6 y el criterio de avisos** en `ESTADO_OPTIMO.md` con la cifra real | Registro fiel | Nulo |
| M-b | **Hecho (2026-11-10)**: `-Wno-unused-dummy-argument` (librería y tests), `-Wno-compare-reals` sólo en la librería y pendiente (PROBLEMAS), tests limpios, `check.sh` falla con cualquier aviso (el completo recompila desde cero). **Avisos:** `-Wno-unused-dummy-argument -Wno-compare-reals` en la librería (intencionales por diseño) o marcar los argumentos; limpiar los `unused-variable/function` de tests; y que `check.sh` cuente los avisos de un build limpio y falle si aumentan | Un aviso nuevo se ve | Bajo; decidir qué avisos se silencian |
| M-c | Cerrar la fila obsoleta de offload en `PROBLEMAS.md`; poner fecha o decisión a los 3 pendientes de solvers | Deuda al día | Nulo |
| M-d | Revisar si la tangente por diferencias finitas de los materiales Hill48 está cubierta en otro test (aspecto 6) | Posible hueco | Nulo |
| M-e | Reconsiderar la integración en `dev` (M2) y la CI (M3) cuando el usuario lo decida: la distancia crece ~40 commits por semana | Menos costo de integrar | Decisión del usuario |
