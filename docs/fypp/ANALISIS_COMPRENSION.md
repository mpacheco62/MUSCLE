# Análisis: cómo mejorar la comprensión del código (2026-11-04, rama `fypp/f9-readability`)

Alcance: aspecto 1 de 4 del análisis del usuario (entendible; quedan código repetido, sobreingeniería y
mantenibilidad). Método: lectura de muestras con cifras, una **prueba del recién llegado** (un agente sin
contexto agrega una ley de endurecimiento siguiendo sólo la documentación pública) y la investigación de la
alternativa «D» (tipos vacíos en lugar de condicionales de variante). Lo VERIFICADO se distingue de lo
ESTIMADO; no se cambió código en esta etapa.

## 1. Diagnóstico

| Zona | Estado | Evidencia |
|---|---|---|
| Modelos (hard, elas, vp/os, visco, yield, kin) | **Se entiende** | Tipo, firma, `intent`, doc, cuerpo y cierre a la vista; 0 llamadas a macros prohibidas |
| Documentación de usuario (README, integration, tensors) | Buena | Ejemplos para Abaqus, ANSYS, LS-DYNA y lotes en GPU |
| Generación de materiales | **Opaca** | 4 plantillas, 1147 líneas, 83 `#:if` de variante; 19 variables externas en el kernel; ver sección 3 |
| `muscle_meta.py` | Menos grave de lo que parecía | 2287 líneas = 971 de funciones (98, la mayor con 46 líneas), 719 de datos, 568 de comentarios |
| Guía de colaboración (`pages/contributing.md`) | **Incompleta** (prueba del recién llegado) | Ver sección 2 |
| Tests | Grandes y monolíticos | `test_yield_kinds` 1057 líneas, `test_hard_laws_kinds` ~440; 74 plantillas + 49 `.f90` |
| Referencias a decisiones en código | Ruido | 384 menciones `D-0xx` en `src/`; las más densas: `muscle_material` (24), `kin_finite` (22), `material_kernel` (20), `muscle_kinds` (20) |
| Nombres del registro | Crípticos sin glosario | `cp_`, `cp2_`, `pz_`, `tm_`, `cut_`, `vm`, `hsv`, `SoA`: ninguna página los define; 80 materiales |
| Mapa de arquitectura | No existe | Sólo hay una cabecera por plantilla; el flujo registro → plantilla → solver → modelos → API no está dibujado |
| Ruido en el build | Estorba | ~52 avisos `-Wunused-dummy-argument` (argumentos de kernel intencionales) tapan los errores reales |
| Contrato de variables de los includes | **Bien documentado** | Cada include lista en su cabecera las variables que espera; sólo `closest_point_2_body` (7 variables) y `cutting_plane_body` (5) omiten algunas |

## 2. Prueba del recién llegado (VERIFICADO)

Tarea: agregar la ley `quadratic` y un material de prueba siguiendo `pages/contributing.md` en la rama
vigente. Resultado: **parcial**. La ley, su test (sp y dp), el módulo de material generado y el build
funcionan; 3 tests existentes fallan porque el material nuevo no tiene golden. ~25 comandos, ~14 archivos,
dos recompilaciones fallidas antes de la tercera. Fricciones, por impacto:

1. **Registrar un material exige pasos no documentados**: una entrada en `tests/materials/golden_params.py`
   (el mensaje de fypp nombra el síntoma pero no el archivo), un material por cada solver que tenga la
   familia (`cp_`, `cp2_` y `cut_`: `tests/tools/test_muscle_materials.py:47`), un `golden_case(...)` en
   `tests/materials/test_mat_*.fypp` y regenerar referencias (`tests/reference/write_golden.sh`).
2. El paso 2 (agregar el modelo a la familia) es fácil de hacer mal: el `use` va en su propia línea y un
   `sed` ingenuo rompe las interfaces. El compilador lo señala bien.
3. El paso 5 (tests) no dice si ampliar `test_hard_laws_kinds.fypp` (monolito de ~440 líneas con leyes
   fijas) o crear un test pequeño por ley.
4. `contributing.md` remite a `.claude/` y `docs/fypp/`; `.claude/` está en `.gitignore`, así que un
   colaborador externo no lo ve.
5. `muscle_meta.py` no dice dónde editar (sin índice de secciones); `props` y el resto de su CLI no
   figuran en la guía.
6. Un fallo temprano del build deja varios tests «Not Run» y parece que fallan muchos; el log de avisos
   dificulta ver `FAILED`.
7. El README no dice que `main` es el layout anterior y que la migración vive en otra rama (una primera
   prueba, hecha sobre `main` por error de encargo, lo mostró); no menciona `build_all_tests` ni Ninja.
(Una octava fricción reportada, la ruta de `muscle_fd.fypp`, resultó un falso positivo.)

## 3. La alternativa «D» (tipos vacíos), investigada

Idea: que la ausencia de una característica sea un tipo vacío y que los genéricos por tipo reemplacen los
`#:if`. Se aplicó de verdad en dos worktrees (parche en `build-snapshots/d_visco_kin.patch`; también el de
sólo VISCO en el scratchpad). VERIFICADO salvo lo marcado:

| Medida | Base | VISCO | VISCO+KIN |
|---|---|---|---|
| `#:if` en plantillas de material | 83 | 69 | 55 |
| Líneas de plantillas de material | 1147 | 1169 | 1110 |
| Líneas nuevas en módulos propios (0 `#:if`) | 0 | 391 | 865 |
| Goldens dp, 8 ejecutables | base | idénticos | idénticos |
| Tests (gfortran, sp;dp) | 181 | 181 | 181 |
| Tests target Intel UHD 770 | 40/40 | 40/40 | no corrido completo |
| Tiempo de compilación de la librería | 5.0 s / 8.6 s (Debug / Release) | sin diferencia | sin diferencia |
| `libmuscle.a` Release | 4.38 MB | +7.5 % | +8.6 % |
| Spill sp en Intel, Release | 0 | 0 | 0 (el aumento en Debug fue artefacto de `-O0`) |

Hallazgos adicionales: un tipo vacío no puede ir como `firstprivate` de una región `target` en Intel
(falla al lanzar); como local de un procedimiento device, con `map(to:)` o con un `integer` dummy funciona.
Con D, VISCO×KIN se combinan registrando un material, sin plantilla nueva (compilan 5 combinaciones de
prueba en un registro aparte; contra un material equivalente sin cinemática, small strain + Maxwell difiere
6e-14 y corotacional + Maxwell 0; las demás no tienen referencia).

**Mi lectura (difiere en parte de la del agente, que recomendó D para VISCO+KIN):**
- El beneficio medible es modesto: −34 % de `#:if` y quedan 55, de los cuales ~8 (firma `eps(6)` frente a
  `F(9)`) son irreducibles por tipos. VISCO solo cuesta 391 líneas por −14 condicionales.
- **Choca con el criterio de D-020**: el eje KIN necesita 64 subprogramas específicos por pareja
  (`flow_enter`, `flow_spin`, `rate_add` por kind) **generados con bucles sobre listas de metadatos**, y el
  tipo de cada componente lo calcula Python (`material_flow`, `material_rate`). Es el patrón que la
  fase acaba de retirar, desplazado a otro archivo.
- Quita la validación «kin + visco no soportado» y por tanto **agrega funcionalidad** sin validación
  física (las combinaciones sólo se compararon contra sí mismas; que TL difiera 1–3 % es una
  afirmación del agente, no comprobada). Eso queda fuera de «sin funcionalidad nueva».
- Cambia el tipo de opacidad: el kernel deja de tener `#:if`, pero «qué función se ejecuta» pasa a depender
  de tipos definidos en otros módulos (`muscle_kin_flow.fypp`, 474 líneas, casi todo documentación).
- Sin medir: ejecución en NVIDIA/AMD y los tests target Intel de la rama KIN completa.
- **Recomendación: no adoptar D ahora.** Primero la documentación y los mensajes de error (sección 4),
  que atacan la fricción real medida. Retomar D (parche listo) sólo si combinar viscoelasticidad con
  cinemática finita pasa a ser un requisito; en ese caso, decidir antes cómo evitar los específicos por
  pareja generados desde metadatos.

## 4. Plan de mejora ordenado por valor/esfuerzo

**Bajo esfuerzo, impacto alto (propuesto como siguiente paso)**
1. `pages/contributing.md`: lista de verificación completa para un modelo y para un material (CMake,
   `golden_params.py`, un material por solver, goldens, tests), el aviso del `use` en su línea, el flujo de
   build (Ninja, `build_all_tests`, `ctest -R`), y quitar o hacer autocontenidas las remisiones a `.claude/`.
2. Mensajes que dicen qué archivo editar: `tests/materials/props_samples.py:60` («agrega una entrada en
   `tests/materials/golden_params.py`») y `tests/tools/test_muscle_materials.py:47` («agrega un material para
   el solver X»).
3. README: mapa de ramas y estado (`main` = layout anterior; migración en `fypp/f9-readability`).
4. Silenciar `-Wunused-dummy-argument` en kernels (o consumirlos con un patrón explícito) para que los
   errores no queden enterrados.
5. Glosario (`cp`, `cp2`, `pz`, `tm`, `cut`, `vm`, `hsv`, `SoA`, `os`, `vp`, `kin`) en una página.

**Esfuerzo medio**
6. `pages/architecture.md`: diagrama registro → plantilla de material → solver → modelos → API de punto y de
   lote, más el `.f90` generado de un material de ejemplo (regenerado en `check.sh` y comparado).
7. Un test pequeño por ley de endurecimiento (patrón copiable) en vez de ampliar el monolito
   `test_hard_laws_kinds`; plantilla con valores, derivadas y diferencias finitas. Lo mismo, a su tiempo,
   para yield y vp.
8. Índice de secciones al inicio de `muscle_meta.py` (con qué plantilla consume cada bloque); evaluar partirlo
   en registro / nombres / generadores de texto. Las funciones ya son cortas: el aporte es de navegación.
9. Reducir las referencias `D-0xx` en `src/` a las que explican un porqué no obvio.
10. Completar las cabeceras de `closest_point_2_body` y `cutting_plane_body` con las variables que usan.

**Esfuerzo alto / a decidir**
11. D (sección 3), B (esqueletos por variante, ~2400 líneas repetidas) o C. Ninguno antes de 1–7.
12. Un andamio (`tools/new_hard_law.py`) que edite los 6 lugares de un modelo nuevo: sólo si el ritmo de
    modelos nuevos lo justifica (es más código que mantener).

## 5. Pendiente de decidir por el usuario
- Si se aplican los puntos 1–5 (bajo esfuerzo) en esta rama.
- Si D se archiva con el parche o se prueba con un caso de uso real (kin×visco).
- Si los materiales nuevos pueden registrarse sin golden (opción de escape) o la regla se mantiene.
