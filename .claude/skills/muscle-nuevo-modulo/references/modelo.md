# Modelo nuevo de una familia existente

Guía completa: `pages/contributing.md`, *How to add a model to an existing family* (ejemplo paso a paso con una
ley de endurecimiento) y la **tabla por familia** al comienzo de esa sección. Esta referencia resume lo que la
prueba del recién llegado (agregar la ley de Perić) mostró que se olvida.

## Dónde va cada cosa (lugares a tocar)

1. **Plantilla** `src/<carpeta de la familia>/muscle_<prefijo>_<modelo>.fypp` (la tabla da la carpeta; `os` y
   `vp` comparten `plasticity/viscoplastic_laws/`, la elasticidad congelada vive en `viscoelasticity/`). Copia el
   vecino más parecido de la familia.
2. **Módulo de familia** `muscle_<prefijo>_laws`: un `use` propio (no lo agregues al `use` de otra ley), el tipo
   en `public` y el específico en **cada** `interface` de la familia. El compilador rechaza un nombre inexistente,
   pero no detecta que falte en una interfaz: revisa las tres o cuatro.
3. **`tools/fypp/muscle_meta.py`:**
   - la lista de modelos de la familia (`HARDENING_MODELS`, `VISCOPLASTIC_MODELS`, `OVERSTRESS_MODELS`,
     `YIELD_MODELS`, `ELASTICITY_MODELS`, `VISCOELASTIC_MODELS`, `KINEMATICS_MODELS`);
   - las listas de rasgos que apliquen: `*_COMPOSED` si compone una ley de endurecimiento (recibe `sy`,
     `dsy_dep`), `VISCOPLASTIC_FLOOR` si tiene piso de tasa, `YIELD_JOINT` si escribe su evaluación conjunta,
     `YIELD_ISOTROPIC`/`ELASTICITY_ISOTROPIC`, `KINEMATICS_SPATIAL`/`_ROTATION`/`_FE_ROTATED`;
   - `MODEL_PARAMS` (mismos nombres, orden, tipo y valor por defecto que el tipo Fortran: lo comprueba
     `tests/tools/test_muscle_model_params.py`), `MODEL_LABELS` (nombre legible), `MODEL_INITS` sólo si tiene
     `init` propio (parámetros derivados o validación), `MODEL_DOMAINS` si `mat_init` debe comprobar una
     condición entre parámetros.
   Cambiar `muscle_meta.py` reconstruye todo (~1 min): es normal.
4. **CMake:** la plantilla en `muscle_add_fypp_sources` del `CMakeLists.txt` de la carpeta, antes del módulo de
   familia.
5. **Tests de la ley** en el test de la familia (columna «Test of the model» de la tabla): una función
   `test_<modelo>_${k}$()` registrada con `call record(...)`. Valores contra forma cerrada calculada en el test
   (en `dp`), derivadas contra diferencias finitas con `check_scaled`; bordes del modelo (pisos, recortes,
   `epd = 0`, `ep = 0`) en el test de límites si la familia lo tiene (`test_muscle_vp_laws_limits`).
6. **Materiales del registro de prueba** `tools/fypp/muscle_materials.py`: un material por **cada solver que usa
   la familia** (tabla; lo exige `tests/tools/test_muscle_materials.py`). Para cada uno:
   - conjuntos de parámetros en `CASE_PARAMS` de `tests/materials/case_params.py` (también los usa el test de
     GPU de la familia, que incluye el modelo solo);
   - un `solver_case(0, ...)` en `tests/materials/test_muscle_material_<solver>.fypp` (campos en
     `mat_sequence_common.fypp`: número de caso 0 = sin registros de referencia, material, elasticidad, nombres de
     los conjuntos, `iter_nw`, `dt`, deformación de fluencia de la secuencia ≈ y0/E). Ese test compara sp con dp
     y la tangente con diferencias finitas en cada paso plástico.
7. **README:** el modelo en la lista de modelos (nada lo comprueba).

## Errores típicos

- Un `sed` que agrega el nombre al `use` de otra ley rompe las interfaces: edita a mano.
- Ley compuesta (recibe `sy`, `dsy_dep` de una ley de endurecimiento): el argumento `ep` igual forma parte del
  contrato aunque no se use (los argumentos sin usar no generan aviso por diseño).
- Divisiones por un parámetro o por `sy`: documenta el dominio; si el usuario quiere validación, va en el `init` o
  en `MODEL_DOMAINS` (decisión de diseño: pregunta).
- Área de mutaciones: la de la columna de la tabla (las leyes `os` están en el área `vp`).
