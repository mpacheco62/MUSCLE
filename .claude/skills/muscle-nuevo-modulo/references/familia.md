# Familia de modelos nueva

Una familia nueva toca la arquitectura (registro, plantilla de material, solvers). Antes de escribir código,
**presenta al usuario un plan** con: el contrato de la familia (funciones genéricas, argumentos y su orden, qué
derivadas), si compone otra familia, qué solvers la usan y cómo, y qué cambia en el estado del material
(`NSTATE`, `PACK`/`UNPACK`: cambiar el layout de `hsv` de materiales existentes requiere decisión). Registra la
decisión en `docs/fypp/DECISIONES.md`.

## Piezas (además de las de un modelo, `references/modelo.md`, para el primer modelo de la familia)

1. **Prefijo** corto y nuevo (`dmg`, `kh`, ...). Módulos `muscle_<prefijo>_<modelo>` y módulo de familia
   `muscle_<prefijo>_laws` (D-027), en una carpeta de `src/` acorde (puede ser nueva; con su `CMakeLists.txt` y
   `add_subdirectory`).
2. **Módulo de familia**: cabecera con la lista de genéricos (firma de cada uno), el párrafo «To add a model» como
   en las demás familias, un `use` por modelo, `public` y una `interface` por función genérica.
3. **`tools/fypp/muscle_meta.py`:**
   - `<FAMILIA>_MODELS` y listas de rasgos que hagan falta;
   - entrada en `MODEL_FAMILIES` (`ModelFamily(prefix=, module='muscle_<prefijo>_{}', family_module=, models=,
     dims=, composition=, composed=)`); `MATERIAL_FIELDS` y el campo nuevo de `Material` si los materiales la usan;
   - `MODEL_PARAMS`, `MODEL_LABELS`, y `MATERIAL_SOLVERS` de los solvers que la llamen.
   Corre `tests/tools` (`ctest -R tools-`): comprueban registro contra archivos, parámetros contra tipos, listas
   escritas a mano contra el registro.
4. **Solvers y plantilla de material**: los cuerpos que la usan llaman sus genéricos; la plantilla de material
   agrega el componente al tipo (automático por `material_components` si `MATERIAL_FIELDS` y `MATERIAL_SOLVERS` la
   incluyen) y, si tiene estado propio, su parte de `PACK`/`UNPACK` en `muscle_material_kernel.fypp` (documentada
   en el índice de la cabecera).
5. **Tests**: `tests/<carpeta>/test_muscle_<prefijo>_laws_kinds.fypp` (genéricos contra específicos bit a bit,
   formas cerradas, derivadas vs diferencias finitas, sp vs dp, un procedimiento `@:DEVICE()`); test de GPU de
   familia: el prefijo en `FAMILY_PREFIXES` y sus firmas en `FAMILY_SIGNATURES` de `tests/gpu/target_cases.py`.
6. **Mutaciones**: un área nueva `tests/mutations/<prefijo>.toml` (`sources`, `unit = "^<carpeta>-"`,
   `reference = ""`, `samples`), corrida con `--samples all`.
7. **Documentación**: fila nueva en la tabla por familia de `contributing.md`, `pages/architecture.md` si cambia el
   flujo, glosario (prefijo y token de materiales), README, `CLAUDE.md` (estructura).
