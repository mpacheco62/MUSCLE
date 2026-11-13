# Módulo de soporte (math, tensores, utilidades de kernel)

Un módulo sin familia: una rutina de álgebra (`src/math/`), una operación de tensores (`src/tensors/`), una utilidad
que usan los kernels.

## Reglas específicas

- **Nombre** `muscle_<área>_<cosa>` (archivo con el mismo nombre). Funciones públicas en **snake_case** con un
  prefijo del área o descriptivo: no sigas los nombres heredados en mayúsculas o camelCase de math (`M33INV`,
  `dEigenvalues_dTensor`), son una excepción histórica documentada (re-análisis de convenciones, V7).
- **Errores**: en kernel, un entero `status` con `STATUS_*` (o, para rutinas de álgebra, el par `ok` + `cond` como
  las inversas de tensores); en host puede devolver un código, documentado. No inventes un quinto canal.
- **Si corre por punto** (lo llama un modelo, solver o material): `pure`, `implicit none`, `@:DEVICE()`, sin
  construcciones prohibidas en GPU; sistemas pequeños con el LU propio (`muscle_math_lu_solve`), nunca LAPACK.
  Sobre argumentos de forma asumida (`A(:,:)`), bucles explícitos en lugar de reducciones con `dim=` o expresiones de
  arrays que creen temporales de tamaño variable (D-030, `NORM_INF`); `fp_scale`/`fp_exponent` en lugar de
  `scale`/`exponent` (D-031, D-019.5b).
  Agrega un caso a un test target (p. ej. `tests/math/test_muscle_math_target.fypp`) para que corra dentro de una
  región `!$omp target` en sp y dp y se compare con la CPU.
- **Host-only** (I/O, LAPACK opcional): marca el archivo con `! lint: host-only` y no lo uses desde kernels.
- **Tensores**: 2D y 3D desde la misma plantilla usando `layout(family, dim)` de `muscle_meta.py`; test de
  embebido 2D→3D (o contra los tensores completos como `test_muscle_tensor_ops_*_full.fypp`); no cambies el orden
  Voigt de ningún tipo; código común a todas las familias va a `include/muscle_tensor_common.fypp` sólo si es
  idéntico carácter a carácter (D-012).
- **Precisión**: umbrales por kind en `muscle_kinds` si dependen de la precisión; tests en sp y dp con casos mal
  condicionados explícitos (tolerancia ≈ 2·κ·ε documentada).
- **Tests**: `tests/<área>/test_muscle_<área>_<cosa>[_<tema>].fypp` con propiedades (inversa por producto,
  ortogonalidad, reconstrucción `V Λ Vᵀ = A`, homogeneidad) y casos borde (matriz singular, autovalores repetidos,
  entradas no cuadradas o de forma inválida → status de error).
- **Mutaciones**: agrega el archivo al área que corresponda en `tests/mutations/` (`math`, `linalg`, `types`,
  `tensors`) o crea una; las rutinas sólo en sp (`#:if k == 'sp'`) no las ve el build dp de mutaciones: dilo en el
  informe.
