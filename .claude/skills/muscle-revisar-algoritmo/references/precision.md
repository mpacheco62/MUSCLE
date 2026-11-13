# Precisión y robustez

## La verdad

- **Construida, no calculada con el código revisado.** Para una descomposición: A = Q diag(λ) Qᵀ con λ y Q conocidos,
  armada en cuádruple (`real128`) y redondeada al kind; para una función: su forma cerrada o una serie en alta
  precisión; si no hay forma cerrada, `mpmath` con 50 dígitos (`mp.mp.dps = 50`; `mp.eigsy`, `mp.log`, ...) a partir
  de **la entrada ya redondeada al kind** (el error de redondear la entrada no es del algoritmo).
- Precedentes: `tests/math/test_muscle_math_eigen_projections.fypp` (verdad en cuádruple),
  `benchmarks/math/bench_spectral.fypp` (banco de precisión contra cuádruple con un generador Park–Miller entero,
  determinista), `tests/kinematics/test_muscle_kin_finite.fypp` (`test_near_identity`: verdad de la F redondeada).
- Un generador aleatorio en cuádruple o con `random_number` cambia entre compiladores: usa Park–Miller entero
  (`muscle_test_helpers`) o semillas fijas documentadas.

## Las métricas

- Error en unidades de ε del kind (`epsilon(1.0_k)`), relativo a la **escala que importa**: max|λ| para autovalores,
  max|λ − media| (el desviador) si lo que interesa es la parte desviadora, |dev| para una deformación casi nula.
  Medir relativo a max|λ| esconde la pérdida del desviador.
- Autovectores y proyecciones: el error natural es ~ε/brecha; una forma cerrada que pierde ~ε/brecha² es peor. Con
  una raíz repetida sólo la suma de las proyecciones de ese autoespacio está definida: compara la suma.
- Separa el algoritmo de la condición: la deformación de Hencky a partir de C = FᵀF pierde ~ε/|F − I| por formar C,
  con cualquier descomposición. Un error del orden del límite del problema no es un defecto; uno 10–100 veces mayor
  sí.
- Reporta máximos por familia de casos (no un máximo global que mezcle escalas) y la mediana si el máximo es un
  caso aislado.

## El catálogo de casos

- Bien separados; pares con brecha 1e-1 … 1e-15 y 0, arriba, abajo y al medio; raíz triple; rotaciones distintas
  (el orden en que un algoritmo iterativo deja los valores depende de la rotación).
- Casi identidad (C = I + 1e-11 en dp): ahí una forma cerrada con umbral absoluto pierde todo el desviador.
- Media grande frente al desviador (c I + D con c/|D| = 1e6 en dp, 1e3 en sp).
- Escalas extremas: dp 1e±30 (y cerca de `tiny`/`huge`), sp 1e±18; homogeneidad con 2^±k es exacta y da una
  verdad gratis (f(2^k x) = 2^k f(x) bit a bit si el algoritmo escala por potencias de 2).
- Ceros, multiplicidad exacta de la identidad, NaN e Inf en cada posición (un NaN debe salir como NaN: `maxval`
  salta los NaN y puede esconderlos), entradas que piden recortes o pisos.
- Lo que el consumidor ve: el tensor transformado de un criterio anisótropo, el C de una cinemática; un caso
  hidrostático con ruido de redondeo.

## Excepciones de punto flotante

Compila el banco o el test con `-ffpe-trap=invalid,zero,overflow` (gfortran) y cuenta las banderas con
`ieee_exceptions` en las escalas extremas. Un overflow intermedio que da el resultado correcto sin trampa (p. ej.
`theta*theta` de una rotación, D-028) es un defecto: aborta con trampas y en GPU es un Inf. Las comparaciones con
NaN levantan `invalid`: prueba la propagación de NaN en un test **sin** trampas.

## sp y dp

Todo en los dos kinds. Si una diferencia aparece sólo en sp, mira primero si la entrada ya perdió la información
(C = I en sp por debajo de ~1e-7) antes de culpar al algoritmo. Los umbrales que dependen del kind van en
`muscle_kinds`.
