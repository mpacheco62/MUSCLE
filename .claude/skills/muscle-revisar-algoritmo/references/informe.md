# Informe de la revisión

Español, con tablas, primero lo que el usuario necesita para decidir. Precedente: `docs/fypp/ESPECTRAL_F0.md`.

## Estructura

1. **Qué se revisó y cómo**: la rutina, sus consumidores, la verdad usada, los casos, los builds, la máquina (CPU,
   GPU, compiladores) y lo que no se pudo medir (y por qué).
2. **Precisión**: tabla por familia de casos, error en ε con su escala; casos que fallan o que levantan excepciones.
3. **Tests y mutaciones**: tests existentes y su verdad; mutaciones detectadas/sobrevivientes, cada sobreviviente
   con su clasificación (equivalente, bajo el redondeo, empate, hueco de test).
4. **Costo**: por llamada; kernel completo en CPU (relativo a hoy, con el control); GPU (tiempo, relativo,
   fallos GPU = CPU).
5. **Kernel en GPU**: memoria privada, spill, stack calls, FP64 emulado, `alloca` dinámicos; riesgos de la tabla de
   `references/kernel_gpu.md` que aplican.
6. **Alternativas** (si se midieron): una fila por candidato con precisión, costo CPU/GPU, kernel y tests.
7. **Recomendaciones**, priorizadas: qué cambiar, por qué (la evidencia de arriba), cuánto cuesta (código, riesgo,
   tests, decisión de API) y qué decisión necesita del usuario. Separa lo **comprobado** de lo que es **hipótesis**
   y di qué medición la confirmaría.

## Reglas

- Ningún número sin su procedimiento (comando, casos, rondas). Si una medición se perdió o salió dudosa (ruido,
  corridas concurrentes en la GPU), repítela o dilo.
- No conviertas una hipótesis en conclusión: «parece ser el tamaño del kernel» no es «es el tamaño del kernel»
  hasta que un experimento lo separe (D-030 se encontró así: la hipótesis del tamaño era falsa; era la pila).
- Los defectos nuevos van a `docs/fypp/PROBLEMAS.md`; una recomendación que el usuario acepte, a una decisión en
  `docs/fypp/DECISIONES.md` antes de tocar `src/`. Si el usuario quiere guardar el informe:
  `docs/fypp/REVISION_<rutina>.md`.
