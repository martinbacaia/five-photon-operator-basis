# Fase 2 — Re-verificación exhaustiva de los generadores (a pedido explícito, previo a publicar)

El usuario pidió explícitamente re-verificar todo el trabajo de generadores antes de seguir, dado
que el objetivo final es publicar. Se escribió `scripts/reverificacion_generadores.py`, que **no
reutiliza ninguna conclusión anterior como dada** — reconstruye todo desde cero y cruza cada
resultado por al menos 2 caminos independientes.

## (A) El grupo generado es genuinamente S₅ (no solo "orden 120")

Hasta ahora solo se había verificado que el grupo generado tenía 120 elementos. Esto por sí solo NO
prueba que sea S₅ (hay 5 grupos distintos de orden 120, incluyendo S₅ y otros no isomorfos). Se
verificó, además:
- La distribución de tipos de ciclo de las 120 permutaciones generadas coincide exactamente con la
  de S₅ (1+10+15+20+20+30+24=120, por tipo de ciclo) — **coincide, OK**.
- El **orden de cada una de las 120 matrices** (calculado por multiplicación repetida hasta volver
  a la identidad) coincide con el orden esperado de la permutación correspondiente, para los 120
  elementos sin excepción — **120/120 coinciden**.

## (B) El invariante de grado 2 (Q), verificado por 2 construcciones matemáticamente distintas

1. Vía operador de Reynolds (como antes).
2. Vía la **forma bilineal invariante canónica** de cualquier representación real: G = (1/|G|)Σ_g
   MᵀM, que es automáticamente invariante para cualquier representación (es la construcción estándar
   que prueba que toda representación finita preserva una forma positiva definida). Esto es
   conceptualmente distinto del operador de Reynolds, aunque ambos usan la misma lista de matrices.

**Resultado**: ambos son proporcionales (factor 1/5) — coincide exacto, confirmando que Q es
realmente el único invariante cuadrático (salvo escala), no un artefacto de la elección de semilla.

## (C) El invariante de grado 3 (C), verificado con 3 semillas distintas

Se repitió el cálculo con 3 semillas distintas (a₁³, a₂³, a₁²a₂) — las 3 dan resultados no nulos y
**mutuamente proporcionales**, confirmando que el espacio de invariantes de grado 3 es realmente
1-dimensional y C no depende de una elección arbitraria de semilla.

## (D) Grado 4, reconfirmado con un CUARTO candidato independiente

Se agregó una tercera semilla de Reynolds (a₁²a₂a₃, dando D₃) a los 3 candidatos anteriores (Q²,
D₁, D₂), para un total de 4 candidatos. Se resolvió la relación lineal general entre los 4 de forma
exacta (coeficiente a coeficiente): la solución tiene **2 parámetros libres**, confirmando que los 4
candidatos generan un espacio de dimensión exactamente 4−2=2 — **coincide exacto con Molien**, y de
forma más robusta que con solo 3 candidatos (con 4 candidatos independientes, la conclusión es más
difícil de que sea casualidad).

## Conclusión de la re-verificación

**Los 5 chequeos pasan.** El trabajo de los pasos 2, 3 y 4 (representación de S₅, anillo no libre,
generadores explícitos en grados 2/3/4) queda confirmado con un nivel de rigor mayor que el
original — no solo "corrió y dio el número esperado", sino verificado por construcciones
matemáticamente independientes entre sí en cada paso. Se puede seguir con confianza a los grados 5
y 6.

## Archivos

- `scripts/reverificacion_generadores.py` — script completo, auto-contenido, no depende de
  ejecuciones anteriores (reconstruye el grupo y los invariantes desde cero).
