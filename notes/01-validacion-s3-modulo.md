# Fase 2 — Paso 1: validación del formalismo S₃ sobre invariantes de Mandelstam

## Qué se hizo

Se implementó en Python (`scripts/s3_partition_functions.py`) una verificación **independiente**
de las ecuaciones (2.40)-(2.44) de Chowdhury-Gadde-Gopalka-Halder-Janagal-Minwalla
(arXiv:1910.14392) — la construcción central que habría que extender de n=4 a n=5 partículas para
atacar el hueco de la pregunta #44 (ver `fase1-mapa/04-profundizacion-crg.md`).

**Qué se verificó en concreto**: para cada grado n=0..8, la dimensión del espacio de polinomios
homogéneos de grado n en (s,t,u), sujetos a la restricción s+t+u=0, que transforman en cada una de
las 3 representaciones irreducibles de S₃ (1S simétrica, 1A antisimétrica, 2M mixta de dimensión 2).

**Dos métodos independientes, sin copiar la fórmula del paper**:
1. **Conteo directo por teoría de caracteres**: se construyeron los 6 elementos de S₃ como mapas
   lineales explícitos sobre (s,t) (con u=-s-t sustituido), se calculó la traza de cada uno sobre el
   espacio de polinomios homogéneos de grado n, y se aplicó la fórmula estándar de multiplicidad de
   representaciones irreducibles (proyección con caracteres).
2. **Fórmulas cerradas del paper** (ecuación 2.41, la versión de S₃ del "cycle index"/serie de
   Molien, más la ecuación 2.43 para imponer la restricción lineal).

**Resultado**: los 2 métodos coinciden exactamente en los 9 grados probados (0 a 8) — ver tabla en
la salida del script.

## Súper-verificación (a pedido explícito del usuario, 16 sept. 2026)

El usuario pidió confirmar de verdad que la reproducción es correcta, no solo que "el script corrió
y dio OK". Se escribió un segundo script, `verificacion_rigurosa.py`, que ataca la pregunta desde 4
ángulos **independientes entre sí y del script original** (ninguno de los 4 asume que el primer
script está bien):

1. **¿Los 6 mapas realmente forman una representación honesta de S₃?** Se verificó por
   composición explícita de funciones (no por confiar en las etiquetas puestas a mano): que cada
   "transposición" tiene orden 2, que cada "ciclo" tiene orden 3, que componer 2 transposiciones da
   un ciclo (propiedad estructural real de S₃, no de un grupo abeliano), y que el grupo generado por
   2 de los mapas cierra en exactamente 6 elementos distintos (ni más ni menos — se generó el grupo
   por composición repetida desde cero y se contó). **Resultado: PASA.** Esto descarta que la
   coincidencia con el paper fuera casualidad por haber armado mal la acción de grupo.
2. **Consistencia de dimensión total, sin usar la fórmula del paper para nada**: la suma
   (dim 1S)+(dim 1A)+2×(dim 2M) debe dar exactamente n+1 en cada grado n (la dimensión conocida,
   por conteo directo de monomios, de polinomios homogéneos de grado n en 2 variables). Verificado
   para n=0 a 12. **Resultado: PASA.**
3. **Comparación palabra por palabra con la prosa del paper** (no con su fórmula matemática): el
   texto de la sección 2.7 dice explícitamente, en inglés, qué representaciones "sobreviven" en los
   grados 0 a 3 tras imponer s+t+u=0 ("one mixed representation" en grado 1, "one symmetric and one
   mixed" en grado 2, "one symmetric, one mixed and one anti-symmetric" en grado 3). Se tradujo esa
   prosa a multiplicidades esperadas y se comparó contra el cálculo. **Resultado: PASA** (coincide
   exacto en los 4 grados).
4. **La fórmula cerrada (2.44) recalculada a mano**, sin usar `sympy.series` (que el script
   original sí usaba): se identificó que D=1/((1-x⁴)(1-x⁶)) es la función generatriz de conteo de
   particiones de n en partes de tamaño 2 y 3 (en la variable y=x²), y se escribió un contador de
   particiones por fuerza bruta en Python puro (sin sympy) para los coeficientes de 1S, 1A y 2M.
   Comparado contra el método de caracteres para n=0 a 12. **Resultado: PASA.**

**Los 4 chequeos pasan.** Esto ya no es "dos métodos que casualmente coinciden" — es: la acción de
grupo está bien armada (chequeo 1), el total de dimensión es consistente de forma independiente
(chequeo 2), coincide con lo que el paper dice en palabras (chequeo 3, no solo en fórmulas), y la
fórmula cerrada fue re-derivada por un tercer camino completamente distinto (chequeo 4, conteo de
particiones a mano). Da confianza real de que el formalismo de n=4 está entendido correctamente
antes de intentar extenderlo a n=5.

## Por qué esto importa para el próximo paso

Esta es la pieza más básica del aparato (representaciones de S₃ = grupo de permutaciones de 3
partículas, para el caso n=4 donde S₄/(ℤ₂×ℤ₂) = S₃). Para atacar la extensión a n=5 partículas hace
falta repetir esta construcción pero con el grupo de permutaciones relevante para 5 partículas (la
estructura exacta del coset — análoga a S₄/(ℤ₂×ℤ₂)=S₃ para n=4 — todavía no se derivó para n=5; es
el siguiente paso concreto, no trivial, de la extensión).

## Próximo paso concreto (no hecho todavía)

1. Derivar la estructura de coset relevante para n=5 partículas (equivalente a la sección 2.3 del
   paper, "Permutations: ℤ₂×ℤ₂ and S₃", pero para el grupo de simetría de 5 partículas idénticas y
   sus invariantes de Mandelstam independientes — a diferencia de n=4, donde solo hay 2 invariantes
   independientes (s,t, con u=-s-t), para n=5 hay 5 invariantes de Mandelstam independientes con
   relaciones más complejas).
2. Identificar las representaciones irreducibles del grupo de permutaciones relevante (S₅, o el
   subgrupo/cociente que corresponda) y repetir el conteo de generadores del "módulo local" (sección
   2.4 del paper) para ese caso.
3. Recién ahí, imponer la condición de Regge growth (CRG, asumida como axioma) para ver qué
   S-matrices de 5 partículas sobreviven — el resultado nuevo, si lo hay, sería comparar contra lo
   que predicen Einstein/tipo II/heterótico para 5 gravitones.

**Honestidad**: este es un paso de trabajo genuino, no una tarde — el paso de n=4 (2 invariantes
independientes, grupo S₃) a n=5 (5 invariantes de Mandelstam independientes con relaciones lineales
no triviales entre ellos) aumenta la complejidad combinatoria significativamente. Es exactamente el
tipo de escalón que en Yang-Mills nos tomó varias sesiones (ver `yang-mills-project`, docs 22-25).

## Cómo correr la verificación

```
cd fase2-objetivo-tecnico/scripts
python s3_partition_functions.py       # verificación original (2 métodos)
python verificacion_rigurosa.py        # súper-verificación (4 ángulos independientes)
```

Sin dependencias más allá de `sympy` (ya en el entorno).
