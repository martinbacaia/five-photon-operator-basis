# Fase 2 — Paso 13: forma cerrada de F^i:F^j (la parte "∥"), verificada

## Resultado

Continuando el Camino B (`16-camino-B-contracciones-directas.md`), se derivó la forma cerrada de la
parte "paralela" de `F^i:F^j := F^i_μν F^j^μν` — el building block gauge-invariante que combina las
polarizaciones de 2 partículas distintas i≠j. Para i=1,j=2 (con la convención: partícula 1 usa
pivote=2,compañeras=(3,4); partícula 2 usa pivote=1,compañeras=(3,4) — misma receta que en las
partes 9-11, aplicada independientemente a cada partícula):

```
F¹:F² = 2(p1.p2)(ε1⊥.ε2⊥)  +  F¹:F²|∥

F¹:F²|∥ = −(α₁α₂ + α₁'α₂')  −  Q·[ α₁α₂'/(2 d₁₃d₂₄)  +  α₁'α₂/(2 d₁₄d₂₃) ]
```

con `Q = d12²+d12d13+d12d14+d12d23+d12d24+d13d24+d14d23` — **la MISMA Q** que ya apareció en la
matriz M de la fórmula inversa (parte 13) y en las relaciones lineales entre los C_jk (parte 16).

**Lo notable**: los coeficientes de α₁α₂ y α₁'α₂' son exactamente **−1, sin ninguna dependencia
cinemática** — solo los términos "cruzados" (α₁α₂' y α₁'α₂) llevan el factor Q. Esto es una
estructura mucho más simple de lo que se podría haber esperado.

## Verificación

1. **Cancelación de gauge** (simbólica, exacta): se verificó que los coeficientes de a₁, a₂ y
   a₁·a₂ (los parámetros puros de gauge de cada partícula) son idénticamente cero en la expresión
   completa antes de fijarlos a cero — confirma que la fórmula es genuinamente gauge invariante en
   ambas polarizaciones simultáneamente, no por casualidad.
2. **Verificación numérica de la parte paralela sola**: se construyeron ε₁,ε₂ *puros* (sin
   componente transversal, ε=ε^∥ con a=0) con la cinemática explícita D=6 ya usada en sesiones
   anteriores, se calculó `F¹:F²` directo por vectores, y se comparó con la fórmula cerrada —
   **coincide exacto** (diferencia 0.0).
3. **Verificación de la descomposición completa, incluyendo ε1⊥.ε2⊥ construido de forma
   independiente** (no solo inferido por resta): se construyeron explícitamente los vectores
   ε1⊥ := ε1 − ε1^∥_puro − a₁p₁ y ε2⊥ := ε2 − ε2^∥_puro − a₂p₂ (con a₁,a₂ obtenidos comparando
   ε_i.p_m contra la parte paralela reconstruida, mismo método que en partes anteriores), y se
   comprobó — **chequeo no trivial, encontró y corrigió un descuido propio en el primer intento**
   (olvidar restar el término a_i·p_i antes de verificar ortogonalidad, lo que daba productos
   ε_i⊥.p_m ≠ 0 para m≠i) — que tras la corrección, **ambos ε1⊥ y ε2⊥ son exactamente ortogonales
   a los 5 momentos p₁..p₅** (∼1e-15), confirmando que viven genuinamente en el espacio
   transversal. Con esto, `ε1⊥.ε2⊥` calculado directamente **coincide exacto** (diferencia ~1e-14)
   con el valor inferido por resta `(F¹:F² − F¹:F²|∥)/(2·p1.p2)` — cierra por completo la
   descomposición `F¹:F² = 2(p1.p2)(ε1⊥.ε2⊥) + F¹:F²|∥` con todas sus piezas verificadas
   independientemente.

## Qué falta

1. Repetir la derivación para los otros 9 pares (i,j) — por simetría de la construcción debería
   dar la misma forma con los índices permutados, pero **hay que verificarlo explícitamente**, ya
   que la convención de pivote usada por cada partícula podría introducir asimetrías (algo que la
   sesión anterior ya mostró que puede pasar de forma no obvia).
3. Con las 10 fórmulas F^i:F^j (i<j) más el anillo de invariantes de Mandelstam (partes 3-8), armar
   productos y sumas S₅-invariantes — candidatos concretos a S-matrices de fotones de n=5. El
   candidato más simple a probar primero: la suma simetrizada `Σ_{i<j} F^i:F^j` (mínima invariante
   de S₅ construible a partir de estos building blocks), y verificar qué tan lejos está de ser el
   S-matrix de mínimo grado derivada esperado.
4. Verificar Regge growth sobre los candidatos que se vayan construyendo (objetivo final).

## Archivos

- Cálculos hechos inline (extendiendo `kinematics_numeric_check.py`); no se generó script nuevo
  separado en esta sesión — ver el historial de comandos para reproducir.
