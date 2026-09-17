# Fase 2 — Paso 14: generalización de F^i:F^j a los 10 pares, corrección del "candidato más
# simple" propuesto al cierre de la sesión anterior, y primer candidato S5-invariante multilineal
# construido y verificado

## Parte A: generalización de F^i:F^j a los 10 pares (i,j) — hecho, sin sorpresas

Se generalizó la forma cerrada de `F¹:F²|∥` (derivada en `17-forma-cerrada-FdotF.md` solo para el
par (1,2)) a los 10 pares posibles de n=5, **sin asumir el patrón** — verificado explícitamente
para cada par, tal como pedía el punto pendiente de la sesión anterior.

**Convención usada** (generalización directa y sistemática de la ya validada para (1,2)): para el
par (i,j), sea `R = sorted({1..5} - {i,j}) = [A, B, leftover]` (orden ascendente); ambas partículas
usan pivote = la otra del par, companeras = (A,B) en ese orden. Para (i,j)=(1,2): R=[3,4,5] →
A=3,B=4,leftover=5 — coincide exactamente con la convención de las partes 9-17.

**Fórmula generalizada** (mismo patrón, índices sustituidos):

```
F^i:F^j|∥ = −(α_i α_j + α_i' α_j') − Q_ij·[ α_i α_j'/(2 d_iA d_jB) + α_i' α_j/(2 d_iB d_jA) ]
Q_ij = d_ij² + d_ij d_iA + d_ij d_iB + d_ij d_jA + d_ij d_jB + d_iA d_jB + d_iB d_jA
```

**Verificación**: script `scripts/FdotF_todos_los_pares.py`. Método: cinemática numérica explícita
D=6 (5 momentos nulos con conservación de momento, igual que `kinematics_numeric_check.py`),
polarizaciones εᵢ *aleatorias* (con componente transversal genuina, no solo la parte paralela) para
que la extracción de (αᵢ,αᵢ') sea una prueba real y no un caso degenerado. Para cada uno de los 10
pares, en 3 semillas distintas: (a) se verificó explícitamente `Q_i == Q_j` (la misma cantidad
calculada independientemente desde cada partícula del par coincide — chequeo no trivial, no
asumido), (b) se reconstruyó `εᵢ^∥_puro` con los w's (sin término de gauge) y se comparó
`F^i:F^j` calculado directo por vectores contra la fórmula cerrada. **Resultado: coincide en los 10
pares × 3 semillas = 30 casos, con el peor error relativo 3.6×10⁻¹⁵** (precisión de máquina). No se
encontró ninguna asimetría inesperada — el patrón generaliza limpiamente.

## Parte B: corrección de un punto estructural no detectado al cierre de la sesión anterior

El "candidato más simple" propuesto en `17-forma-cerrada-FdotF.md` para probar primero,
`Σ_{i<j} F^i:F^j`, **no puede ser una S-matriz de 5 fotones**: cada término de esa suma es lineal en
εᵢ,εⱼ pero de **grado cero** en las otras 3 polarizaciones (no las contiene en absoluto). Una
amplitud de n fotones debe ser **exactamente multilineal** — grado 1 en cada una de las n
polarizaciones — porque así entra cada pata externa de fotón (así están construidos todos los
building blocks del caso n=4 en el paper: para n=4, con 4 patas, los pares (εᵢ.εⱼ)(εₖ.εₗ) cubren
las 4 polarizaciones exactamente una vez cada una, porque 4 es par y se puede apantallar
completamente en parejas). **Para n=5 (impar) esto no se puede hacer solo con parejas** — hace
falta un objeto adicional lineal en una sola polarización (el "singlete"), contraído solo con
momentos.

**Building block para el singlete**: `C^k_{ab} := p_a.F^k.p_b` (ya estudiado en las partes 12-16;
lineal en εₖ solamente, no involucra εₐ,ε_b).

## Parte C: primer candidato multilineal — hallazgo real de cancelación exacta, corregido

**Ansatz inicial** (grado (1,1,1,1,1) en las 5 polarizaciones, la estructura mínima que cubre las 5
patas: 2 parejas vía F:F + 1 singlete vía C):

```
T(i,j,l,m,k) := F^i:F^j · F^l:F^m · C^k_{ij}
S_naive := Σ_{σ∈S5} T(σ(1),σ(2),σ(3),σ(4),σ(5))     (suma sobre la órbita completa — operador
                                                        de Reynolds, S5-invariante por construcción)
```

**Resultado: S_naive ≡ 0 idénticamente** (verificado numéricamente en varias semillas: valores
∼10⁻¹² a 10⁻¹⁰, ruido de punto flotante de un cero exacto — no un "casi cero" dependiente de gauge).

**Diagnóstico** (no se aceptó el cero sin explicarlo): `C^k_{ab}` es **antisimétrico** en sus 2
índices de momento (`C^k_{ab} = −C^k_{ba}`, verificado explícitamente:
`C(1;2,3)=14.049…, C(1;3,2)=−14.049…`), mientras que `F^i:F^j` es **simétrico** en (i,j). El
producto `F^i:F^j · C^k_{ij}` es entonces **impar bajo el intercambio i↔j** — y esa transposición
está incluida en la suma sobre la órbita completa de S5. Sumar sobre ambos órdenes (i,j) y (j,i)
cancela exactamente esa contribución, para cualquier partición — de ahí el cero idéntico. Es un
hallazgo estructural real (el tipo de "resultado demasiado bueno/limpio que hay que desconfiar",
ya varias veces en este proyecto), no un bug de cómputo.

**Corrección**: multiplicar por una cantidad también antisimétrica bajo i↔j construida con los
invariantes de Mandelstam ya calculados (partes 3-8), de modo que el producto completo sea par:

```
T(i,j,l,m,k) := F^i:F^j · F^l:F^m · C^k_{ij} · (d_ik − d_jk)
S := Σ_{σ∈S5} T(σ(1),…,σ(5))
```

**Verificado** (`scripts/primer_candidato_s5_invariante.py`):
1. **Invariancia de gauge**: shift simultáneo aleatorio εᵢ→εᵢ+ζᵢpᵢ en las 5 polarizaciones deja S
   invariante (diferencia relativa ∼2×10⁻¹⁵).
2. **Invariancia S5**: recalculado tras relabeling aleatorio de las 5 partículas — coincide exacto
   (diferencia relativa ∼1×10⁻¹⁶; esto es garantizado por construcción — suma sobre la órbita
   completa — pero se verificó igual para atrapar bugs de implementación, no para "descubrir" la
   propiedad).
3. **No trivialidad**: S ≠ 0 en 5 semillas distintas, con valores de orden 10³–10⁶ (sin cancelación
   patológica).

## Estado y evaluación honesta

Se tiene ahora un primer candidato genuino, no trivial, multilineal en las 5 polarizaciones,
gauge-invariante y S5-invariante — construido enteramente a partir de building blocks ya
verificados (F^i:F^j de las partes 16-17, C^k_ab de la parte 16, y el anillo de Mandelstam de las
partes 3-8). **Esto todavía no confirma que sea una S-matriz físicamente aceptable**: falta (a)
determinar su grado total y compararlo contra la serie de Hilbert/Molien esperada para el caso con
spin (no calculada aún — la ya computada es solo para escalares), (b) verificar que no sea
idénticamente proporcional a algo ya conocido/trivial de menor grado, (c) el paso final del
proyecto: verificar la condición de Regge growth.

## Qué falta (próximo paso concreto)

1. Determinar si S (grado (1,1,1,1,1) en polarizaciones, y algún grado homogéneo en momentos —
   contar: F:F es grado 2 en momentos, C es grado 2 en momentos, el factor (d_ik−d_jk) es grado 1
   en momentos ⇒ S tiene grado 2+2+2+1=7 en momentos) es el candidato de **grado mínimo** o si hay
   uno más simple con menos factores de Mandelstam extra (el factor (d_ik−d_jk) fue la corrección
   mínima que se probó, pero puede no ser la única ni la más económica — explorar alternativas,
   p.ej. multiplicar por otras combinaciones antisimétricas de grado más bajo si existen).
2. Enumerar sistemáticamente TODOS los candidatos de grado mínimo multilineales en las 5
   polarizaciones (no quedarse con el primero que sobrevive) — construir la base completa, análoga
   a como el paper original arma la base completa de invariantes para n=4.
3. Verificar la condición de Regge growth (objetivo final) sobre los candidatos construidos.

## Archivos

- `scripts/FdotF_todos_los_pares.py` — generalización y verificación de F^i:F^j a los 10 pares.
- `scripts/primer_candidato_s5_invariante.py` — construcción y verificación del primer candidato
  S5-invariante multilineal (incluye la demostración numérica de la cancelación de la versión
  ingenua y la versión corregida).
