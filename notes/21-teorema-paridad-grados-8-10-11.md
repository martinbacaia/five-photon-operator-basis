# Fase 2 — Paso 17: teorema de paridad (grados 8 y 10 = 0, sin cómputo) y extensión de la tabla
# numérica de la serie de Hilbert con spin a grado 11

## Teorema de paridad: TODO invariante multilineal de n=5 fotones tiene grado TOTAL impar en
## momentos — grados 8 y 10 son idénticamente cero, y esto se puede probar, no solo observar

**Motivación**: en vez de lanzar cómputo a ciegas en grados 8 y 10, conviene preguntarse primero si
son siquiera alcanzables con los átomos disponibles.

**Argumento** (generaliza la observación de la parte 19 sobre por qué el corrector mínimo es grado
2, no grado 1): la invariancia de gauge fuerza toda la dependencia de cada εᵢ a entrar
exclusivamente a través de `F^i_{μν}` (hecho estructural ya usado en toda la sesión, desde el
Camino B). Un "átomo" (contracción de Lorentz irreducible) que es lineal en exactamente `m`
polarizaciones distintas se construye contrayendo `m` field-strengths `F^{i_1},...,F^{i_m}` entre
sí y con momentos externos hasta saturar sus `2m` índices libres. Cada contracción **F-F** (entre 2
field-strengths) usa 2 índices y **0** momentos extra; cada contracción **F-momento** usa 1 índice y
**1** momento extra. Si hay `x` contracciones F-F e `y` contracciones F-momento, entonces
`2x + y = 2m` (saturación de índices), luego `y = 2(m − x)` — **siempre par**. El grado total del
átomo es `m` (momentos intrínsecos de cada F, uno por F) `+ y` (momentos externos) `= m + par`, que
tiene **la misma paridad que `m`**.

**Consecuencia**: cualquier invariante completo, construido como producto de átomos que cubren las
5 polarizaciones exactamente una vez cada una (más cualquier multiplicador de Mandelstam, que
siempre tiene grado par porque se construye con productos de `d_ab = p_a.p_b`), tiene grado total
`≡ Σm_atomo (mod 2) = 5 (mod 2) = 1` — **siempre impar**, sin excepción, para CUALQUIER combinación
de átomos y correctores.

**Verificado por conteo explícito** (`scripts` inline): las 3 topologías mínimas posibles
(2 parejas+1 singlete → 7; 1 pareja+3 singletes → 11; 5 singletes → 15) dan efectivamente grado
base impar en los 3 casos, consistente con el teorema.

**Consecuencia práctica**: **grados 8 y 10 tienen dimensión 0, por imposibilidad estructural** — no
existe ninguna forma de construir un invariante multilineal de grado par para n=5. No hace falta
computar nada: no hay ningún candidato que evaluar. (Chequeo de consistencia: para n=4, la misma
cuenta da `Σm_atomo = 4 (mod 2) = 0` — grado siempre PAR, consistente con que las construcciones
conocidas del paper original para n=4 —pares (εᵢ.εⱼ)(εₖ.εₗ)— tienen grado 4, 6, 8... todas pares.
Esto es un chequeo de sanidad no trivial: la teoría predice una diferencia estructural real entre
n par y n impar, y coincide con lo que ya se sabía para n=4.)

## Grado 11: extensión de la tabla numérica

Grado 11 es el primer grado impar disponible después de 9 (7 y 9 ya resueltos: 7→dimensión 0
—parte 19—, 9→dimensión 1 —parte 20—). Se probaron 2 familias de construcción
(`scripts/hilbert_series_spin_grado11.py`):
- **(A)** "2 parejas F:F + 1 singlete C" (grado base 7) × un multiplicador de Mandelstam de grado 4
  (55 monomios posibles, producto de 2 de los 10 productos punto entre los 5 slots) × 3 formas de
  contraer el singlete → 165 semillas.
- **(B)** "1 pareja F:F + 3 singletes C" (grado 11 exacto, sin corrector) — 2 variantes de
  contracción representativas (no exhaustivo).

### Problema encontrado: el rango numérico NO convergió con pocas muestras — diagnosticado antes
### de reportar cualquier número

A diferencia de grado 9 (donde 25 muestras ya daban un salto de 14 órdenes de magnitud, sin
ambigüedad), en grado 11 el rango estimado **cambió con el número de muestras y con qué muestras se
usaban**: 5 (25 muestras, semillas 1-25), 5 (40 muestras, semillas 1-40), 9 (40 muestras, semillas
101-140), 3 (40 muestras, semillas 201-240). Esta inestabilidad **se reportó explícitamente antes de
sacar conclusiones** — no se aceptó ningún número hasta entender la causa.

**Diagnóstico**: no es error de precisión de punto flotante por cancelación catastrófica (los
valores singulares "pequeños" en los primeros intentos, del orden de 10⁻² a 10⁻¹, están muy por
encima del piso de ruido de máquina, ∼10⁻¹⁵ tras la normalización) — es un problema de **muestreo
insuficiente**: con solo 25-40 puntos cinemáticos aleatorios, algunas direcciones genuinamente
independientes del espacio de invariantes tienen una "huella" muy débil en esas muestras concretas
(su varianza a través de esos puntos es chica en comparación con las direcciones dominantes), y por
eso su singular value aparece artificialmente pequeño y se confunde con ruido — hasta que se usan
suficientes muestras para que esas direcciones se manifiesten con claridad.

**Primer intento de resolución (INCORRECTO, corregido más abajo)**: al aumentar a 200 muestras
independientes (float64), apareció un salto limpio de 13 órdenes de magnitud entre la 10ª (0.032) y
la 11ª (3.7×10⁻¹⁵) valor singular, sugiriendo rango 10. **Pero al probar más lotes grandes e
independientes de muestras (300, 150, 250) el resultado NO se sostuvo**: rango 12, 3, y 10
respectivamente, según el lote — es decir, la inestabilidad **empeoró** en vez de resolverse con más
muestras. Esto descartó la hipótesis de "solo hace falta más muestreo" y señaló un problema más
profundo.

## Diagnóstico correcto: NO es un problema de muestreo — es precisión numérica insuficiente en
## float64 para objetos de grado 11 (11 potencias de momento por término, sumados sobre 120
## permutaciones, con cancelaciones internas severas)

Se rehizo el cálculo completo con **aritmética de precisión arbitraria** (`mpmath`, 60 dígitos
decimales — `scripts/hilbert_series_grado11_altaprecision.py`, con generador de cinemática D=6
propio en alta precisión, sin depender de `numpy`/`scipy`) para las mismas 167 semillas.

**Resultado, con 30 muestras cinemáticas independientes**: el rango es **exactamente 3**, con un
residual tras el 3er pivote de la eliminación gaussiana de **1.7×10⁻⁶⁰** — es decir, cero hasta el
límite de la propia precisión de 60 dígitos usada. Sin ambigüedad alguna, sin necesidad de elegir un
umbral de tolerancia arbitrario: el 4º pivote candidato es indistinguible de cero al nivel de
1 parte en 10⁶⁰.

**Conclusión metodológica importante para el resto del proyecto**: la oscilación observada en
float64 (3, 5, 9, 10, 12 según la muestra) fue **enteramente un artefacto de precisión numérica
insuficiente** — no reflejaba ninguna propiedad real del espacio de invariantes ni "direcciones
débiles que necesitan más muestras". A partir de grado ~11 (production de ~11 potencias de momento,
sumadas sobre 120 términos con cancelaciones internas), **float64 deja de ser confiable** para este
tipo de cálculo — hace falta precisión arbitraria. **Se recomienda para el resto del proyecto usar
`mpmath` (o aritmética racional exacta) en cualquier cálculo de grado ≥10**, y tratar con sospecha
cualquier resultado de rango en float64 que no muestre un salto de valores singulares de al menos
~10 órdenes de magnitud, estable frente a variar el lote de muestras.

**Re-chequeo de grado 9 con el mismo método de alta precisión** (por las dudas, dado que fue
calculado originalmente en float64): confirmado — rango exactamente 1, residual tras el pivote
4.3×10⁻⁶⁰. El resultado de la parte 20 se sostiene sin cambios.

## Resultado final de grado 11

**Dimensión = 3**, dentro de la familia de construcción probada (topología "2 parejas + 1 singlete"
× corrector Mandelstam grado 4, en sus 3 formas de contracción, más 2 semillas representativas de
la topología alternativa "1 pareja + 3 singletes"). Mismo alcance honesto que en grado 9: no es
necesariamente exhaustivo sobre todas las topologías posibles en grado 11, pero es un resultado
numérico ahora sólido (no ambiguo) dentro de la familia identificada.

## Tabla acumulada de la serie de Hilbert con spin (dentro de las topologías exploradas)

| grado | dimensión | cómo se sabe |
|-------|-----------|--------------|
| 7     | 0         | teorema general (parte 19): cualquier construcción de esta topología se anula |
| 8     | 0         | imposible por paridad (parte 21) — ningún átomo puede sumar grado par |
| 9     | 1         | cómputo numérico, confirmado en alta precisión (partes 20 y 21) |
| 10    | 0         | imposible por paridad (parte 21) |
| 11    | 3         | cómputo en alta precisión (parte 21) — dimensión saltó de 1 a 3 |

## Archivos

- `scripts/hilbert_series_spin_grado11.py` — versión float64 (histórica; sirvió para descubrir el
  problema de precisión, no usar sus números de rango como definitivos).
- `scripts/hilbert_series_grado11_altaprecision.py` — versión de alta precisión (mpmath, 60
  dígitos), la que da el resultado confiable: rango 3 en grado 11, rango 1 en grado 9 (re-chequeo).
