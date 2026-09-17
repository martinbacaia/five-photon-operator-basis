# Fase 2 — Paso 25: el objeto "sandwich" revela que grado 9 tiene dimensión ≥2 (corrección de la
# parte 20) — y la cancelación dentro de grado 9 baja el crecimiento de s⁴ a s^3.5

## Corrección importante: la serie de Hilbert de grado 9 calculada en la parte 20 estaba incompleta

La parte 20 concluyó "dimensión 1" en grado 9, pero esa búsqueda **solo cubrió la topología
"2 parejas F:F + 1 singlete C"** — nunca se probó el objeto "sandwich"
`M^{ij}_{ab} := p_a.F^i.F^j.p_b` (mencionado como pendiente desde la parte 16, nunca implementado
hasta ahora). Este objeto es una contracción de Lorentz distinta a `F^i:F^j` (una es la "traza
completa" de los 2 field-strengths, la otra es un "producto de matrices" sandwich entre 2 momentos
externos) — bilineal en `ε_i,ε_j` igual que `F^i:F^j`, pero de grado 4 en momentos (no 2).

**Se implementó y verificó gauge-invariante** (shift `εᵢ→εᵢ+ζpᵢ` deja `M^{ij}_{ab}` invariante,
diferencia ~1e-15):

```
M^{ij}_{ab} = [(p_a.p_i)(ε_i.p_j) - (p_a.ε_i)(p_i.p_j)]·(ε_j.p_b)
              - [(p_a.p_i)(ε_i.ε_j) - (p_a.ε_i)(p_i.ε_j)]·(p_j.p_b)
```

**Topología de grado 9 usando el sandwich**: `M^{ij}_{kl} · F^l:F^m · C^k_{lm}` (grado 4+2+3=9),
simetrizada sobre la órbita de S₅. La elección de contracción del sandwich importa (igual mecanismo
de antisimetría de siempre): `M^{ij}_{lm}` (ambos argumentos externos en la misma pareja del F:F) se
anula idénticamente; `M^{ij}_{kl}` (uno de los argumentos es el propio singlete k) **sobrevive**,
verificado no-trivial en varias semillas.

## Es un invariante genuinamente NUEVO, no proporcional al ya conocido

Se comparó `T_sandwich` contra el `T_grado9` ya conocido (parte 18) en 5 semillas de cinemática
genérica: el cociente `T_grado9/T_sandwich` da `-8.36, 34.9, 14.4, -7.1, 11.7` — **NO es constante**,
confirmando que son 2 invariantes **linealmente independientes**. **Corrección**: el espacio de
invariantes de grado 9 tiene dimensión **al menos 2**, no 1 como se concluyó en la parte 20 (esa
conclusión sigue siendo válida DENTRO de la familia limitada que se probó ahí, pero no es la
dimensión real del espacio completo).

## Consecuencia inmediata: SÍ hay margen para cancelar dentro de grado 9

Con dimensión ≥2, la corrección conceptual de la parte 24 ("solo se puede cancelar DENTRO del mismo
grado, y hace falta dimensión >1 para eso") ahora aplica directamente a grado 9. Se evaluaron ambos
invariantes en el límite de Regge de n=5 ya validado (parte 23):

```
T_grado9   ~ E^(7.9999)  ~ s^4.0000
T_sandwich ~ E^(8.0001)  ~ s^4.0000
```

Mismo orden líder — exactamente el escenario donde la cancelación es posible. Se extrajo la razón de
cancelación por extrapolación de Richardson (usando pares consecutivos de `E`, convergiendo
establemente): `r_∞ = -19.742974` (4+ cifras estables). La combinación
`T_grado9 - r_∞·T_sandwich`:

```
E=    1000  combo=-5.38e+22
E=    3000  combo=-1.14e+26
E=   10000  combo=-5.14e+29
E=   30000  combo=-1.12e+33
E=  100000  combo=-5.13e+36

pendientes consecutivas: 6.968, 6.990, 6.999, 7.000  ->  E^7 EXACTO
```

**El término líder `E^8~s^4` se cancela por completo, y lo que queda crece como `E^7 ~ s^3.5`** —
una mejora real y limpia (pendiente convergiendo a 7.0000 exacto), no un artefacto numérico.

## Estado y qué significa

- `s^3.5` sigue **excediendo** el bound `s²`, así que grado 9 (con estos 2 invariantes conocidos)
  **todavía no da un candidato válido** — pero el hallazgo es significativo: la combinatoria de
  grado 9 tiene más estructura de la que se había mapeado, y hay margen de mejora real dentro del
  grado (no hacía falta ir a grados más altos).
- **Corrección adicional a anticipar**: si grado 9 tenía una dimensión oculta por no incluir la
  topología "sandwich", es muy probable que **grado 11 (parte 21, "dimensión 3") también esté
  incompleto** por la misma razón — no se probó el sandwich ahí tampoco. Recalcular grado 11
  incluyendo el sandwich es ahora una prioridad, no una curiosidad.
- Pregunta abierta inmediata: ¿hay una TERCERA dimensión en grado 9 (otra variante del sandwich, u
  otra topología) que permita cancelar también el remanente `s^3.5` y bajar más, quizás hasta `s²`?
  Se probó una sola variante alternativa del sandwich (`M^{ij}_{km}`) y resultó ser exactamente
  `-M^{ij}_{kl}` (no una dirección nueva) — pero no se agotó la búsqueda de otras topologías
  (p.ej. sandwich con distinta asignación de índices, o construcciones de 3 field-strengths tipo
  traza ya descartadas en grado bajo pero quizá viables acá con un corrector).

## Próximo paso concreto

1. Recalcular la dimensión real de grado 9 y grado 11 incluyendo sistemáticamente la topología
   sandwich (no solo la variante encontrada) — repetir el cálculo tipo "serie de Hilbert numérica"
   de las partes 20-21 con el sandwich agregado al conjunto de semillas.
2. Si aparece una 3ª dimensión en grado 9, repetir el procedimiento de cancelación (ahora con 3
   parámetros libres en vez de 1) y ver si se puede bajar por debajo de `s^3.5`, idealmente hasta
   `s²`.
3. Mantener la honestidad: el objetivo es no quedarse en "encontramos una mejora" sino perseguir
   hasta el final si existe o no un candidato que satisfaga el bound.

## Archivos

- Cálculos hechos inline (extendiendo `regge_limit_n5.py` con la función `Msandwich` y las nuevas
  variantes de contracción); no se generó un script nuevo separado dado que son extensiones cortas.
