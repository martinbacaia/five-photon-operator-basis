# Fase 2 — Paso 16: primer punto de la "serie de Hilbert con spin" — cálculo numérico riguroso
# (no fórmula cerrada) — dimensión exacta = 1 en grado 9, resolviendo la pregunta abierta de la
# parte 19

## Por qué numérico y no una fórmula cerrada

El propio paper de referencia (Henning-Lu-Melia-Murayama 2017) admite que para el caso con spin en
n=5 no se sabe siquiera si el anillo de invariantes es Cohen-Macaulay — la maquinaria de Molien
usada para el caso escalar (partes 3-8) no se generaliza de forma directa cuando los datos de
polarización (αᵢ,αᵢ') no forman una representación lineal simple de S₅ sino una estructura de
fibrado (parte 14-15). Intentar una fórmula cerrada tipo Molien en una sesión no es realista. En su
lugar, se calculó **numéricamente** la dimensión real del espacio de invariantes en un grado
concreto — la misma técnica de "operador de Reynolds + rango numérico" ya usada con éxito en las
partes 4-6 para el anillo de Mandelstam puro, aplicada ahora a construcciones que sí incluyen
polarización.

## Método

Se armó una familia de 30 "semillas" de grado 9 (la topología ya establecida en la parte 19: 2
parejas `F^i:F^j` + 1 singlete `C^k_{ab}` + 1 factor de Mandelstam de grado 2), variando
sistemáticamente:
- **3 formas de contraer el singlete**: con la primera pareja (`C^k_{ij}`), con la segunda
  (`C^k_{lm}`), o de forma mixta (`C^k_{il}`).
- **los 10 productos punto posibles** entre los 5 "slots" (i,j,l,m,k) como factor multiplicador de
  grado 2.

Cada semilla se simetrizó sumando sobre la **órbita completa de S₅** (exactamente como en las
partes 18-19), y se evaluó el resultado en **25 cinemáticas numéricas independientes** (D=6,
momentos nulos aleatorios con conservación de momento resuelta, mismo generador que todos los
scripts anteriores). Se armó la matriz (25 muestras × 30 semillas) y se calculó su **rango
numérico** — sin asumir independencia lineal de nada, calculándola directamente.

## Resultado

- **12 de las 30 semillas se anulan idénticamente** (consistente con el mecanismo ya identificado en
  la parte 19: un multiplicador que es *simétrico* bajo la transposición que ya mata a la
  construcción base — p.ej. `d(ij)` o `d(lm)` para la contracción `C_P` — no compensa la
  antisimetría y el término sigue cancelándose).
- De las 18 semillas no nulas, el **rango numérico es exactamente 1** — con una separación de
  valores singulares completamente inequívoca: `4.24, 2.7×10⁻¹⁴, 2.0×10⁻¹⁵, …` (un salto de 14
  órdenes de magnitud entre el primer valor singular y el resto, que son ruido numérico puro).

**Conclusión: dentro de esta familia de construcción (la topología "2 parejas + 1 singlete", con
cualquier elección de contracción del singlete y cualquier corrector de Mandelstam de grado 2), el
espacio de invariantes de grado 9 tiene dimensión exactamente 1** — no hay una familia más rica
escondida entre estas variantes. Esto responde la pregunta que había quedado abierta en
`19-enumeracion-grado-minimo.md` ("¿dimensión 1 o más?"): dentro de este sector, es 1. (El candidato
`T_grado9_v1` ya encontrado y verificado en la parte 18 es, hasta normalización, ESE único
invariante — es una combinación lineal de 2 de las 18 semillas no nulas de esta lista,
`C_P__d(ik)` y `C_P__d(jk)`.)

## Alcance honesto de este resultado — qué NO dice

1. **No es una serie de Hilbert completa.** Solo cubre la topología "2 parejas F:F + 1 singlete C +
   1 corrector Mandelstam grado 2" en grado 9. No incluye otras topologías del mismo grado que
   podrían existir (p.ej. usando el objeto "sandwich" `p_a.F^i.F^j.p_b`, de grado 4, en vez de una
   de las parejas `F^i:F^j` — combinación `1 F:F(2) + 1 sandwich(4) + 1 C(3) = 9`, no probada
   todavía).
2. **No prueba que grado 9 sea el grado mínimo absoluto** — solo que dentro de la topología ya
   identificada como la más económica conocida, es el piso (parte 19), y que en ESE piso la
   dimensión es 1.
3. Es un cálculo numérico (rango de una matriz de valores flotantes), no una prueba simbólica
   exacta — aunque la separación de valores singulares (14 órdenes de magnitud) hace la conclusión
   extremadamente robusta, no es una demostración formal de independencia lineal exacta como sí se
   pudo hacer con álgebra simbólica en el caso escalar (partes 3-8, verificado además con Singular).

## Qué falta (próximo paso concreto)

1. Extender la tabla a más grados (7 ya sabido =0; probar 8, 10, 11) para empezar a armar una
   verdadera tabla tipo "serie de Hilbert con spin", aunque sea numérica y restringida a las
   topologías conocidas — el objetivo sería ver si emerge un patrón reconocible (como pasó con la
   serie de Molien escalar).
2. Incluir la topología del objeto "sandwich" `p_a.F^i.F^j.p_b` en grado 9 para saber si agrega una
   dimensión extra al espacio ya encontrado (dimensión 1) o si termina siendo proporcional también
   al mismo invariante único.
3. Con el único invariante de grado 9 ya confirmado como la pieza más simple disponible, sería
   razonable adelantar el chequeo de Regge growth sobre él en paralelo a seguir extendiendo la
   tabla — no hace falta esperar la clasificación completa para hacer el primer chequeo del
   objetivo final del proyecto.

## Archivos

- `scripts/hilbert_series_spin_grado9.py` — construcción de las 30 semillas, simetrización,
  evaluación numérica en 25 cinemáticas, cálculo de rango y valores singulares.
