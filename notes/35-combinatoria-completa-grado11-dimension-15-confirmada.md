# Fase 2 — Paso 35: combinatoria COMPLETA de las familias C y D de grado 11 — dimensión 15
# confirmada (no era un artefacto de los subconjuntos representativos usados en la parte 28)

## Motivación

La parte 28 encontró dimensión ≥15 para grado 11 usando las familias A+B+C+D, pero las familias C y
D se habían restringido a subconjuntos representativos (144 y 27 semillas respectivamente) "para no
explotar combinatoriamente" — dejando abierta la pregunta de si la combinatoria completa revelaría
más direcciones independientes. Se decidió (recomendación propia, aprobada por el usuario) ampliar la
búsqueda a la combinatoria COMPLETA antes de seguir con otras líneas.

**Combinatoria completa**: familia C, 6 (sandwich) × 6 (singlete) × 10 (corrector, TODOS los pares,
antes 4) = 360 semillas (antes 144); familia D, 6×6×10 (TODOS los pares para el singlete, antes 3) =
360 semillas (antes 27). Total A+B+C+D = 887 semillas (antes 338).

## Problema de infraestructura encontrado y resuelto (importante, no cosmético)

Un primer intento de correr esto con `run_in_background` (proceso totalmente desatendido) resultó
POCO CONFIABLE en este entorno: el proceso quedaba con uso de CPU casi nulo durante largos períodos
(verificado con `Get-Process | Select CPU`, no solo "se ve lento") sin haber terminado ni fallado —
un cuelgue de facto, aparentemente por scheduling del sistema operativo con procesos backgroundeados
largos, NO un bug del código (la MISMA computación, corrida en foreground con un timeout explícito
generoso, o en llamadas cortas repetidas, terminó siempre limpia, sin colgarse, con tiempos por
muestra de 25-135s). **Lección para el resto del proyecto**: para cómputos largos en este entorno,
preferir llamadas foreground con timeout explícito (o dividir en varias llamadas cortas con guardado
incremental) en vez de `run_in_background` desatendido.

Se implementó `grado11_combinatoria_completa_incremental.py` (guarda resultados en un pickle,
resumible entre llamadas) para poder acumular las 20+ muestras necesarias en varias llamadas cortas.

## Resultado: dimensión 15 confirmada, EXACTA, en 2 lotes independientes

- **Lote 1** (semillas de cinemática 1-20): 656 de 887 semillas no nulas, **rango = 15**, residual
  tras el último pivote `4.4×10⁻⁶¹` (esencialmente exacto).
- **Lote 2** (semillas de cinemática 201-220, completamente independiente): también 656 no nulas,
  **rango = 15**, residual `1.4×10⁻⁶¹`.

**La combinatoria completa de C y D NO encontró ninguna dirección nueva** más allá de las 15 ya
encontradas con los subconjuntos representativos de la parte 28 — el rango se mantiene exactamente
en 15 pese a más que duplicar el número de semillas (338→887). Esto es una confirmación fuerte (no
solo "no se buscó lo suficiente antes") de que **la dimensión real de grado 11, dentro de las 4
topologías ya identificadas (A: 2×F:F+singlete; B: F:F+3 singletes; C: F:F+sandwich+singlete; D: 2
sandwiches+singlete), es exactamente 15** — no ≥15 como una cota provisional, sino 15 como resultado
firme dentro de este espacio de construcciones.

## Qué significa esto para el proyecto

- El resultado de las partes 30/32/33/34 (base de 9/2 vectores universales, evaluados en un espacio
  de "dimensión ≥15") queda sobre un terreno más firme: 15 es la dimensión real de las topologías
  conocidas, no una subestimación por muestreo incompleto.
- Para encontrar generadores adicionales (si los hay) haría falta una topología GENUINAMENTE NUEVA,
  no cubierta por A, B, C, D (p.ej. estructuras con 3 sandwiches, o contracciones con el
  antisimétrico de Levi-Civita si D lo permite, u otras combinaciones de field-strengths no
  exploradas) — no simplemente ampliar la combinatoria de contracciones dentro de las 4 topologías
  ya conocidas, que ya se demostró agotada.

## Archivos

- `scripts/grado11_combinatoria_completa.py` — construye las 887 semillas (combinatoria completa).
- `scripts/grado11_combinatoria_completa_incremental.py` — versión resumible (guarda en pickle,
  pensada para correr en varias llamadas foreground cortas en vez de background desatendido).
- `scripts/combinatoria_completa_rows.pkl`, `combinatoria_completa_rows_batch2.pkl` — los 2 lotes
  independientes (20 muestras cada uno), rango 15 en ambos.
