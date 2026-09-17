# Fase 2 — Paso 7: estructura COMPLETA del anillo de invariantes de S₅, calculada con Singular

## Resultado (el más fuerte de todo el proyecto hasta ahora)

Usando `finvar.lib` de Singular (algoritmo de Decker-Heydtmann-Schreyer para la normalización de
Noether + Kemper-Steel para invariantes secundarios), sobre las mismas matrices generadoras de S₅
ya derivadas y verificadas en Python (`n5_mandelstam_representation.py`), se obtuvo la estructura
**completa** del anillo de invariantes:

- **5 invariantes primarios** (una normalización de Noether — siempre hay exactamente
  dim(V)=5 de ellos, sea o no libre el anillo), de grados **2, 3, 4, 5, 6**.
- **6 invariantes secundarios**, de grados **0, 6, 7, 8, 9, 15** (el de grado 0 es simplemente la
  constante 1, siempre presente). Complejidad creciente: 1, 158, 283, 495, 638 y 3565 términos
  respectivamente — el de grado 15 es un polinomio genuinamente enorme.

Esto da la "descomposición de Hironaka" completa del anillo:
R^{S₅} = ⊕ᵢ ηᵢ · ℂ[θ₁,θ₂,θ₃,θ₄,θ₅], con θ (primarios) de grados 2,3,4,5,6 y η (secundarios) de
grados 0,6,7,8,9,15.

**Esto coincide y extiende exactamente lo encontrado a mano** en los pasos 4 y 6
(`04-generadores-explicitos-s5.md`, `06-generadores-grado5-6.md`): los generadores de grados 2, 3,
4 y 5 que fuimos encontrando uno por uno con el operador de Reynolds son, en efecto, 4 de los 5
invariantes primarios; el quinto primario (grado 6) y los 6 secundarios son lo que hacía falta para
cerrar la estructura completa — exactamente donde nuestro método artesanal empezaba a volverse
lento y propenso a errores (grado 6, donde encontramos el 3er bug de la sesión).

## Verificación (no se aceptó el resultado sin cruzarlo)

Se reconstruyó la serie de Hilbert (dimensión de invariantes por grado) a partir de la estructura
de Singular, con la fórmula estándar:

H(x) = (Σᵢ x^{deg ηᵢ}) / (Πⱼ (1 − x^{deg θⱼ})) = (1 + x⁶ + x⁷ + x⁸ + x⁹ + x¹⁵) /
((1−x²)(1−x³)(1−x⁴)(1−x⁵)(1−x⁶))

y se comparó, término a término, contra la **serie de Molien calculada independientemente**
(sección `03-anillo-invariantes-s5.md`, vía autovalores exactos de las matrices, sin ninguna
relación con este cálculo de Singular). **Coincide exacto en los 12 grados probados** (0 a 11):
1, 0, 1, 1, 2, 2, 5, 4, 8, 9, 13, 15. Esta es la confirmación más fuerte posible: dos cálculos
completamente independientes (uno de teoría de representaciones/caracteres hecho a mano y en
Python, otro de un algoritmo de Gröbner/normalización de Noether en Singular) dan exactamente el
mismo resultado.

## Cómo reproducir

```bash
wsl -d Ubuntu-24.04 -- bash -c "cd '/mnt/c/Users/marti/Documents/repositorios/cuerdas/fase2-objetivo-tecnico/scripts' && Singular -q s5_invariant_ring_summary.sing"
```

para el resumen (grados y tamaños), o `s5_invariant_ring.sing` para los polinomios completos
(salida enorme, guardada en `scripts/s5_invariant_ring_output.txt` — el invariante secundario de
grado 15 por sí solo tiene 3565 términos, no apto para inspección manual).

## Qué significa esto para el objetivo final (clasificar S-matrices de 5 partículas)

Ahora tenemos la base algebraica **completa y rigurosa** sobre la que construir el "módulo local"
de S-matrices para n=5 (el análogo de la sección 2.4 de arXiv:1910.14392, que para n=4 usaba el
anillo libre de S₃ sobre (s,t,u)). La diferencia real con n=4:
- n=4: módulo libre sobre un anillo de polinomios de 2 generadores (grados 2 y 3).
- n=5: módulo sobre un anillo con 5 generadores primarios (grados 2,3,4,5,6) MÁS una estructura de
  6 "capas" secundarias (grados 0,6,7,8,9,15) — significativamente más rico, con el generador
  secundario de grado 15 siendo un objeto genuinamente complejo (miles de términos).

**Evaluación honesta**: esto confirma con certeza (no solo sospecha) que la extensión de n=4 a n=5
para clasificar S-matrices vía este método es un problema de escala mucho mayor de lo que sugería
la analogía inicial con S₃. Es factible en principio (Singular ya nos dio la base algebraica
completa), pero construir el módulo de S-matrices sobre esta base y después imponer la condición de
Regge growth va a requerir trabajar con polinomios de miles de términos — la siguiente etapa debería
hacerse también en Singular (o similar), no volviendo a sympy artesanal.

## Próximo paso concreto (no hecho, para decidir con el usuario)

1. Repetir para n=5 la construcción de "S-matrices locales de 4 fotones/gravitones" (sección 2.4-2.9
   del paper original) usando esta base algebraica ya calculada — esto requiere primero traducir
   la parte de "datos de polarización" (los ε_i, α_i del paper) a la representación de 5 partículas,
   que es un paso conceptual nuevo no cubierto todavía en este proyecto.
2. Evaluar, dado el tamaño de los invariantes encontrados (miles de términos), si conviene primero
   buscar un resultado más modesto/acotado (ej. contar dimensiones sin construir explícitamente los
   S-matrices completos) antes de intentar la clasificación completa.

## Archivos de esta sesión

- `scripts/s5_invariant_ring.sing` — script completo (genera los polinomios enteros).
- `scripts/s5_invariant_ring_summary.sing` — script de resumen (grados y tamaños, sin los
  polinomios completos).
- `scripts/s5_invariant_ring_output.txt` — salida completa guardada (incluye los polinomios).
