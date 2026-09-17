# Fase 2 — Paso 8 (inicio): estructura de datos de polarización para fotones, n=5

## Contexto y decisión

Tras confirmar que el caso escalar de n=5 ya está publicado (ver `10-hallazgo-ya-publicado-escalares.md`),
se decidió (junto con el usuario, sopesando honestamente la dificultad) seguir con el caso **con
spin**, empezando por **fotones** (spin 1) en vez de ir directo a gravitones (spin 2) — es el
escalón natural: más simple que gravitones, pero ya requiere generalizar el aparato de datos de
polarización que en el paper original (arXiv:1910.14392) solo se hizo para n=4.

## Cómo funciona para n=4 (releído a fondo, sección 2.1 del paper)

1. Los 4 momentos pᵢ (masa nula, con conservación de momento) abarcan genéricamente un subespacio
   de **3 dimensiones** de Minkowski D-dimensional (el "plano de scattering") — porque hay
   n−1=3 momentos independientes tras imponer conservación.
2. Cada vector de polarización εᵢ (D componentes, sujeto al gauge de Lorentz εᵢ·pᵢ=0) se separa en:
   - **εᵢ⊥**: la parte transversal al plano de scattering, viviendo en el espacio de
     **D−3 dimensiones** ortogonal a ese plano. Aquí actúa el "little group" SO(D−3) — la simetría
     de rotaciones transversales, que no afecta a los momentos.
   - **εᵢ∥**: la parte dentro del plano de 3 dimensiones.
3. La restricción εᵢ·pᵢ=0 fuerza a εᵢ∥ a vivir en un subespacio de **2 dimensiones** del plano de 3D.
4. La invariancia de gauge residual (shift εᵢ → εᵢ + ζᵢpᵢ) permite eliminar 1 dimensión más:
   el resultado es que εᵢ∥ queda parametrizado, salvo por la redundancia de gauge, por **un único
   número complejo αᵢ** por partícula (ecuación 2.11-2.12 del paper).
5. Entonces, los datos físicos gauge-invariantes por partícula son el par (εᵢ⊥, αᵢ): un vector en
   SO(D−3) más un escalar complejo.

## El hallazgo estructural para n=5 (nuevo, no está en ningún paper visto hasta ahora)

Repitiendo el mismo conteo con 5 partículas:

1. Los 5 momentos abarcan genéricamente un plano de scattering de **4 dimensiones** (n−1=4
   momentos independientes, para D≥4 — coincide con el régimen "sin condiciones de Gram" que ya
   usamos para los invariantes de Mandelstam).
2. εᵢ⊥ (transversal) vive ahora en un espacio de **D−4 dimensiones** — el little group relevante
   para n=5 es **SO(D−4)**, no SO(D−3).
3. εᵢ∥ (dentro del plano de 4D), tras imponer εᵢ·pᵢ=0 (1 restricción), vive en un subespacio de
   **3 dimensiones** del plano de 4D.
4. La invariancia de gauge residual (mismo shift εᵢ→εᵢ+ζᵢpᵢ) elimina 1 dimensión más.
5. **Resultado**: a εᵢ∥ le quedan **2 dimensiones** físicas gauge-invariantes por partícula — NO
   un solo número complejo αᵢ como en n=4, sino **un par de números** (llamémoslos provisionalmente
   (αᵢ, βᵢ), sin nombre definitivo todavía) por partícula.

**Esto es una diferencia estructural real, no solo de tamaño**: para n=4, el paso de "momentos" a
"datos gauge-invariantes en el plano" colapsa 2 dimensiones a 1 (un escalar complejo). Para n=5,
colapsa 3 dimensiones a 2 — el resultado NO es un escalar, es un objeto con 2 componentes, cuya
estructura de transformación bajo permutaciones de partículas y bajo el grupo de invariantes de
Mandelstam todavía hay que trabajar (análogo a la ecuación 2.16 del paper original, que da αᵢ
explícitamente en términos de contracciones del field strength Fᵢ con los otros momentos).

## Estado: esto es un punto de partida, no una solución

Este documento deja registrado el conteo dimensional (verificado a mano, siguiendo exactamente el
mismo método del paper original aplicado a n=5) pero **no** deriva todavía:
- La fórmula explícita análoga a (2.16)/(2.17) para estos 2 parámetros por partícula en términos de
  contracciones de momentos/field-strength.
- Cómo actúa S₅ sobre este nuevo objeto de 2 componentes (el análogo de la sección 2.2 del paper,
  que para n=4 mostró que la acción de permutaciones sobre αᵢ tiene una sutileza de signo —
  representación proyectiva de S₄).
- Cómo se combina esto con el anillo de invariantes de Mandelstam ya calculado (parte 7,
  `08-anillo-completo-via-singular.md`) para formar el módulo completo de S-matrices de fotones.

## Próximo paso concreto (para la siguiente sesión)

1. Derivar explícitamente la parametrización de los 2 parámetros gauge-invariantes por partícula en
   n=5 (el análogo de eq. 2.11 del paper), usando combinaciones de momentos apropiadas para el plano
   de 4 dimensiones.
2. Encontrar la fórmula tipo (2.16) que exprese estos parámetros en términos de contracciones
   Lorentz-invariantes de los momentos y el field strength — necesaria para poder escribir
   S-matrices manifiestamente invariantes.
3. Estudiar la acción de S₅ sobre estos nuevos datos (verificar si hay una sutileza de signo
   análoga a la de n=4, o algo distinto).

## Fuentes

- Chowdhury, Gadde, Gopalka, Halder, Janagal, Minwalla, arXiv:1910.14392, sección 2.1 (releída a
  fondo esta sesión, páginas 17-22 del PDF).
