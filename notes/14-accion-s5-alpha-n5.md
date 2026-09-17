# Fase 2 — Paso 10: acción de S₅ sobre (α₁, α₁') — hallazgo revisado (v2): el subgrupo S₃ que fija
# 1 y un "leftover" auxiliar actúa por la representación irreducible 2M genuina, sin proyectividad

## Resumen del giro en esta sección

Un primer intento (v1, ver más abajo "Historial") pareció mostrar que ninguna elección de pivote da
una acción covariante de S₅ sobre (α₁,α₁'). **Ese primer intento tenía un error de método real**
(la misma clase de error ya documentado en `02-representacion-s5-mandelstam.md` de esta sesión: usar
la matriz M "vieja" en vez de la matriz M evaluada en la configuración permutada). Al corregirlo,
aparece un resultado limpio y positivo: **el subgrupo S₃ ⊂ S₄ que permuta {2,3,4} (fijando 1 y 5)
actúa sobre (α₁,α₁') exactamente como la representación irreducible 2-dimensional "2M" de S₃** —
la misma representación que actúa sobre (s,t) en el caso n=4 del paper original — **sin ningún
signo ±1 de por medio** (a diferencia de la anomalía de fase encontrada en la sección 2.2 del paper
para n=4).

## Construcción corregida

(α₁,α₁') se definen resolviendo `(C_{2,3}, C_{2,4}) = M·(α₁,α₁')`, con M la matriz 2×2 hallada en
`13-formula-inversa-polarizacion-n5.md`, construida con pivote=partícula 2, compañeras=3,4, y
partícula 5 como "leftover" (no usada en la construcción). Al aplicar una permutación σ que fija 1
y 5 (permuta solo {2,3,4}), el pivote/compañeras de la config *relabeled* son (σ(2),σ(3),σ(4)) —
**el error del primer intento fue usar la matriz M vieja (con el pivote 2 original) en vez de la
matriz M_σ evaluada en el NUEVO pivote/compañeras** (M_σ = M_para_pivote(σ(2),σ(3),σ(4))). Con la
matriz correcta:

```
rho(sigma) := M_sigma^{-1} · N(sigma),   N(sigma) = coeficientes de (C_{σ(2)σ(3)}, C_{σ(2)σ(4)})
                                                      en la base (alpha1,alpha1') vieja
```

**Resultado: para las 6 permutaciones de {2,3,4} (fijando 5), ρ(σ) es una matriz CONSTANTE**
(verificado simbólicamente Y numéricamente en 2 puntos cinemáticos distintos para descartar
coincidencia con simplificación incompleta).

## Corrección de convención: homomorfismo vs anti-homomorfismo

Al verificar que ρ define una representación genuina (chequeo de cierre de grupo, no solo
individual), se encontró que **ρ, tal como se definió, es un ANTI-homomorfismo**
(`ρ(σ)ρ(τ) = ρ(τ∘σ)` para los 36 pares probados, en vez de `ρ(σ)ρ(τ)=ρ(σ∘τ)`) — exactamente la
misma clase de error de convención (σ vs σ⁻¹) ya encontrada y corregida en
`02-representacion-s5-mandelstam.md` de esta sesión. **Definiendo `ψ(σ) := ρ(σ⁻¹)` se obtiene un
homomorfismo genuino**, verificado explícitamente en los 36 pares (σ,τ) de S₃: `ψ(σ)ψ(τ)=ψ(σ∘τ)`
exacto (no aproximado, sin usar simplify — comparación de matrices con entradas racionales
exactas).

## Identificación de la representación

Las 6 matrices ψ(σ) (para las 6 permutaciones de {2,3,4}):

| σ (en {2,3,4}, fija 1,5) | tipo    | ψ(σ)                | traza | det |
|---|---|---|---|---|
| identidad                | —       | [[1,0],[0,1]]       | 2     | 1   |
| (2 3)                    | transp. | [[-1,-1],[0,1]]     | 0     | -1  |
| (2 4)                    | transp. | [[0,1],[-1,-1]]     | 0     | -1  |
| (3 4)                    | transp. | [[0,1],[1,0]]       | 0     | -1  |
| (2 3 4)                  | 3-ciclo | [[-1,-1],[1,0]]     | -1    | 1   |
| (2 4 3)                  | 3-ciclo | [[1,0],[-1,-1]]     | -1    | 1   |

Estos caracteres (2, 0, 0, 0, -1, -1) coinciden **exactamente** con la tabla de caracteres de la
representación irreducible 2-dimensional **2M** de S₃ (ec. 2.33 del paper: χ(e)=2, χ(transposición)=0,
χ(3-ciclo)=-1) — la misma representación que actúa sobre (s,t) para n=4 (sección 2.7 del paper). Se
verificó irreducibilidad directamente: los autovectores de dos transposiciones distintas (p.ej.
(3 4) → autovectores (1,1),(1,-1); (2 3) → autovectores (1,-2),(1,0)) no coinciden, luego no hay
subespacio invariante común — la representación es irreducible, no es suma de 2 representaciones
de dimensión 1.

**Conclusión positiva de esta sección**: para el subgrupo que fija tanto la partícula "activa" (1)
como una partícula auxiliar "leftover" (5, en esta convención), el par (α₁,α₁') transforma
**exactamente igual que (s,t) transforma bajo S₃ en el caso n=4** — sin ninguna anomalía de signo.
Esto es una fuerte confirmación de que la construcción (pivote + compañeras + 1 leftover) es la
generalización correcta, al menos para esta parte del grupo.

## Continuación: función de transición explícita entre distintos "leftover" — resuelto el punto 1

Se calculó explícitamente la función de transición entre la trivialización con leftover=5
(pivote=2, compañeras=3,4 — la usada en todo lo anterior) y la trivialización con leftover=4
(pivote=2, compañeras=3,5), **sin permutar nada todavía** — es simplemente el mismo vector físico
ε₁^∥ expresado en dos bases distintas del mismo plano 2D.

**Método**: se expresó w₂₅ (el vector que reemplaza a w₂₄ cuando el leftover pasa de 5 a 4) como
combinación lineal de la base {w₂₃, w₂₄, p₁} — encontrando (ver script)

```
w₂₅ = λ₁ w₂₃ + λ₂ w₂₄ + μ p₁,     λ₁ = −s₁₃/s₁₅,   λ₂ = −s₁₄/s₁₅,   μ = 1/s₁₅
```

(usando p₁ = −(p₂+p₃+p₄+p₅) por conservación de momento, y s₁ⱼ := −2 p₁.pⱼ). Igualando las dos
expresiones de ε₁^∥ (una en cada base) y despejando:

```
alpha1^(4)  = alpha1^(5) − (s13/s14) alpha1'^(5)
alpha1'^(4) = −(s15/s14) alpha1'^(5)
```

es decir, en forma matricial, T = [[1, −s₁₃/s₁₄],[0, −s₁₅/s₁₄]] (triangular superior,
det(T) = −s₁₅/s₁₄).

**Verificación**: se recalculó (α₁^(4), α₁'^(4)) de forma completamente independiente —
resolviendo el sistema 2×2 con la matriz M para el pivote/compañeras (2,3,5) y los datos
C₂₃, C₂₅ escritos en la base (α₁^(5),α₁'^(5)) — y coincide **exactamente** (diferencia simbólica
cero, sin redondeo) con la fórmula de arriba.

**Por qué esto es una buena noticia**: la dependencia cinemática encontrada en el primer intento
(v1) no es un desastre incontrolable — es una función de transición GL(2) **explícita y simple**
(solo 2 cocientes de invariantes s₁ⱼ, sin raíces cuadradas ni fracciones complicadas), exactamente
del mismo tipo de objeto que aparece en la construcción de fibrados vectoriales con trivializaciones
locales. Esto resuelve el punto 1 de la lista de pendientes de la sección anterior.

## Lo que falta (honesto, no resuelto)

Falta entender la acción de las permutaciones que **mueven la partícula "leftover" (5)** — es
decir, el resto de S₄ (que tiene orden 24, y el subgrupo hallado tiene orden 6: falta un factor de
índice 4). Estas permutaciones cambian CUÁL partícula juega el rol de "leftover", y como tal
requieren comparar dos construcciones con leftovers distintos (p.ej. leftover=5 vs. leftover=4) —
la comparación directa (fijando de antemano una sola matriz M) sí da dependencia cinemática (esto
es real, no un bug: son bases genuinamente distintas del mismo plano 2D). La estructura sugerida es
la de un **fibrado**: para cada elección de "leftover" j ∈{2,3,4,5} hay una trivialización local
(α_1^{(j)}, α_1'^{(j)}) del plano físico, con funciones de transición kinematics-dependientes entre
trivializaciones distintas, y S₄ actúa permutando las trivializaciones Y mezclando las coordenadas
según esas funciones de transición.

1. ~~Calcular explícitamente las funciones de transición entre trivializaciones con leftover
   distinto~~ — **hecho, ver sección anterior**: son simples (cocientes de s₁ⱼ, sin raíces).
2. Usando la función de transición ya calculada, construir la matriz completa ρ(σ) para las
   permutaciones de S₄ que SÍ mueven la partícula 5 (p.ej. (4 5), (2 5), el 4-ciclo), combinando
   (i) la reetiqueta cruda de contracciones (como en v1) con (ii) la función de transición T entre
   leftovers — análogo a como se combinó M_σ en la sección "Continuación" de arriba, pero ahora
   con dos "saltos": uno de re-pivoteo (dentro del mismo leftover) y otro de cambio de leftover.
   Verificar si el resultado, así construido con cuidado, sí da matrices constantes (cerrando la
   representación completa de S₄), o si hace falta aceptar que la única descripción limpia es la
   acción combinada (permutación de leftover) + (transición GL(2)), es decir una representación de
   S₄ que actúa en un fibrado de rango 2 sobre las 4 elecciones de leftover en vez de en un espacio
   vectorial fijo de dimensión 2.
3. Verificar si el conjunto total (4 trivializaciones × 2 componentes = 8 números, con relaciones)
   admite una descripción limpia como una única representación de S₄ de dimensión mayor (candidato
   natural: 8 = 2×4, quizás relacionado con inducir la representación 2M de S₃ hasta S₄ — el
   carácter inducido se puede calcular con la fórmula de Frobenius y comparar).
4. Solo después de esto, extender a permutaciones que mueven la partícula 1 misma (el resto de
   S₅), análogo pleno de la sección 2.2 del paper.

## Historial (v1, error corregido)

La primera versión de este análisis (antes de la corrección) concluía —incorrectamente— que
NINGUNA elección de pivote fijo era covariante, ni siquiera para el subgrupo S₃ que ahora sabemos
que SÍ funciona. El error: al calcular ρ(σ), se resolvía el sistema con la matriz M **vieja** (con
pivote fijo en la partícula "2" original) en vez de reevaluar M en el pivote/compañeras de la
config permutada (M_σ). Al mezclar una matriz M de una convención con datos C_jk de otra, el
resultado espurio dependía de la cinemática — un artefacto del error, no un hecho físico. Lección
metodológica (ya anotada una vez en esta sesión, y ahora reconfirmada): **al construir la acción de
un grupo sobre datos derivados de una convención dependiente de índices (aquí: "pivote"), hay que
re-evaluar TODA la convención (incluida la matriz que define los parámetros) en la configuración
permutada, no solo los datos crudos.**

## Continuación 2: intento de cerrar el S₄ completo (v3/v4) — resultado "24/24 constante" es un
## FALSO POSITIVO tautológico, detectado y corregido antes de darlo por bueno

Se intentó extender la representación más allá del S₃ que fija 5, agregando una función de
transición extra T (de la trivialización que induce σ de vuelta a la trivialización de referencia)
sobre el cálculo de `rho_v2`. El resultado (script `s5_accion_alpha_n5_v4.py`) daba matrices
constantes para los **24/24** elementos de S₄ — a primera vista, un cierre completo de la
representación.

**Al investigar por qué, antes de reportarlo como logro, se encontró que es un artefacto
tautológico, no un resultado físico real.** Razón: `M_σ⁻¹·N(σ)` (la construcción de `rho_v2`, SIN
ninguna T adicional) ya es, por pura álgebra lineal, *idéntica* a la función de transición directa
`transition_triple(REF, (σ(2),σ(3),σ(4)))` — ambas calculan lo mismo (expresar el mismo vector
físico εₚ₁∥, fijo, en el chart destino) por dos caminos matemáticamente equivalentes. Esto se
verificó explícitamente: para σ=(2 3), `M_σ⁻¹·N(σ)` y `transition_triple((2,3,4),(3,2,4))` dan
**exactamente la misma matriz** `[[-1,-1],[0,1]]`. Por lo tanto, aplicarle a `rho_v2` una T
*adicional* de vuelta al chart de referencia (`T · M_σ⁻¹·N(σ)`) es multiplicar una transición por su
propia inversa — el resultado es **la identidad para cualquier σ, sin importar la física real**
(se comprobó explícitamente: para σ=(2 3) y σ=(3 4), ambas dieron `rho = I`, lo cual contradice
directamente el resultado ya verificado de la sección anterior, donde esas mismas transposiciones
dan matrices no triviales de traza 0 — la contradicción fue la señal de alarma que llevó a
encontrar el error).

**Conclusión correcta (reemplaza la v3/v4)**: NO existe, con esta construcción, una manera de
obtener una matriz constante de S₄ completo actuando en un único espacio vectorial fijo de 2
dimensiones. Lo que sí existe, y es genuino:

1. El subgrupo S₃ (fija el leftover=5) actúa por la representación 2M genuina — **esto sigue en
   pie, ya verificado como homomorfismo, sección anterior**.
2. Las permutaciones que mueven la partícula leftover cambian de "fibra" (de trivialización), y la
   comparación entre trivializaciones distintas está dada por funciones de transición GL(2)
   genuinamente kinematics-dependientes — **calculadas explícitamente para las 3 necesarias**:

```
T(5->4) = [[1, -d13/d14], [0, (d12+d13+d14)/d14]]     = [[1,-s13/s14],[0,-s15/s14]]  (ya verificada antes)
T(5->3) = [[-d14/d13, 1], [(-d12-d13-d14)/d13, 0]]     (nueva)
T(5->2) = [[d14/d12, (d12+d14)/d12], [(-d12-d13-d14)/d12, (-d12-d13-d14)/d12]]  (nueva)
```

(dᵢⱼ:=pᵢ.pⱼ, notación de siempre; recordar d12+d13+d14=-d15 por conservación de momento). Estas
completan el mapa de las 4 trivializaciones (leftover=2,3,4,5) del mismo plano físico 2D.

**Interpretación honesta y precisa**: la estructura NO es una representación lineal ordinaria de
S₄ sobre un espacio fijo de dimensión 2 (ni de ninguna dimensión fija construida ingenuamente
concatenando las 4 trivializaciones, porque las funciones de transición no son constantes). Es un
**fibrado vectorial de rango 2, S₄-equivariante, sobre el espacio de cocientes S₄/S₃ (4 puntos —
las 4 elecciones de leftover)**, con S₃ actuando en la fibra por la rep. 2M y las funciones de
transición dadas arriba conectando fibras distintas. Este es un objeto matemático perfectamente
respetable (más rico que una representación lineal simple), pero significa que "encontrar la
acción de S₅ sobre (α₁,α₁')" en el sentido ingenuo (una matriz fija por cada permutación) **no
tiene solución** con esta construcción — hay que o bien (a) trabajar directamente con el fibrado,
o (b) buscar una renormalización kinematics-dependiente de (α₁,α₁') — análoga al factor √(st/u)
que ya usa el paper en su ecuación (2.11) para n=4 — que trivialice las funciones de transición y sí
dé una representación lineal genuina. La pista más concreta para (b): las 3 funciones de transición
tienen determinantes muy simples (`(d12+d13+d14)/d1j` para j=2,3,4, es decir `-s15/s1j`) — sugiere
que reescalando α^{(j)} por alguna potencia de s1j (o de s15) los determinantes podrían volverse 1,
un primer paso hacia trivializar.

**Lección metodológica (tercera vez en la sesión que aparece esta clase de error)**: al construir
la acción de un grupo sobre datos definidos vía una convención con múltiples "grados de libertad"
(aquí: elección de pivote Y elección de leftover), **hay que verificar que la matriz resultante NO
sea una tautología antes de aceptar un resultado "demasiado bueno" (24/24 constante)** —
contrastando contra un caso ya conocido y verificado independientemente (aquí: las transposiciones
de S₃ ya verificadas, que deberían seguir dando lo mismo, y no dieron — esa fue la señal).

## Archivos

- `scripts/s5_accion_alpha_n5_v2.py` — construcción corregida (M_σ), matrices ψ(σ) para S₃⊂S₄,
  verificación de homomorfismo genuino en los 36 pares, identificación de caracteres. **Válido.**
- `scripts/s5_accion_alpha_n5_v3.py`, `s5_accion_alpha_n5_v4.py` — intento de cerrar el S₄
  completo agregando una transición T; contienen las 3 funciones de transición (5→2,5→3,5→4)
  **que sí son válidas y reutilizables**, pero su conclusión final ("24/24 constante") es un
  **falso positivo tautológico** — ver explicación arriba. No usar `rho_v3`/`rho_v4` como
  representación de S₄; sí usar `transition_matrix`/`transition_triple` para las funciones de
  transición.
- `scripts/s5_accion_alpha_n5.py` — versión v1 (con el error original), conservada por
  trazabilidad; NO usar sus conclusiones (ver "Historial" arriba).
