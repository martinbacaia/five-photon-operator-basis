# Fase 2 — Paso 15: enumeración sistemática de candidatos de grado mínimo — corrección de conteo
# de grado, obstrucción general de antisimetría probada, y confirmación de grado 9 como piso real

## Corrección: error de conteo de grado en `18-generalizacion-FdotF-y-primer-candidato.md`

Ese archivo afirmaba `C^k_{ab}` tiene grado 2 en momentos. **Es falso**: contando potencias de
momento explícitamente en `C^k_{ab} = p_a.F^k.p_b = (p_k.p_a)(εₖ.p_b) − (p_k.p_b)(εₖ.p_a)`, cada
término tiene **3** momentos (p_k, p_a, p_b) — grado 3, no 2. `F^i:F^j` sí es grado 2 (dos
momentos). Esto corrige el grado del candidato ya encontrado: `S = F:F · F:F · C · (d_ik−d_jk)`
tiene grado `2+2+3+2 = 9`, no 7. Se deja constancia por honestidad — es un error propio, detectado
al re-hacer la cuenta con más cuidado antes de generalizar.

## El piso teórico de la topología "2 parejas F:F + 1 singlete C" es grado 7 — y se prueba que
## SIEMPRE se anula, para cualquier elección de contracción del singlete

Con el conteo correcto, la topología mínima que cubre las 5 polarizaciones exactamente una vez
(usando los 2 únicos "átomos" disponibles — `F^i:F^j`, grado 2, cubre 2 patas; `C^k_{ab}`, grado 3,
cubre 1 pata) tiene grado piso `2+2+3=7`.

**Se probaron 2 formas distintas de grado 7** (`scripts/enumeracion_candidatos_grado_minimo.py`):
- shape (a): singlete contraído con AMBOS momentos de la MISMA pareja, `C^k_{ij}` (la ya probada en
  la parte 18) — se anula idénticamente (∼10⁻¹¹, ruido de punto flotante).
- shape (c): singlete contraído con UN momento de CADA pareja, `C^k_{il}` — **también se anula
  idénticamente** (∼10⁻¹¹), resultado nuevo de esta sesión, no anticipado.

**Diagnóstico general (no solo observado — demostrado)**: `C^k_{ab}` es intrínsecamente
antisimétrico bajo el intercambio de sus 2 argumentos (`C^k_{ab}=−C^k_{ba}`, verificado
numéricamente: `C(1;2,3)=14.05, C(1;3,2)=−14.05`), mientras que `F^i:F^j` es simétrico en (i,j).
Para **cualquier** elección de {a,b} ⊂ {i,j,l,m} (los 4 índices que no son el singlete k), existe
una permutación de S₅ que (1) preserva la partición en 2 parejas + 1 singlete (dejando
`F^i:F^j·F^l:F^m` invariante) y (2) intercambia exactamente a↔b (invirtiendo el signo de `C^k_{ab}`):
- si {a,b}={i,j}: la transposición (i j) sola alcanza.
- si {a,b}={l,m}: la transposición (l m) sola alcanza.
- si {a,b} tiene un elemento de cada pareja (p.ej. {i,l}): la doble transposición (i l)(j m)
  intercambia a↔b y además intercambia las 2 parejas entre sí — pero como el producto
  `F^i:F^j·F^l:F^m` es simétrico bajo ese intercambio de parejas (conmutatividad del producto), el
  factor de parejas queda invariante mientras `C^k_{il}→C^k_{li}=−C^k_{il}` cambia de signo. Se
  verificó explícitamente esta identidad (`T2(τ') = −T2(τ)` para τ'=τ∘(il)(jm)).

Como esa permutación pertenece a la órbita completa que se suma, cada término se cancela
exactamente con su pareja bajo esa permutación → **la suma sobre toda la órbita de S₅ es
idénticamente cero para CUALQUIER elección de contracción del singlete, en la topología de grado
7**. Esto no es un intento fallido más — es un resultado estructural general (probado, no solo
observado en 2 casos) que descarta por completo el grado 7 para esta topología.

## Por qué la corrección mínima da exactamente grado 9 (y no 8)

Para "arreglar" la cancelación hace falta multiplicar por un factor que sea también antisimétrico
bajo la MISMA permutación que mata el término (p.ej. i↔j para shape (a)). Un escalar de Lorentz
construido solo con productos punto de momentos siempre tiene grado **par** (cada contracción
empareja 2 vectores) — no existe un invariante escalar de grado impar (en particular grado 1) que
pueda servir de compensador, salvo que se recurra al tensor de Levi-Civita (dependiente de D, fuera
del espíritu "D genérico" del paper original). El compensador de grado mínimo posible es entonces
grado 2: `(d_ik − d_jk)` (antisimétrico bajo i↔j, construido con el anillo de Mandelstam ya
calculado). Esto da grado total `7+2=9` — el mismo candidato ya encontrado y verificado no-trivial
en la parte 18. **Conclusión: grado 9 es el piso real para esta topología**, no una elección
arbitraria.

## Los 2 candidatos de grado 9 encontrados son EL MISMO invariante (proporcionales) — no una base
## de dimensión 2

Se probó una segunda variante de grado 9, usando la shape (c) corregida:
`T_grado9_v2 = F^i:F^j·F^l:F^m·C^k_{il}·(d_lk−d_mk)`, sumada sobre la órbita completa de S₅, y
comparada contra la ya verificada `T_grado9_v1 = F^i:F^j·F^l:F^m·C^k_{ij}·(d_ik−d_jk)`.

**Resultado**: el cociente `S1/S2` da **exactamente 2.000000** en 7 semillas distintas (cinemática
numérica independiente cada vez) — son el **mismo invariante** (hasta normalización), no 2
elementos independientes de una base. Es decir, la construcción por sumas de órbita, aplicada a
"semillas" que difieren solo en la elección interna de contracción, colapsa al mismo resultado final
(fenómeno esperado del operador de Reynolds: puede proyectar semillas distintas al mismo elemento
del espacio de invariantes si están relacionadas por la acción del grupo).

**Conclusión actual**: se tiene **1 solo** invariante multilineal S₅-invariante confirmado (no una
base) de grado 9, hasta ahora.

## Chequeo de una topología alternativa de grado naive más bajo — también se anula, y su corrección
## costaría más caro que grado 9

Se exploró si una topología distinta podría dar grado piso menor: en vez de "F:F pareja + C
singlete", usar una traza de 3 field-strengths `Tr(F^i F^j F^k) := F^i_{μν}F^j^{νρ}F^k_ρ^μ` (un
objeto que cubre 3 polarizaciones a la vez, con grado piso 3 — no necesita momentos externos, solo
usa los momentos propios ya contenidos en cada F), combinada con 1 pareja `F^l:F^m` (grado 2) para
las 2 patas restantes: grado naive total `3+2=5`, mucho menor que 7.

**Verificado** (numéricamente, con matrices D=6 explícitas): esta traza es efectivamente
antisimétrica bajo el intercambio de cualquier par de sus 3 índices (`Tr(F¹F²F³)=−Tr(F¹F³F²)`,
diferencia ∼10⁻¹⁵ — propiedad estándar de trazas de 3 matrices antisimétricas, ya que
`tr(ACB)=−tr(ABC)` para matrices antisimétricas reales). Por el mismo mecanismo que arriba, la suma
sobre la órbita completa de S₅ de `Tr(F^i F^j F^k)·F^l:F^m` se anula idénticamente — **verificado
numéricamente** (∼10⁻¹¹, ruido de punto flotante, en 3 semillas).

**Costo de la corrección**: para compensar la antisimetría total bajo las 6 permutaciones de
(i,j,k) (no solo una transposición) hace falta un factor construido con Mandelstam que sea también
totalmente antisimétrico en 3 índices — el candidato de grado mínimo para eso es un determinante
3×3 de invariantes tipo `d_ab` (análogo a un discriminante de Vandermonde), que tiene grado
`3×2=6` (determinante de una matriz 3×3 de entradas de grado 2). Esto daría un total de
`3+2+6=11` — **más caro que 9**. Por lo tanto esta topología alternativa no mejora el resultado ya
encontrado (no se llevó la verificación numérica de esta corrección de grado 11 hasta el final, dado
que ya es peor que la opción de grado 9 ya confirmada — no vale la pena el esfuerzo de verificarla
en detalle).

## Estado honesto

- **Grado 9 confirmado como el mínimo real** (no arbitrario) alcanzable con las construcciones
  exploradas hasta ahora, con exactamente 1 invariante no-trivial encontrado (no una base de
  dimensión mayor).
- **No se ha probado que grado 9 sea el mínimo ABSOLUTO** sobre todas las topologías posibles
  (solo se descartaron 2: "2 parejas + C singlete" sin corrección, que es imposible en cualquier
  grado 7; y la traza de 3 F's, que costaría grado ≥11) — quedan otras combinaciones de átomos por
  explorar (p.ej. el objeto "sandwich" `p_a.F^i.F^j.p_b`, grado 4, bilineal en εᵢ,ε_j, con una
  estructura de contracción distinta a F^i:F^j, que podría combinarse de otras formas).
- **No se sabe todavía si el espacio de invariantes de grado 9 tiene dimensión 1 o más** — falta
  buscar semillas genuinamente distintas (no solo variar la elección de contracción del singlete,
  que ya demostró dar el mismo resultado).

## Qué falta (próximo paso concreto)

1. Buscar más "semillas" de grado 9 genuinamente distintas (variando la elección de qué pareja
   lleva el factor Mandelstam corrector, o usando combinaciones distintas como
   `(d_il-d_jl)` en vez de `(d_ik-d_jk)`) para determinar si el espacio de invariantes de grado 9
   tiene dimensión 1 o más.
2. Considerar seriamente calcular la serie de Hilbert/Molien para el caso CON SPIN (little group
   SO(D−4) para n=5, ver parte 11) — sin esto, seguir "adivinando" candidatos uno por uno es
   ineficiente; con ella sabríamos de antemano cuántos invariantes independientes existen en cada
   grado, como ya se hizo para el caso escalar (partes 3-8). Es un paso más difícil (la teoría de
   invariantes con spin no tiene la maquinaria de Molien tan directa como el caso escalar) pero
   sería el cambio de estrategia más rentable en este punto.
3. Verificar Regge growth sobre el candidato de grado 9 ya encontrado, aunque no se sepa aún si es
   parte de una base completa — al menos como primer chequeo del objetivo final.

## Archivos

- `scripts/enumeracion_candidatos_grado_minimo.py` — candidatos de grado 7 (ambas shapes, ambas se
  anulan) y grado 9 (v1, v2 — proporcionales, mismo invariante).
- Verificación de la traza de 3 field-strengths hecha inline (no se generó script separado, cálculo
  corto reutilizando `gen_kinematics` de `enumeracion_candidatos_grado_minimo.py`).
