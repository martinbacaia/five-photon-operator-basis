# Fase 2 — Paso 2: la representación de S₅ sobre los invariantes de Mandelstam de 5 partículas

## Resultado nuevo de esta sesión

**Para scattering de 5 partículas sin masa, el espacio de invariantes de Mandelstam independientes
(dimensión 5, por la fórmula general n(n-3)/2) es una representación IRREDUCIBLE de S₅ de dimensión
5** — no una simple permutación de 5 objetos entre sí (a diferencia de n=4, donde el espacio de
2 invariantes independientes (s,t) sí resultó ser la restricción de una representación irreducible
de S₄ de dimensión 2, que a su vez restringe a la "2M" de S₃ ya usada en el paper de referencia).

Esto es exactamente la pieza que el paper original (arXiv:1910.14392) no necesitaba derivar para
n=4 (porque ahí S₄/(ℤ₂×ℤ₂)=S₃ actúa de forma simple, permutando literalmente (s,t,u) entre sí) pero
que SÍ hace falta derivar para n=5, porque no existe una base natural de 5 invariantes que S₅
simplemente permute entre sí — la acción es genuinamente más rica (es la representación de Specht
correspondiente a la partición [3,2] de S₅, aunque para el uso práctico no hace falta ese nombre,
alcanza con tener las matrices explícitas, que sí se calcularon).

## Cómo se derivó (método, no solo resultado)

1. Los C(n,2) invariantes "crudos" s_ij = (p_i+p_j)² forman la representación de permutación de S_n
   sobre los 2-subconjuntos de {1,...,n} (las "aristas" del grafo completo K_n).
2. La conservación de momento (Σp_i=0, p_i²=0) implica, para cada partícula k, que
   Σ_{j≠k} s_kj = 0 — esto es exactamente el **mapa de incidencia** arista→vértice de K_n, que es
   equivariante bajo S_n (permutar los vértices conmuta con tomar la suma de aristas incidentes).
3. El espacio de invariantes independientes es el **núcleo** de ese mapa. Por ser el núcleo de un
   mapa equivariante, es automáticamente un sub-espacio invariante bajo S_n — es decir, es una
   representación genuina, no solo un espacio vectorial cualquiera.
4. Se calculó explícitamente (código, no a mano) la base del núcleo (vía eliminación gaussiana
   exacta con `sympy`) y la matriz de acción de cada permutación sobre esa base, para n=4 y n=5.

Implementado en `scripts/n5_mandelstam_representation.py`.

## Verificación (3 chequeos independientes, ninguno de ellos "confiar en que corrió")

1. **Fórmula de carácter derivada a mano, verificada contra las matrices explícitas**: se dedujo
   analíticamente que χ(g) = C(fix(g),2) + c₂(g) − fix(g), donde fix(g) = número de puntos fijos de
   g y c₂(g) = número de 2-ciclos en g (viene de restar el carácter de la representación de vértices
   del carácter de la representación de aristas). Se comparó esta fórmula, para las 5 clases de
   conjugación de S₄ y las 7 de S₅, contra la traza real de las matrices calculadas por código —
   **coincide exacto en las 12 clases probadas**.
2. **Irreducibilidad confirmada por norma de carácter**: ‖χ‖² = (1/n!)Σ|clase|·χ(g)² dio
   **exactamente 1** tanto para n=4 como para n=5 — condición necesaria y suficiente de
   irreducibilidad. No es una representación reducible con piezas más chicas escondidas.
3. **Test de homomorfismo por fuerza bruta**: se generaron pares aleatorios de permutaciones de S₅
   (40 pares) y se verificó que ρ(g₁)ρ(g₂) = ρ(g₁g₂) para las matrices calculadas — **0 fallos**.
   Esto confirma que lo que se calculó es una representación de grupo genuina y no una asignación
   de matrices que por casualidad tiene las trazas correctas.

## Un bug real encontrado y corregido durante la verificación cruzada (importante, dejar constancia)

Al intentar contrastar el resultado de n=4 (recalculado hoy vía el núcleo de incidencia) contra la
representación "2M" ya validada la sesión anterior, se armó una segunda derivación "física" de la
acción de S₄ sobre (s,t) = (s₁₂,s₁₃), usando la conservación de momento explícita
(s₁₂=s₃₄, s₁₃=s₂₄, s₁₄=s₂₃=-s₁₂-s₁₃, derivada a mano y verificada esta sesión). Esta segunda
derivación, en su primera versión, usaba la fórmula `nuevo_s_ij = viejo_s_{σ⁻¹(i),σ⁻¹(j)}`.

Al buscar una matriz de cambio de base P que hiciera **ambas** representaciones (la del núcleo, y
la "física") coincidir para los 6 elementos de S₃ simultáneamente, P funcionaba para (12) y (123)
pero **fallaba para (13) y (23)** — un resultado real, no un error de tipeo, que había que
diagnosticar en vez de ignorar.

**Diagnóstico** (verificado por computo, no supuesto): se probó explícitamente si la asignación
"física" cumplía ρ(g₁)ρ(g₂)=ρ(g₁g₂) (0 de 18 pares cumplía) o si en cambio cumplía la regla de
**anti-homomorfismo** ρ(g₁)ρ(g₂)=ρ(g₂g₁) (18 de 18 pares cumplía). Confirmado: la fórmula usaba la
convención de acción por la DERECHA (con σ⁻¹) en vez de la izquierda, mientras que la construcción
del núcleo de incidencia (`action_on_kernel`) usa la convención estándar de acción por la izquierda
P(σ)e_i = e_{σ(i)}. **Corrección**: cambiar `σ⁻¹(i),σ⁻¹(j)` por `σ(i),σ(j)` en la fórmula física.
Tras la corrección, el test de homomorfismo pasó 100%, y la búsqueda de P encontró una matriz válida
que conjuga los 6 elementos simultáneamente (P=[[1,2],[2,1]], det=-3≠0) — **confirmando que las 2
representaciones (núcleo de incidencia, y acción física directa sobre s₁₂,s₁₃) son genuinamente
isomorfas**, no solo coincidentes en traza/determinante (que por sí solo no prueba nada, porque
traza y determinante son constantes en cada clase de conjugación y no distinguen relabelings).

**Lección metodológica** (aplicable a partir de ahora): cuando se comparan 2 construcciones de una
misma representación, coincidencia de traza/determinante elemento por elemento NO es suficiente —
hace falta o bien un test de homomorfismo explícito, o encontrar una matriz de cambio de base única
que funcione simultáneamente para un conjunto generador del grupo.

Ver `scripts/crossverificacion_n4_s3.py` (versión final, corregida) para el código completo.

## Qué significa esto para el próximo paso

Ya tenemos, para n=5, las matrices explícitas de la acción de S₅ sobre el espacio de invariantes de
Mandelstam independientes (una representación irreducible de dimensión 5, verificada). El siguiente
paso (todavía no hecho) es la parte que de verdad clasifica S-matrices: repetir para n=5 la
construcción de la sección 2.4 del paper (arXiv:1910.14392) — el "módulo local" de S-matrices como
módulo sobre el anillo de polinomios de los invariantes de Mandelstam. La diferencia clave con n=4:
para n=3 invariantes (s,t,u permutados simplemente por S₃), el anillo de invariantes de S₃ es un
anillo de polinomios libre (generado por 2 invariantes elementales, ver `s3_partition_functions.py`,
consistente con que S₃ actuando así es un **grupo de reflexiones** en el sentido de
Chevalley-Shephard-Todd). **Para la representación de 5 dimensiones de S₅ encontrada hoy, hace falta
determinar primero si el anillo de invariantes de S₅ sobre esta representación específica es
también un anillo de polinomios libre (grupo de reflexiones en esa representación) o si tiene
relaciones no triviales** — esto no puede asumirse por analogía, hay que verificarlo. Es el próximo
paso concreto de trabajo.

## Fuentes / archivos de esta sesión

- `scripts/n5_mandelstam_representation.py` — deriva y verifica la representación para n=4 y n=5.
- `scripts/crossverificacion_n4_s3.py` — cross-validación n=4 vs. la "2M" de la sesión anterior
  (incluye el bug de convención encontrado y corregido, documentado en el propio código).
