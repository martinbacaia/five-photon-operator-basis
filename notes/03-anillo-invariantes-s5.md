# Fase 2 — Paso 3: el anillo de invariantes de S₅ NO es libre (hallazgo que aumenta la dificultad real)

## Pregunta que se atacó

La que quedó pendiente al cerrar el paso anterior (`02-representacion-s5-mandelstam.md`): el
"módulo local" de S-matrices para n=4 partículas (arXiv:1910.14392, sección 2.4) se construye sobre
el anillo de polinomios de (s,t,u), que es un anillo **libre** (polinomios sin relaciones,
generado por 2 invariantes elementales de S₃: s²+t²+u² en grado 2, y stu en grado 3 — ver la serie
`D=1/((1-x⁴)(1-x⁶))` ya validada en `01-validacion-s3-modulo.md`). ¿Es el anillo de invariantes de
S₅ sobre el espacio de 5 invariantes de Mandelstam (dimensión 5, encontrado en el paso anterior)
también un anillo libre, o tiene relaciones?

## Método: criterio de Chevalley-Shephard-Todd

Teorema clásico de teoría de invariantes: el anillo de polinomios G-invariantes sobre un espacio
vectorial V es un anillo de polinomios libre **si y solo si** G, actuando en V, está generado por
**pseudo-reflexiones** (elementos que fijan un hiperplano — es decir, cuyo autovalor 1 tiene
multiplicidad exactamente dim(V)−1).

- Para S₃ actuando sobre (s,t) (la "2M", dim 2): las 3 transposiciones SÍ son reflexiones genuinas
  (fijan una recta, autovalores (1,−1)) — S₃ en esa representación es literalmente el grupo diédrico
  de 6 elementos actuando como grupo de reflexiones del plano. Por eso el anillo es libre.
- Para S₅ actuando sobre el espacio de 5 invariantes de Mandelstam (dim 5): había que verificar si
  ALGÚN elemento de S₅ actúa como reflexión en ESA representación específica (no en la
  representación estándar de permutación de 5 puntos, que es una cosa distinta).

## Resultado: NINGÚN elemento de S₅ es una reflexión en esta representación

Se calcularon los autovalores exactos (con `sympy`, sin aproximación numérica) de la matriz de
acción de un representante de cada una de las 7 clases de conjugación de S₅ sobre el espacio de
5 dimensiones. Para ser una reflexión, un elemento necesita autovalor 1 con multiplicidad 4 (y el
quinto autovalor, por ser una representación real de un elemento de orden finito, tendría que ser
necesariamente −1). Resultado:

| clase | autovalores | ¿reflexión? |
|---|---|---|
| identidad | 1(×5) | no aplica |
| transposición | −1(×2), 1(×3) | **no** (fija solo 3 dim, no 4) |
| doble transposición | −1(×2), 1(×3) | **no** |
| 3-ciclo | 1(×1), 2 pares de complejos conjugados | **no** |
| (3,2) | 1(×1), 4 raíces de la unidad complejas | **no** |
| 4-ciclo | −1(×2), 1(×1), ±i | **no** |
| 5-ciclo | 1(×1), 4 raíces complejas de orden 5 | **no** |

**Ningún elemento no trivial tiene autovalor 1 con multiplicidad 4.** El caso más cercano son las
transposiciones, que fijan un subespacio de dimensión 3 (no 4) — actúan más "profundamente" en esta
representación que en la de permutación estándar.

**Conclusión, por Chevalley-Shephard-Todd**: el anillo de polinomios S₅-invariantes sobre el espacio
de 5 invariantes de Mandelstam **no es un anillo de polinomios libre** — tiene relaciones no
triviales entre sus generadores. Esto es cualitativamente distinto de n=4, y significa que el
"módulo local" análogo para n=5 no se puede construir simplemente enumerando generadores libres
como en la sección 2.4 del paper — hace falta encontrar tanto los generadores COMO las relaciones
del anillo de invariantes (teoría de invariantes computacional genuina, probablemente con bases de
Gröbner, más allá de lo que hicimos hasta ahora).

## Serie de Molien (Hilbert series del anillo de invariantes, calculable siempre)

Aunque el anillo no sea libre, su serie de Hilbert (dimensión del espacio de invariantes en cada
grado) se puede calcular con la fórmula clásica de Molien:
H(x) = (1/|G|) Σ_g 1/det(I − x·ρ(g)).

Resultado (grados 0 a 11):

| grado | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dim invariantes | 1 | 0 | 1 | 1 | 2 | 2 | 5 | 4 | 8 | 9 | 13 | 15 |

**Verificado de forma independiente** (sin confiar en la fórmula de Molien de `sympy`): el
coeficiente de grado 2 (=1) se recalculó a mano por teoría de caracteres, usando la fórmula estándar
de multiplicidad del trivial en Sym²(V) = (1/|G|)Σ_clases |clase|·½(χ(g)²+χ(g²)), con los valores de
χ(g) y χ(g²) obtenidos de la tabla de caracteres ya validada. **Coincide exacto (=1).**

Nota honesta: el grado 1 da 0, lo cual es exactamente lo esperado (una representación irreducible no
trivial no tiene vectores invariantes) — otro chequeo de consistencia gratuito.

## Qué significa esto para el proyecto (evaluación honesta del camino que sigue)

Este es un hallazgo real que **aumenta la dificultad del problema**, no la disminuye. La extensión
de n=4 a n=5 no es "repetir la misma receta con más índices" — el paso de S₃ (grupo de reflexiones)
a S₅ actuando en esta representación de 5 dimensiones (que NO es de reflexiones) cambia la
naturaleza del problema de "álgebra de anillos de polinomios libres" a "teoría de invariantes de
grupos finitos no-reflectantes", que es genuinamente más difícil (comparable, en espíritu, al salto
de "resolver por trazas cortar-y-pegar" a "bucket elimination" que tuvimos que dar en Yang-Mills).

**Próximo paso concreto** (no hecho todavía): usar la serie de Molien ya calculada para identificar
candidatos a generadores del anillo de invariantes (2 invariantes en grados bajos — 1 en grado 2, 1
en grado 3, igual que en S₃, es una coincidencia interesante en los primeros 2 grados; hay que ver
si en grado 4 aparecen 2 invariantes nuevos independientes o si son productos de los anteriores) y
verificar computacionalmente (construyendo los invariantes explícitos, no solo contando dimensiones)
si son o no algebraicamente independientes — eso revelaría dónde exactamente aparece la primera
relación no trivial.

## Archivos de esta sesión

- `scripts/anillo_invariantes_s5.py` — chequeo de reflexiones (autovalores exactos) + serie de Molien.
