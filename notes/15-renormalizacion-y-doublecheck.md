# Fase 2 — Paso 11: double-check con cinemática numérica explícita + renormalización que
# trivializa 3 de las 4 transiciones de leftover (la 4ª, con cambio de pivote, resiste — con
# obstrucción demostrada, no solo "no encontrada")

## Parte A: double-check de todo lo hecho hasta ahora, con cinemática explícita (no solo álgebra formal)

Todo el trabajo de las partes 9-10 (`13-formula-inversa-polarizacion-n5.md`,
`14-accion-s5-alpha-n5.md`) se verificó hasta ahora **solo** tratando los productos punto pᵢ.pⱼ
como símbolos formales sujetos a las 5 ecuaciones de conservación de momento — nunca se había
comprobado que estos símbolos corresponden genuinamente a momentos D-dimensionales reales
(vectores concretos). Esto es una brecha real: un error en el paso de "resolver conservación de
momento" (sección con `sp.solve`) podría no manifestarse en los chequeos puramente simbólicos si
fuera autoconsistente pero no correspondiera a cinemática física real.

**Se cerró esta brecha.** Se construyó cinemática explícita de 5 momentos no nulos en D=6
(`p_i^2=0`, `Σp_i=0`, métrica mostly-plus) usando un método numérico: 4 momentos nulos
aleatorios (energía + dirección espacial unitaria aleatoria) más un 5º momento cuya energía se
ajusta (`scipy.optimize.brentq`) para que `p1=-(p2+p3+p4+p5)` sea también nulo — verificado
`p_i²=0` para los 5 y `Σp_i=0` a precisión de máquina (~1e-16).

Con esta cinemática real (y un vector de polarización ε₁ aleatorio sujeto solo a ε₁.p₁=0), se
verificó, recalculando **desde cero con vectores concretos** (nunca reutilizando las fórmulas
simbólicas ya derivadas, solo la definición cruda `F¹_μν=p1_μillon ε1_ν−p1_ν ε1_μ`):

1. **Las variables "dependientes" resueltas simbólicamente** (d15,d25,d34,d35,d45 en términos de
   las 5 libres) **coinciden exactamente** con los productos punto reales calculados directamente
   de los vectores — confirma que la resolución de conservación de momento (usada en TODOS los
   pasos previos de la sesión) es correcta y corresponde a cinemática físicamente realizable, no
   solo a álgebra formal autoconsistente.
2. **La fórmula ε₁^∥ = α₁w₂₃+α₁'w₂₄+a₁p₁** (parte 12): se reconstruyó ε₁^∥ con los α's obtenidos
   vía la fórmula inversa (parte 13) y se comparó con ε₁ real componente a componente (vía
   ε₁.pₘ para m=2,3,4,5) — la diferencia es **exactamente proporcional a pᵢ.pₘ con la MISMA
   constante a₁** en los 4 casos (0.272845... en el punto probado) — confirma que la diferencia es
   pura gauge (dirección p₁) y nada más, tal como predice la fórmula.
3. **La representación 2M de S₃** (parte 14, v2): se recalcularon α₁,α₁' relabeleando físicamente
   los vectores momento (no solo sus símbolos) para las permutaciones (2 3) y el 3-ciclo (2 3 4), y
   coinciden **exactamente** (diferencia ~1e-15) con lo predicho por las matrices ψ(σ) ya
   tabuladas.
4. **La función de transición T(5→4)** (parte 14): recalculada con leftover=4 usando los mismos
   vectores físicos, coincide exactamente con la predicción `T·(α₁,α₁')`.

**Conclusión de la parte A**: los 4 resultados centrales de las partes 12, 13 y 14(v2) están
confirmados de forma independiente y más fuerte que antes — con vectores físicos concretos, no
solo álgebra simbólica formal. No se encontró ningún error.

## Parte B: búsqueda de la renormalización que trivialice las transiciones de leftover

### Resultado positivo: Φₓ := αₓ/s₁ₓ trivializa las 3 transiciones "mismo pivote"

Para las 3 elecciones de leftover que comparten pivote=2 (leftover∈{3,4,5}), se derivó (a mano,
usando la identidad general `w_{2,m} = −(s₁ₐ/s₁ₘ)w_{2,a} − (s₁_b/s₁ₘ)w_{2,b} + (1/s₁ₘ)p₁` para
cualquier terna {a,b,m}⊂{3,4,5}) la transición general entre "leftover=z" (compañeras x,y) y
"leftover=y" (compañeras x,z):

```
alpha_x^(y) = alpha_x^(z) − (s1x/s1y) alpha_y^(z)
alpha_z^(y) = −(s1z/s1y) alpha_y^(z)
```

**Definiendo Φₓ := αₓ/s₁ₓ (dividir cada componente por el invariante s₁ del PARTÍCULA a la que
corresponde, sin importar en qué base/leftover aparece), la transición se vuelve**:

```
Phi_x^(y) = Phi_x^(z) − Phi_y^(z)
Phi_z^(y) = −Phi_y^(z)
```

**una matriz literalmente constante `[[1,−1],[0,−1]]`, la MISMA para cualquier elección de
{x,y,z}⊂{3,4,5}** — verificado simbólicamente (identidad exacta, ver script) y numéricamente con
la cinemática explícita de la parte A (coincide a 1e-15 para las 3 transiciones entre leftover=3,4,5).

**Esto es un resultado limpio y nuevo**: análogo estructural directo del factor √(st/u) del paper
para n=4 — una renormalización kinematics-dependiente que convierte una construcción básica
(pivote fijo) en una cantidad que transforma con una matriz constante bajo el subgrupo relevante.

### Resultado negativo (con demostración, no solo intento fallido): no se puede extender a leftover=2

La 4ª transición (hacia leftover=2, que **obliga a cambiar el pivote** de 2 a 3, porque "2" es
precisamente la partícula que se deja afuera) tiene una estructura distinta:

```
T(pivote2,companeras(3,4) -> pivote3,companeras(4,5)) = (1/s12) * [[s14, s12+s14], [s15, s15]]
```

(se usó s₁₂+s₁₃+s₁₄+s₁₅=0 para simplificar). **Se probó, con un cálculo simbólico exhaustivo (no
solo prueba y error), que NINGÚN reescalado monomial de la forma**

```
Phi_x^(L) := alpha_x^(L) * s1(pivote_L)^a / s1x^b
```

**(para ningún par de exponentes a,b, incluyendo a=−1,b=1 que sí funciona para las otras 3
transiciones) puede volver esta matriz constante.** Razón estructural identificada: la entrada
(1,2) de T es `(s₁₂+s₁₄)/s₁₂` — una **suma de un término homogéneo de grado 0 y otro no
homogéneo** (1 + s₁₄/s₁₂) que ningún reescalado diagonal monomial puede simplificar a una
constante (un reescalado monomial multiplica toda la entrada por un factor único; no puede
"partir" una suma de grados distintos). Se verificó exigiendo que el coeficiente del primer término
(Φₓ) sea constante — lo cual fija unívocamente a=−1,b=1 — y comprobando que, con esos valores YA
FIJOS (sin más libertad), el coeficiente del segundo término (Φ_y) queda `(s₁₂+s₁₄)/s₁₃`,
genuinamente dependiente de la cinemática.

**Esto no es "no lo encontramos" — es una obstrucción demostrada** dentro de la familia de
renormalizaciones diagonales/monomiales más natural. Quedan 2 caminos abiertos para la próxima
sesión:

1. Buscar una renormalización **no diagonal** (mezclando Φ₃,Φ₄ con una matriz general, no solo un
   reescalado por componente) — existe en principio (el "fibrado" es trivial topológicamente,
   al vivir sobre solo 4 puntos discretos), la pregunta es si existe una con forma cerrada simple.
2. **Abandonar la búsqueda de una base explícita para (α,α') y trabajar directamente con las
   contracciones C_jk = p_j.F¹.p_k (o el propio tensor F¹) como generadores del módulo** — evitando
   por completo la necesidad de "invertir" a coordenadas α, que es precisamente el paso que
   introdujo la dependencia de pivote/leftover. Esta es la estrategia que ya se había anotado como
   opción 1 en `14-accion-s5-alpha-n5.md`, y esta sesión aporta evidencia de que puede ser la más
   prometedora, dado que la dificultad parece inherente a la elección de coordenadas, no a la
   física.

## Archivos

- `scripts/kinematics_numeric_check.py` *(construido inline en esta sesión, ver comandos)* —
  construcción de cinemática explícita D=6, verificación de las partes 12/13/14 con vectores reales.
- Cálculos de la parte B (renormalización Φₓ, prueba de obstrucción) hechos inline; ver esta
  sección para las fórmulas exactas — no dejan un script separado nuevo, están documentados aquí
  íntegramente porque son derivaciones cortas y autocontenidas.
