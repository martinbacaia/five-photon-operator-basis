# Fase 2 — Paso 12: cambio de estrategia (Camino B) — trabajar directamente con contracciones de
# field-strength en vez de coordenadas (α,α'), evitando por completo el problema del pivote/leftover

## Motivación

Las partes 9-11 mostraron que reducir los datos de polarización de fotones de n=5 a un par de
coordenadas explícitas (α₁,α₁') introduce una dependencia espuria del "pivote" elegido — con una
obstrucción demostrada (no solo no encontrada) para trivializarla en general (ver
`15-renormalizacion-y-doublecheck.md`). Esta sección documenta un cambio de estrategia: **evitar
completamente la necesidad de coordenadas explícitas**, trabajando con las contracciones de
field-strength directamente — que resultan ser honestamente covariantes bajo S₅ sin ningún ajuste.

## Hallazgo 1: C^i_jk := p_j.F^i.p_k es perfectamente S₄-covariante, sin ambigüedad de signo ni pivote

Se verificó — simbólica y numéricamente, con la cinemática explícita D=6 ya construida en
`kinematics_numeric_check.py` — que la familia completa de contracciones
`C^1_{jk} := p_j.F¹.p_k` (para {j,k}⊂{2,3,4,5}, 6 valores) transforma bajo permutaciones de
{2,3,4,5} **exactamente como la relabeling directa de sus índices**, sin ningún factor de
ambigüedad:

```
C^1_new_{j,k} = C^1_{sigma(j),sigma(k)}   (evaluado en la configuracion vieja)
```

verificado con vectores físicos concretos (diferencia 0.00e+00 en los 6 pares probados bajo la
transposición (2 3)). Esto es fundamentalmente distinto de (α₁,α₁'), que requieren *elegir* 2 de
los 6 pares como base — ese paso de elección es precisamente el origen de toda la dificultad de
las partes 9-11. C^i_jk, al no requerir elección alguna, no tiene ese problema.

## Hallazgo 2: las 4 relaciones lineales entre los 6 C_jk, calculadas explícitamente

Dado que el espacio físico es de rango 2, los 6 `C_jk` (para partícula 1, pares en {2,3,4,5}) están
sujetos a 4 relaciones lineales independientes, con coeficientes que son funciones racionales
simples de los invariantes de Mandelstam (calculadas resolviendo el sistema 2×2 de
`13-formula-inversa-polarizacion-n5.md` y sustituyendo):

```
C_25 = -C_23 - C_24
C_34 = -(s14/s12) C_23 + (s13/s12) C_24
C_35 = ((s12+s14)/s12) C_23 - (s13/s12) C_24
C_45 = -(s14/s12) C_23 + ((s12+s13)/s12) C_24
```

(dᵢⱼ:=pᵢ.pⱼ, sᵢⱼ:=−2dᵢⱼ). Estas relaciones son la versión "sin elegir base" de la información que
antes intentábamos capturar invirtiendo a (α₁,α₁') — **no requieren pivote, no tienen ambigüedad**,
y describen completamente el subespacio físico de 2 dimensiones dentro del espacio de 6
antisimétricos, para cualquier cinemática dada.

## Hallazgo 3 (el más importante): el building block correcto para combinar 2 partículas distintas
## es F^i:F^j (contracción completa de los 2 tensores field-strength), NO ε_i^∥.ε_j^∥

Al buscar cómo construir invariantes que combinen los datos de polarización de **2 partículas
distintas** (necesario para cualquier S-matriz real, que debe ser lineal en CADA εᵢ), se descartó
la idea ingenua de usar ε_i^∥.ε_j^∥ directamente: **no es gauge invariante** en general (contiene
una pieza proporcional a a_i·(p_i.ε_j^∥), que no se cancela salvo que ε_j^∥ sea ortogonal a p_i,
lo cual no está garantizado).

**La cantidad correcta es la contracción completa de los 2 tensores field-strength**:

```
F^i : F^j := F^i_{mu nu} F^j^{mu nu} = 2[(p_i.p_j)(eps_i.eps_j) - (p_i.eps_j)(eps_i.p_j)]
```

Esta cantidad es **manifiestamente gauge invariante en ambas polarizaciones simultáneamente**
(F^i es invariante bajo εᵢ→εᵢ+ζᵢpᵢ por construcción — ya usado toda la sesión —, y lo mismo para
F^j; la contracción de 2 objetos gauge-invariantes es gauge invariante). **Verificado numéricamente**
con la cinemática explícita: se calculó F¹:F² con polarizaciones ε₁,ε₂ aleatorias (sujetas solo a
εᵢ.pᵢ=0), y se comprobó que un shift de gauge arbitrario simultáneo (ε₁→ε₁+ζ₁p₁, ε₂→ε₂+ζ₂p₂) deja
F¹:F² **invariante a precisión de máquina** (diferencia ~7e-15).

**Por qué esto resuelve el problema de fondo**: F^i:F^j no requiere elegir NINGÚN pivote ni base —
es una contracción de Lorentz completamente estándar entre 2 tensores ya gauge-invariantes. No hay
ninguna elección arbitraria escondida, y por tanto no puede aparecer ninguna dependencia espuria de
cinemática del tipo encontrado en las partes 9-11.

**Estructura**: expandiendo,

```
F^i:F^j = 2(p_i.p_j)(eps_i^perp.eps_j^perp) + 2[(p_i.p_j)(eps_i^par.eps_j^par) - (p_i.eps_j^par)(eps_i^par.p_j)]
```

el primer término es el análogo directo del building block "ε_i^⊥.ε_j^⊥" que ya usa el paper
original para n=4 (sección 3), y el segundo término (el que involucra las partes ∥, o sea los datos
tipo α) es estructuralmente del mismo tipo que los C^i_{jk} ya estudiados — de hecho se puede
verificar que reduce a una combinación de C-type contracciones. Falta (próxima sesión): derivar la
forma cerrada explícita de este segundo término en términos de (α_i,α_i',α_j,α_j') y los invariantes
de Mandelstam, generalizando el cálculo ya hecho para C_jk individuales.

## Camino a seguir (próxima sesión, concreto)

1. Derivar la forma cerrada de la parte "∥" de F^i:F^j en términos de contracciones C-type (sin
   pasar por α,α' explícitos si es posible, o al menos con una convención de pivote fija POR
   PARTÍCULA, ya que F^i:F^j en sí mismo no tiene el problema de pivote — solo lo tendría si
   insistimos en reexpresarlo en términos de α's explícitos, lo cual ya no es necesario).
2. Repetir para F^i:F^j con i,j variando sobre las 5 partículas, y para el análogo de grado más
   alto en momentos (contracciones con más potencias de p, para construir toda la torre de
   S-matrices polinomiales, como hace el paper en su sección 5 para n=4).
3. Construir combinaciones S₅-invariantes (sumas sobre el grupo, o proyecciones a la
   representación trivial) de productos de estas F^i:F^j (y de invariantes de Mandelstam, cuyo
   anillo completo ya se calculó en las partes 3-8) — esto da directamente candidatos a
   S-matrices de fotones para n=5, sin necesitar resolver el problema de pivote de las partes 9-11
   en absoluto.
4. Verificar la condición de Regge growth (el objetivo final del proyecto) sobre estos candidatos.

**Evaluación honesta**: este cambio de estrategia parece resolver la dificultad estructural
encontrada en las partes 9-11 (el problema no era la física, era la elección de coordenadas). Queda
bastante trabajo de construcción explícita (puntos 1-4), pero la vía ahora parece clara y sin
obstrucciones conocidas.

## Archivos

- Verificaciones de esta sección hechas inline (extendiendo `kinematics_numeric_check.py`); no se
  generó un script nuevo separado, dado que son chequeos cortos. Cálculo simbólico de las 4
  relaciones lineales entre C_jk: ver comandos en el historial de esta sesión (reutiliza el
  armazón de `formula_inversa_n5.py`).
