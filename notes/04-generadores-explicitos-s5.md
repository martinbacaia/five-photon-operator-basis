# Fase 2 — Paso 4: generadores explícitos del anillo de invariantes de S₅ (grados 2, 3 y 4)

## Qué se hizo

Se construyeron **explícitamente** (no solo se contaron dimensiones) los invariantes de bajo grado
del anillo de S₅ sobre el espacio de 5 invariantes de Mandelstam, usando el **operador de Reynolds**
R(f) = (1/|G|) Σ_g f(ρ(g)·a) — promediar un polinomio sobre los 120 elementos del grupo, generados
por multiplicación de matrices a partir de 2 generadores ((12) y (12345)), con la homomorfía
verificada de nuevo por muestreo aleatorio antes de confiar en el resultado.

**Encontrado**:
- Grado 2: un invariante Q(a) no nulo (coincide con la única dimensión que predice Molien).
- Grado 3: un invariante C(a) no nulo (ídem).
- Grado 4: se generaron 3 candidatos (Q², y 2 vía Reynolds sobre monomios distintos) — ver la
  sección del bug de abajo antes de tomar esto por bueno.

Ver `scripts/generadores_invariantes_s5.py`.

## Un segundo bug real, encontrado y corregido en esta misma sesión

Al chequear si los 3 candidatos de grado 4 eran linealmente independientes, se usó primero
evaluación numérica en puntos aleatorios racionales (con `sympy.nsimplify` sobre los valores ya
evaluados) — dio **rango 3**, lo cual contradecía las 4 fórmulas de Molien ya calculadas y
verificadas (que dicen dimensión 2 en grado 4). En vez de ignorar la contradicción, se investigó:

1. Se confirmó por fuerza bruta que los 3 candidatos SÍ son invariantes bajo los 120 elementos del
   grupo, sin excepción (0 fallos) — descarta que el problema fuera falta de invariancia genuina.
2. Se recalculó la dimensión en grado 4 por un CUARTO método independiente: suma numérica directa
   (autovalores vía `numpy`) sobre los 120 elementos del grupo, sin usar tabla de caracteres ni
   representantes de clase — **dio 2 de nuevo**, reforzando que Molien estaba bien.
3. Se abandonó la evaluación en puntos aleatorios (que usaba `nsimplify`, riesgoso sobre valores ya
   exactos) y se resolvió el problema de independencia lineal de la forma correcta: comparación
   **simbólica exacta, coeficiente a coeficiente**, de la combinación c₁Q²+c₂D₁+c₃D₂=0.

**Resultado correcto**: SÍ existe una relación lineal exacta:
`D₁ = 10·Q² − 9·D₂` (con D₁, D₂ los 2 candidatos obtenidos vía Reynolds). Es decir, el espacio real
es de dimensión 2 (no 3) — **coincide exactamente con Molien**. El error estaba en el método de
verificación (uso indebido de `nsimplify` sobre puntos ya exactos), no en la matemática de fondo.

**Lección metodológica** (segunda de esta sesión, después de la del anti-homomorfismo): verificar
independencia lineal de polinomios exactos por evaluación numérica en puntos aleatorios es
razonable, pero pasar los valores por `nsimplify` (pensado para "adivinar" formas cerradas a partir
de floats, no para limpiar racionales exactos) puede introducir errores silenciosos. La comparación
exacta de coeficientes es la forma correcta cuando ya se trabaja con expresiones simbólicas exactas.

## Conclusión sobre los generadores (después de la corrección)

El espacio de invariantes de grado 4 (dimensión 2, confirmado por 4 métodos) se descompone como:
- 1 dimensión = Q² (el cuadrado del generador de grado 2 — un "descendiente", no un generador nuevo).
- 1 dimensión genuinamente **nueva** (ninguna combinación de Q y C, los únicos generadores
  conocidos hasta grado 3, puede llegar a grado 4 salvo Q² — no hay invariante de grado 1 para
  combinar con C de grado 3).

**Esto confirma y hace concreto** lo que ya se sabía cualitativamente desde `03-anillo-invariantes-s5.md`
(el anillo no es libre en el sentido de S₃): el anillo de invariantes de S₅ sobre este espacio de
5 dimensiones necesita **al menos 3 generadores** (grados 2, 3 y 4) — a diferencia de S₃ sobre
(s,t,u), que necesita solo 2 (grados 2 y 3, sin nada nuevo después). El próximo generador nuevo
esperado, según la serie de Molien (grados 0-11: 1,0,1,1,2,2,5,4,8,9,13,15), habría que buscarlo
donde la dimensión salte más de lo que explican los productos de generadores ya conocidos — grado 5
da dimensión 2 (que ya se explica por Q·C, un producto de los generadores de grado 2 y 3 — hay que
verificar si es 1-dimensional o si también hace falta ahí un generador nuevo) y grado 6 da 5
(bastante más alto, candidato fuerte a otro generador nuevo o a la primera RELACIÓN entre los ya
encontrados).

## Próximo paso concreto (no hecho)

1. Verificar en grado 5 si Q·C (producto de los generadores de grado 2 y 3) agota la dimensión 2
   predicha por Molien, o si hace falta otro generador nuevo ahí también.
2. En grado 6 (dimensión 5 según Molien), contar cuántas dimensiones cubren los productos de los
   generadores ya encontrados (Q³, C², y el nuevo generador de grado 4 combinado con Q) — la
   diferencia con 5 revelaría si hace falta un generador nuevo en grado 6, o si en cambio aparece la
   primera RELACIÓN (cuando los productos ya conocidos exceden la dimensión real del espacio).

## Archivos de esta sesión

- `scripts/generadores_invariantes_s5.py` — construcción de Q, C, y candidatos de grado 4 vía
  Reynolds, más la verificación (corregida) de independencia lineal.
