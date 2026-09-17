# Fase 2 — Paso 6: grados 5 y 6 del anillo de invariantes — hace falta más generadores de lo previsto

## Resultado

- **Grado 5** (Molien predice dimensión 2): el producto Q·C (de los generadores de grado 2 y 3)
  aporta 1 dimensión, pero **no alcanza** — se confirmó con 3 semillas de Reynolds distintas, cada
  una genuinamente independiente de Q·C. **Hace falta un generador nuevo en grado 5** (llamado E₅
  en la documentación, sin forma cerrada aislada todavía).
- **Grado 6** (Molien predice dimensión 5): los productos de generadores ya conocidos
  {Q³, C², E₄·Q} solo alcanzan dimensión 3 (verificado, rango exacto = 3 — no hay forma de usar E₅,
  que es de grado 5, para llegar a grado 6 sin un generador de grado 1, que no existe). **Faltan 2
  dimensiones más**, que se confirmó están genuinamente presentes agregando 5 semillas nuevas de
  Reynolds — el conjunto completo (8 candidatos) alcanza **exactamente** la dimensión 5 predicha
  por Molien.

**Conclusión honesta**: el anillo de invariantes de S₅ sobre este espacio de 5 dimensiones necesita
generadores en (al menos) los grados 2, 3, 4, 5, y aparentemente 2 más en grado 6 — bastante más
rico que el caso de S₃ (que se agota en 2 generadores, grados 2 y 3). Esto es consistente con lo ya
establecido en el paso 3 (`03-anillo-invariantes-s5.md`): al no ser un grupo de reflexiones en esta
representación, no hay garantía de que el número de generadores sea chico ni de que aparezcan en
pocos grados.

## Un TERCER bug real, encontrado y corregido (nivel de rigor pedido explícitamente por el usuario)

Al calcular el rango del conjunto {Q³, C², E₄·Q} + 5 semillas nuevas en grado 6, el primer intento
dio **rango 6** — que es matemáticamente imposible si Molien dice que la dimensión real es 5 (no
puede haber 6 vectores independientes en un espacio de dimensión 5). En vez de aceptar el número, se
investigó sistemáticamente:

1. Se re-confirmó el coeficiente de grado 6 de Molien por un método completamente fresco (suma
   numérica de autovalores sobre los 120 elementos, sin tabla de caracteres) — dio exactamente 5.
2. Se verificó que los 8 candidatos SÍ eran genuinamente invariantes bajo los 120 elementos (0
   fallos) — descartando que el problema fuera falta de invariancia real.
3. Se verificó el cálculo del rango en sí (auto-test con casos triviales conocidos, más
   confirmación vía forma escalonada reducida — RREF con 6 pivotes) — la función de cálculo de
   rango estaba bien implementada.
4. **Causa real encontrada**: una de las 5 semillas nuevas (`a₁²·a₂·a₃·a₄`) tiene grado **5**, no
   6 — un error de tipeo simple al escribir la lista de monomios semilla (faltaba un factor `a₅`
   para llegar a grado 6). Al mezclar un invariante de grado 5 genuino con invariantes de grado 6
   genuinos en el mismo cálculo de rango, aparece como "independiente" — pero solo porque vive en
   un subespacio de monomios completamente disjunto (grado 5 vs. grado 6), no porque aporte una
   dimensión real al espacio de invariantes de grado 6 que nos interesaba. Esto infla el rango en
   exactamente 1 de forma espuria — coincide perfecto con el 6 encontrado vs. el 5 correcto.

**Corregido**: se arregló la semilla (`a₁²·a₂·a₃·a₄·a₅`, grado 6 genuino) y además se agregó un
**guardrail permanente** a la función `rank_of_span()`: ahora exige explícitamente que todos los
polinomios pasados sean homogéneos del mismo grado, y frena con un error claro si no lo son — en vez
de devolver silenciosamente un número que mezcla grados distintos. Tras la corrección, el resultado
es exactamente 5, coincidiendo con Molien.

**Tercera lección metodológica de esta sesión** (después del anti-homomorfismo y el mal uso de
`nsimplify`): al construir listas de "semillas" a mano para probar independencia lineal, verificar
SIEMPRE la homogeneidad de grado de cada semilla antes de mezclarla con otras — un error de conteo
de grado al escribir un monomio a mano puede inflar un resultado de forma silenciosa y
convincente (el número resultante "parece" razonable, no da un error obvio).

## Estado del conteo de generadores (actualizado, honesto)

| grado | dim. Molien | cubierto por productos de generadores previos | generador(es) nuevo(s) necesario(s) |
|---|---|---|---|
| 2 | 1 | — | Q (grado 2) |
| 3 | 1 | — | C (grado 3) |
| 4 | 2 | 1 (Q²) | 1 nuevo (E₄) |
| 5 | 2 | 1 (Q·C) | 1 nuevo (E₅) |
| 6 | 5 | 3 (Q³, C², E₄·Q) | 2 nuevos (sin nombrar todavía) |

## Evaluación honesta del camino que sigue

Encontrar el conjunto COMPLETO de generadores y las relaciones (la "presentación" completa del
anillo) por este método de "adivinar semillas grado por grado" es viable pero cada vez más costoso
computacionalmente (los grados 7+ requieren manejar polinomios con cientos de términos). Para hacer
esto de forma sistemática y confiable —especialmente si el objetivo es publicar un resultado
correcto— **el paso más razonable ahora es usar software de teoría de invariantes computacional
dedicado** (ej. paquetes de invariant theory en Macaulay2, Singular, Magma, o el paquete
`sympy`/`sage` de invariant rings si existe con soporte adecuado para grupos de permutación no
reflectantes) en vez de seguir escribiendo semillas a mano — el riesgo de errores humanos tipo
"semilla de grado equivocado" crece con la complejidad, como ya se vio.

## Próximo paso concreto (no hecho, para decidir con el usuario)

Evaluar si conviene:
(a) seguir con el método artesanal (Reynolds + rank_of_span, ya con el guardrail de grado) unos
    grados más, o
(b) invertir tiempo en instalar/aprender una herramienta de teoría de invariantes computacional
    dedicada (más robusta para esta escala de problema, pero con curva de aprendizaje).

## Archivos de esta sesión

- `scripts/generadores_grado5_6.py` — corregido, con guardrail de homogeneidad de grado agregado.
