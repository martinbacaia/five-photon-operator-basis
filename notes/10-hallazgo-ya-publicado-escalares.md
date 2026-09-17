# Fase 2 — Hallazgo crítico: el caso de escalares en n=5 YA ESTÁ PUBLICADO (2017)

## Qué se encontró

Antes de escribir el resultado de "S-matrices de 5 escalares" como hallazgo propio, se hizo la
búsqueda de literatura correspondiente (el método de Hilbert series para contar operadores de EFT
es estándar en la comunidad de amplitudes). Se encontró y leyó el paper central de esta técnica:

**Henning, Lu, Melia, Murayama, "Operator bases, S-matrices, and their partition functions"**,
arXiv:1706.08520 (JHEP, 2017).

La sección 5.6.1 de ese paper ("The n=5 operator basis") calcula EXACTAMENTE lo mismo que
calculamos hoy con Singular: la serie de Hilbert del anillo de invariantes de S₅ sobre la
representación de 5 dimensiones de los invariantes de Mandelstam para 5 escalares idénticos, en
dimensión genérica d≥4 (sin condiciones de Gram — el mismo régimen que asumimos nosotros). Su
ecuación (5.50a):

```
H = (1 + t⁶ + t⁷ + t⁸ + t⁹ + t¹⁵) / [(1-t²)(1-t³)(1-t⁴)(1-t⁵)(1-t⁶)]
```

con **5 invariantes primarios de grados 2,3,4,5,6** y **6 invariantes secundarios de grados
0,6,7,8,9,15** (que además reducen a 4 secundarios "irreducibles" porque η₁₅=η₆η₉=η₇η₈) —
**coincide exacto, número por número**, con lo que calculamos hoy en
`08-anillo-completo-via-singular.md`. Incluso usan `finvar.lib` de Singular con sintaxis
prácticamente idéntica a la nuestra (lo muestran explícitamente en su texto).

## Lectura de esto: 2 cosas distintas

### 1. Validación externa fuerte (la parte buena)

Que un paper publicado en JHEP (con casi 300 citas) llegue exactamente al mismo número, por un
camino de derivación que — hasta donde pudimos comparar — es distinto en la exposición pero
matemáticamente equivalente, es una confirmación externa muy fuerte de que **todo el trabajo de
esta sesión sobre S₅ está bien hecho**: la representación de S₅ sobre los invariantes de Mandelstam
(parte 2), que no es libre (parte 3), sus generadores (partes 4 y 6), y la estructura completa vía
Singular (parte 7) — todo consistente con literatura ya revisada por pares.

### 2. El caso de escalares NO es un resultado nuevo (la parte a tener en cuenta)

El propio paper de 2017 ya resuelve completamente el caso escalar de n=5 en dimensión genérica.
Publicar esto como "nuestro" resultado sería, en el mejor de los casos, redundante, y en el peor,
parecer desconocer literatura básica del campo — exactamente el tipo de error que no nos podemos
permitir dado el objetivo de publicar algo real.

## Lo que SIGUE abierto (confirmado por el propio paper de 2017)

En la sección 5.4.4 ("Generalizations with spin") del mismo paper, los autores dicen explícitamente:

> "Moving beyond scalar fields [...] Gram conditions now involve the polarization tensors in
> addition to the momenta [...] the rings for operators composed of spinning particles in general
> will *not* be free as ℂ[θ₁,...,θₘ]-modules, i.e. they will not be Cohen-Macaulay."

Es decir: **el caso con spin (fotones/gravitones) NO está resuelto en este paper** — lo dejan
explícitamente para trabajo futuro, y advierten que es estructuralmente más difícil (ni siquiera
garantizan que sea Cohen-Macaulay, una propiedad que sí usamos implícitamente para el caso escalar).//
Esto es justo el objetivo real del proyecto (la pregunta #44 de Minwalla, sobre S-matrices de n
gravitones/fotones) — sigue genuinamente abierto, y es más grande de lo que pensábamos, porque ni
siquiera el caso escalar (más simple) tenía una estructura de anillo libre para n=5, y el caso con
spin es reconocido por los propios expertos como todavía más complicado.

## Conclusión y recomendación

**No hay que escribir el caso escalar como hallazgo propio.** Todo el trabajo técnico de hoy sigue
siendo valioso — como validación del método y como base de partida real para el caso con spin — pero
hay que ser honestos en cualquier escrito futuro: la parte escalar reproduce (correctamente) un
resultado ya publicado en 2017, y el aporte, si lo hay, tiene que estar en la extensión a partículas
con spin (fotones/gravitones), que es donde el propio campo reconoce que no hay solución todavía.

**Decisión pendiente con el usuario**: dado que el caso con spin es reconocido como estructuralmente
más difícil por los propios expertos (Henning-Lu-Melia-Murayama, 2017) y no se ha intentado en n=5
en ningún lado hasta donde sabemos, hay que decidir conscientemente si seguir por ese camino (el
verdadero hueco, pero más grande) o buscar otra línea distinta dentro del campo de cuerdas.

## Fuentes

- Henning, Lu, Melia, Murayama, "Operator bases, S-matrices, and their partition functions",
  arXiv:1706.08520 (JHEP 2017) — guardado en `fuentes/operator-bases-hilbert-series.pdf`.
