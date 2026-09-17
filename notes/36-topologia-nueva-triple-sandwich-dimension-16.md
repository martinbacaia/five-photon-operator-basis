# Fase 2 — Paso 36: topología GENUINAMENTE NUEVA en grado 11 (triple sandwich) — dimensión real
# sube de 15 a 16, verificado en 2 lotes independientes

## Motivación y razonamiento previo (antes de programar nada)

Tras confirmar (parte 35) que la combinatoria completa de las 4 topologías conocidas (`F^i:F^j`,
`C^k_{ab}`, `M^{ij}_{ab}` sandwich doble) agota la dimensión en 15, se buscó una topología
GENUINAMENTE nueva, no una recombinación de lo mismo.

**Argumento de conteo (hecho a mano antes de escribir código)**: cada átomo conocido ocupa un número
fijo de "slots" de polarización (uno por partícula que usa linealmente) y tiene un grado fijo en
momentos: `F^i:F^j` (2 slots, grado 2), `C^k_{ab}` (1 slot, grado 3), `M^{ij}_{ab}` (2 slots, grado
4). Para cubrir exactamente los 5 slots de un invariante multilineal de n=5, las únicas particiones
posibles de `5` en partes `{1,2}` son `2+2+1` (con los 2-bloques siendo `F:F` o `M`, dando 3
combinaciones: `FF+FF+C`, `FF+M+C`, `M+M+C` — familias A, C, D) y `2+1+1+1` (con el 2-bloque `F:F` o
`M`: `FF+C+C+C` grado 11 exacto = familia B; `M+C+C+C` grado 13, no cabe en grado 11). **Esto agota
formalmente todas las combinaciones posibles de los 3 átomos conocidos en grado 11** — consistente
con (y explica de raíz) el resultado empírico de la parte 35.

**Átomo nuevo propuesto**: "sandwich triple" `T(i,j,k;a,b) := p_a.F^i.F^j.F^k.p_b` — trilineal en 3
polarizaciones (3 slots), grado 5 en momentos (generaliza `M` de 2 a 3 field-strengths encadenados).
Con este átomo, aparece una partición nueva: `3+1+1` (`T+C+C`, grado `5+3+3=11` EXACTO, sin
corrector — el tipo de partición "exacta" que ya había resultado ser la más rica en A/D).

## Implementación y verificación

Se implementó `T` propagando un vector a través de `F^i, F^j, F^k` sucesivamente (misma lógica que
`Msandwich`, generalizada): `v → v.F^m = (v.p_m)·ε_m − (v.ε_m)·p_m`, representando el vector
únicamente por sus productos punto con todos los momentos/polarizaciones (nunca componentes
explícitas) — evita cualquier elección de base o gauge.

**Verificación 0 (obligatoria antes de confiar en la cadena)**: la MISMA cadena genérica, aplicada 2
veces, reproduce exactamente `Msandwich` ya validado en partes anteriores — coincide en todos los
`(i,j,a,b)` probados. Esto valida el mecanismo de propagación antes de usarlo para un átomo nuevo.

**Verificación de no-anulación** (lección de la parte 13: varias topologías se anulan idénticamente
por antisimetría al simetrizar sobre S₅ — no aceptar nada sin este chequeo): se probaron 6
variantes de `T(i,j,k;a,b)·C(l;c,d)·C(m;e,f)` con distintas elecciones de patas de momento. **3 de
6 se anulan idénticamente** (mismo mecanismo de siempre: antisimetría de `C` combinada con una
permutación de la órbita S₅) — **3 sobreviven**, no nulas, con la variante `a,b=(l,m)`,
`(c,d)=(i,j)`, `(e,f)=(i,k)` dando el valor más limpio.

## Resultado: dimensión sube de 15 a 16

Se construyó una familia representativa de 90 semillas (`grado11_familia_T_triple_sandwich.py`:
10 elecciones de patas `(a,b)` × 3×3 elecciones de los pares de `C` dentro de `{i,j,k}`) y se
combinó con las filas YA CALCULADAS de A+B+C+D (`combinatoria_completa_rows.pkl`, mismas 20
muestras cinemáticas — reduce cómputo reutilizando trabajo ya hecho).

- **Lote 1** (muestras cinemáticas 1-20): rango combinado (A+B+C+D+T) = **16** (antes 15 sin T),
  residual `4.4×10⁻⁶¹`.
- **Lote 2** (muestras 201-220, independiente): rango combinado = **16** también, residual
  `1.2×10⁻⁶¹`.

**La topología T agrega exactamente 1 dimensión nueva, confirmado en 2 lotes independientes con
alta precisión (60 dígitos, residual ~10⁻⁶¹).** No es un artefacto de muestreo — reproducido idéntico
con cinemáticas completamente distintas.

## Estado honesto — qué falta

1. Sólo se probó una familia REPRESENTATIVA de 90 semillas de la topología T (variando patas de
   momento, pero restringiendo los pares de `C` a los 3 dentro de `{i,j,k}`, no los 10 pares
   posibles entre los 5 slots). **No se descartó que T aporte más de 1 dimensión** — sólo se probó
   que aporta AL MENOS 1. Ampliar la combinatoria de T (análogo a lo hecho para C/D en la parte 35)
   es el paso obvio siguiente.
2. No se probaron las otras 2 particiones nuevas que T habilita (`T+F:F+corrector`,
   `T+M+corrector`, ambas grado 11 con corrector Mandelstam) — podrían aportar dimensiones
   adicionales.
3. No se evaluó todavía esta nueva dirección en ninguna de las 3 familias de límites de Regge (partes
   22-23, 32, 34) — es el paso que realmente importa para el objetivo del proyecto (¿ayuda a
   satisfacer el bound `s²`?). Dado el patrón ya visto (partes 30-34), no hay que apresurarse a
   evaluar esto antes de entender bien la dimensión completa que aporta la topología T.

## Actualización (misma sesión, continuación): dimensión de T totalmente acotada — sigue en +1

Se ampliaron los 3 chequeos pendientes:
1. **Familia T, combinatoria completa de pares de `C`** (los 10 pares posibles entre los 5 slots,
   no sólo los 3 dentro de `{i,j,k}`, con `a,b` fijo en el mejor caso encontrado): 100 semillas
   adicionales — **el rango se queda en 16**, no sube más dentro de esta rama.
2. **Partición `T+F:F+corrector`** (grado `5+2+4=11`): probadas variantes, 1 de 4 no nula. Familia
   representativa de 60 semillas.
3. **Partición `T+M+corrector`** (grado `5+4+2=11`): probadas variantes, 2 de 3 no nulas. Familia
   representativa de 160 semillas.

**Resultado combinado final** (A+B+C+D + las 4 familias de `T`, 1297 semillas totales, 794 no
nulas): **rango = 16 exacto**, mismo resultado que con sólo la familia T original (90 semillas) —
ninguna de las ampliaciones (más pares de `C`, ni las otras 2 particiones que habilita `T`) aportó
una dimensión adicional. **La contribución del átomo "sandwich triple" está ahora bien acotada: +1
dimensión, no más, dentro de lo explorado.**

**Nota de infraestructura** (relevante para el resto del proyecto): el cómputo combinado de las
familias `T+F:F` y `T+M` (220 semillas) se atascó 2 veces corriendo como un solo proceso largo
(CPU casi nula durante minutos, sin terminar ni fallar — mismo patrón ya visto en la parte 35).
Diagnosticado: NO es un seed específico patológico (verificado por separado, cada familia sola
corre en <1 minuto por muestra) — es un problema de acumulación/scheduling en procesos largos de
este entorno. **Solución que funcionó**: dividir en llamadas foreground cortas (5 muestras por
llamada) con guardado incremental (`familia_TFFTM_incremental.py`) — terminó limpio en 4 llamadas
de ~5-7 minutos cada una. Reafirma la lección de la parte 35: evitar depender de un solo proceso
largo desatendido en este entorno.

## Próximo paso concreto

1. Repetir el análisis de rango/núcleo en las 3 familias de Regge ya construidas
   (`regge_n5_simbolico_grado11.py`, `regge_n5_familia2_grado11.py`, `fast_grado11_familia3.py`)
   extendidas para incluir la dirección nueva de `T`, para ver si el subespacio universal deja de
   ser `{0}` al contar con ella — es el paso que conecta este hallazgo con el objetivo real del
   proyecto (CRG).
2. Si eso no alcanza, buscar UN átomo genuinamente distinto de `T` (no una ampliación de la misma
   idea) — candidatos: Levi-Civita, o estructuras con 4+ field-strengths encadenados (aunque el
   patrón de esta parte sugiere que "más de lo mismo" tiende a saturar rápido).

## Archivos

- `scripts/topologia_nueva_triple_sandwich.py` — define y verifica el átomo `T` (sandwich triple).
- `scripts/grado11_familia_T_triple_sandwich.py` — familia representativa de 90 semillas, combinada
  con A+B+C+D, rango 16 en 2 lotes independientes.
- `scripts/familia_T_rows.pkl` — filas calculadas de la familia T (20 muestras, lote 1).
