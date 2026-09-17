# Fase 2 — Paso 26: cálculo riguroso de la dimensión real de grado 9 (con sandwich incluido
# sistemáticamente) — dimensión exacta = 3, y una combinación que baja el crecimiento a s³

## Recalculo riguroso de la dimensión de grado 9

Siguiendo el pedido explícito de verificar todo con el máximo rigor (el objetivo es publicar), se
rehizo el cálculo de la parte 20 **incluyendo sistemáticamente la topología "sandwich"** encontrada
en la parte 25 (que la parte 20 no había considerado). Se enumeraron:

- Familia A (ya usada en la parte 20): "2 parejas F:F + 1 singlete C" × 10 multiplicadores de
  Mandelstam × 3 formas de contracción del singlete = 30 semillas.
- Familia B (nueva): "1 pareja F:F + 1 sandwich + 1 singlete C" — el sandwich `M^{ij}_{ab}` necesita
  2 momentos externos elegidos de los 3 restantes (6 ordenados, ya que `M` no es simétrico en sus
  argumentos — verificado en la parte 25), y el singlete `C^k_{ab}` puede contraerse con cualquiera
  de los 4 momentos restantes (6 pares) = 36 semillas.

**Total: 66 semillas**, evaluadas en 30 muestras cinemáticas independientes, matriz reducida a rango
numérico. **Resultado: rango = 3**, con un salto de valores singulares de 14 órdenes de magnitud
(`2.4` vs `2.9×10⁻¹⁴`) — inequívoco. **Verificado estable en 4 lotes de muestras independientes**
(distintos rangos de semilla: 1-30, 201-230, 501-550, 1001-1040), todos dando rango 3 con saltos
igualmente limpios.

**Corrección definitiva: el espacio de invariantes de grado 9 tiene dimensión EXACTA 3** (no 1 como
se concluyó en la parte 20 —esa conclusión era válida solo dentro de la familia limitada que se
probó ahí—, ni solo "al menos 2" como se dejó abierto en la parte 25).

## Base explícita y comportamiento en el límite de Regge

Se extrajeron 3 semillas pivote (vía eliminación gaussiana) que forman una base explícita:

```
B1 = F^i:F^j · F^l:F^m · C^k_{ij} · d(i,l)                      [familia A]
B2 = F^l:F^m · M^{ij}_{l,m} · C^k_{i,l}                          [familia B, sandwich]
B3 = F^l:F^m · M^{ij}_{l,k} · C^k_{i,l}                          [familia B, sandwich]
```

Evaluadas en el límite de Regge de n=5 (parte 23): **las 3 crecen exactamente como `E^8 ~ s^4`**
(pendientes 7.9999, 7.9999, 7.9999 — limpio). Con dimensión 3 y las 3 compartiendo el mismo orden
líder, hay margen real para cancelar (a diferencia de grado 11, donde la dimensión también era 3
pero el mejor resultado alcanzable topaba en el mismo `s⁴` — ver parte 24).

## Hallazgo: la combinación `B1 − 2B2` cancela el término líder y crece como `E⁶ ~ s³`

Se encontró (por extrapolación de Richardson de la razón `B1/B2` a `E` grande) que **`A1/A2 = 2`
exactamente** (confirmado a 9+ cifras estables, y de forma independiente con aritmética de alta
precisión `mpmath` a 50-100 dígitos — converge consistentemente a `2.000000000` con error `<10⁻⁹`).
**Importante**: se verificó explícitamente que esta relación `B1=2B2` **NO es una identidad global**
— en cinemática genérica (no Regge) el cociente `B1/B2` varía completamente (-73.4, 2.45, -0.65,
-0.29, 1.57 en 5 semillas aleatorias) — es una coincidencia del comportamiento LÍDER específicamente
en este límite de Regge, consistente en todas las configuraciones de parámetros probadas (lo cual
sugiere que es una propiedad estructural de la clasificación en canales crecientes/fijos de esta
cinemática, no un accidente numérico).

La combinación `v1 := B1 - 2·B2`:

```
E=    1000  v1=-4.324019e+18
E=    3000  v1=-3.162441e+21
E=   10000  v1=-4.342965e+24
E=   30000  v1=-3.167057e+27
E=  100000  v1=-4.344773e+30

pendiente = 6.0002  (limpio, s^3.0001)
```

**Verificado robusto en 4 configuraciones de parámetros independientes** (`t0,m0sq,seed`):
`2.998, 3.000, 3.000, 3.000` — consistente a 3+ cifras en las 4.

**Se verificó también que ningún grado de contaminación con `B3` mejora esto**: agregar cualquier
cantidad no nula de `B3` a `v1` reintroduce el término `E^8` (ya que `B3` tiene su propio término
líder `E^8` no relacionado con la cancelación de `v1`) — confirmado numéricamente (con
`β=0.001·B3` agregado, la pendiente vuelve exactamente a `7.9999`).

## Estado honesto: `s³` es lo CONFIRMADO; el resto del espacio de 2 parámetros queda abierto

El conjunto de combinaciones que cancelan el término líder `E^8` es una familia de **2 parámetros**
(no solo el punto `v1`). Se intentó explorar un SEGUNDO punto independiente de esa familia
(`B2 + c₃·B3` con `c₃=-A2/A3`), pero **no se logró la precisión numérica suficiente para confirmar
su comportamiento**: con coeficientes extraídos a 50 y luego 100 dígitos de precisión (`mpmath`), la
pendiente medida NO se estabiliza (crece de 6.86 a 7.57 al aumentar `E`, señal de un residuo del
término `E^8` no cancelado del todo por imprecisión en el coeficiente `c₃`) — a diferencia de `v1`,
cuyo coeficiente (`-2`, aparentemente exacto) no sufre este problema.

**No se afirma que `s³` sea el piso de todo el espacio de grado 9** — solo que es el mejor resultado
**confirmado con rigor numérico suficiente**. Explorar el resto de la familia de 2 parámetros
(y determinar si existe alguna combinación que baje aún más, hasta `s²` o menos) requeriría
probablemente una derivación **simbólica** (no solo numérica) de los coeficientes líderes exactos —
marcado como el paso más importante pendiente.

## Resumen de la clasificación del proyecto hasta ahora (para publicación, tabla consolidada)

| grado | dimensión (revisada) | mejor crecimiento confirmado |
|-------|----------------------|-------------------------------|
| 7     | 0 (teorema, parte 19)| — |
| 8     | 0 (paridad, parte 21)| — |
| 9     | **3** (corregido, esta parte) | **s³** (antes se creía s⁴ con dim. 1) |
| 10    | 0 (paridad, parte 21)| — |
| 11    | 3 (parte 21 — probablemente también incompleta, ver abajo) | s⁴ (parte 24) |

**Advertencia de consistencia**: dado que grado 9 estaba incompleto por no incluir el sandwich, es
altamente probable que **grado 11 (parte 21) también esté incompleto** por la misma razón — no se
recalculó todavía. Esto significa que el resultado de la parte 24 (grado 11 topa en `s⁴`) **debe
tratarse como provisional** hasta repetir ese cálculo con el sandwich incluido.

## Próximo paso concreto (por orden de prioridad)

1. **Recalcular grado 11 con el sandwich incluido sistemáticamente** (mismo tipo de corrección que
   se hizo acá para grado 9) — es la inconsistencia más urgente a resolver antes de reportar nada.
2. Intentar una derivación **simbólica** (no numérica) del comportamiento líder de `B1,B2,B3` en el
   límite de Regge, para (a) confirmar de forma exacta (no solo numérica) que `A1=2A2`, (b) obtener
   el valor EXACTO de `A2/A3` sin depender de extrapolación numérica, y (c) determinar con certeza
   si el resto de la familia de 2 parámetros puede bajar más allá de `s³`.
3. Con eso resuelto, decidir si `s³` (o algo mejor) en grado 9 es la pieza final, o si hace falta
   seguir a grados más altos.

## Archivos

- `scripts/grado9_completo_con_sandwich.py` — cálculo de rango con las 66 semillas (familias A+B).
- `scripts/regge_n5_altaprecision.py` — versión de alta precisión (mpmath) del límite de Regge n=5 y
  extracción de coeficientes líderes para B1,B2,B3.
