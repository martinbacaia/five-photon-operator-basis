# Fase 2 — Paso 8 (continuación): fórmula explícita de polarización para fotones, n=5

## Resultado

Generalizando la ecuación (2.11) del paper original (arXiv:1910.14392) — que para n=4 da
ε₁∥ = α₁√(st/u)(p₂/s − p₃/t) + a₁p₁, un único parámetro físico α₁ — se derivó y **verificó** la
fórmula análoga para n=5:

```
ε₁∥ = α₁ · w₂₃ + α₁' · w₂₄ + a₁ · p₁,    donde   w_jk := p_j/s₁ⱼ − p_k/s₁ₖ
```

con **dos** parámetros físicos gauge-invariantes (α₁, α₁') por partícula — confirmando el conteo
dimensional del paso anterior (`11-estructura-polarizacion-fotones-n5.md`). El término `a₁p₁` es
la parte pura de gauge (no física, se cancela en cualquier S-matrix invariante).

## Cómo se verificó (no se aceptó por conteo dimensional solamente)

1. **Identidad de ortogonalidad**: se confirmó algebraicamente que p₁·w_jk=0 para los 6 pares
   posibles (j,k) entre las 4 partículas restantes {2,3,4,5} — misma identidad que en n=4
   (p₁·p_j/s₁ⱼ = −1/2 siempre, por definición de s₁ⱼ=−2p₁·p_j).
2. **Rango real del espacio "crudo"** (antes de eliminar la dirección de gauge p₁): se construyeron
   los 6 vectores w_jk como combinaciones lineales de {p₂,p₃,p₄,p₅} (usando los coeficientes
   1/s₁ⱼ, expresados en términos de las 5 variables de Mandelstam independientes ya validadas en
   pasos anteriores) y se calculó el rango exacto — **da 3**, no 2 como se había estimado
   ingenuamente al principio.
3. **Calibración del método en el caso n=4 ya conocido**: se descubrió el error de estimación
   repitiendo el mismo cálculo para n=4 (donde la respuesta ya se conoce: 1 parámetro físico). El
   rango "crudo" ahí también da 2 (no 1), y el vector que representa a p₁ mismo YA está contenido
   en ese espacio de rango 2 — la dimensión física (1) surge recién al eliminar (cocientar) la
   dirección de p₁. Aplicando la misma lógica a n=5: rango crudo 3, y se confirmó que el vector p₁
   también está contenido en ese span (rango se mantiene en 3 al agregarlo) — la dimensión física
   es 3−1=2, coincidiendo con el conteo geométrico original.
4. **Verificación de que {w₂₃, w₂₄, p₁} es una base completa** del espacio crudo de 3 dimensiones
   (rango exacto 3), y que el vector restante (w₂₅, representante del resto de las combinaciones
   por simetría) se expresa exactamente como combinación lineal de esos 3 — confirmando que la
   fórmula propuesta (usando solo w₂₃ y w₂₄) captura toda la información física, sin perder ni
   sobrar ningún grado de libertad.

## Nota metodológica

La primera estimación ("se espera rango 2") resultó ser un error de razonamiento (confundir el
rango del espacio antes vs. después de cocientar por la dirección de gauge), detectado al obtener
un resultado que no coincidía con lo esperado. En vez de forzar la conclusión, se recalibró el
método contra el caso n=4 ya conocido (donde la respuesta correcta está publicada) antes de confiar
en el resultado de n=5 — mismo principio de verificación aplicado toda la sesión.

## Qué falta todavía

1. La fórmula "inversa" (análoga a la ecuación 2.16 del paper), que exprese (α₁,α₁') en términos de
   contracciones Lorentz-invariantes de field-strength — necesaria para escribir S-matrices
   manifiestamente invariantes de Lorentz y gauge.
2. La acción de S₅ sobre el par (α₁,α₁') por partícula — el análogo de la sección 2.2 del paper
   (que para n=4 mostró una sutileza de fase/signo bajo permutaciones). Con 2 componentes en vez de
   1, es de esperar que la acción de S₅ sea más rica (posiblemente mezclando las 2 componentes
   entre sí bajo ciertas permutaciones, no solo un signo global).
3. Combinar esto con ε₁⊥ (que vive en SO(D−4)) y con el anillo de invariantes de Mandelstam ya
   calculado (parte 7) para construir el módulo completo de S-matrices de fotones.

## Archivos de esta sesión

- `scripts/polarizacion_n5_exploracion.py` — construcción y verificación de la fórmula, incluye la
  calibración contra el caso n=4 conocido.
