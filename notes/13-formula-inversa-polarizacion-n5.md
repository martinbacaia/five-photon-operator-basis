# Fase 2 — Paso 9: fórmula inversa (análoga a la ec. 2.16) para n=5

## Resultado

Se derivó y verificó la fórmula que expresa (α₁, α₁') — los 2 parámetros gauge-invariantes de
polarización de la partícula 1 en n=5 (ver `12-formula-explicita-polarizacion-n5.md`) — en
términos de contracciones Lorentz-invariantes del field-strength

```
F¹_μν = p1_μ ε1_ν − p1_ν ε1_μ        (misma definición (2.18) del paper, sin cambios)
```

con pares de momentos de las otras 4 partículas. Definiendo `C_jk := p_j.F¹.p_k` (antisimétrico:
`C_kj = −C_jk`), y usando el par `C₂₃ := p2.F¹.p3`, `C₂₄ := p2.F¹.p4` como datos:

```
C₂₃ = A α₁ + B α₁'                      A = −s₂₃/2 (es decir, A=−d23)
C₂₄ = D α₁ + E α₁'                      E = −s₂₄/2 (es decir, E=−d24)
                                         B = D = −Q/(2 d₁₄)  resp. −Q/(2 d₁₃)
```

con `Q = d₁₂² + d₁₂d₁₃ + d₁₂d₁₄ + d₁₂d₂₃ + d₁₂d₂₄ + d₁₃d₂₄ + d₁₄d₂₃` (dᵢⱼ := pᵢ.pⱼ, y las 5
variables independientes son d₁₂,d₁₃,d₁₄,d₂₃,d₂₄, ya fijadas en pasos anteriores). Invirtiendo el
sistema 2×2:

```
(α₁, α₁') = M⁻¹ (C₂₃, C₂₄),   M = [[A,B],[D,E]],   det M = AE − BD
```

Esta es la generalización genuina de (2.16) — para n=4 el sistema era 1×1 y se resolvía
trivialmente; para n=5 hace falta *un par* de contracciones (no una sola) porque hay 2 parámetros
físicos, y la matriz M no es diagonal (B=D≠0): las dos contracciones "se mezclan".

## Cómo se verificó

1. **Calibración del método contra (2.16) (n=4)**: antes de tocar n=5, se reprodujo *exactamente*
   la fórmula publicada α₁ = 2 p2.F¹.p3/√(stu) repitiendo el mismo procedimiento con el ansatz
   (2.11) del paper para n=4. Se usó una parametrización cuidadosa de las raíces cuadradas
   (`sqrt(s)=rs`, `sqrt(t)=rt`, `sqrt(u)=ru` como generadores independientes, en vez de
   `sqrt(producto)` o `sqrt(cociente)` directos) para evitar falsos negativos de simplificación
   simbólica por ambigüedad de rama — un error real que apareció en el primer intento (dio `-α₁`
   en vez de `α₁` hasta corregir esto). **Con la corrección, coincide exacto.**
2. **Consistencia de gauge**: se verificó que ambas contracciones C₂₃, C₂₄ son independientes de
   `a₁` (el parámetro puro de gauge) — condición necesaria para que sean cantidades físicas.
3. **Round-trip**: sustituyendo la solución (α₁(C), α₁'(C)) de vuelta en las definiciones de
   C₂₃, C₂₄ se recupera la identidad exacta (no aproximada).
4. **Redundancia de una tercera contracción**: se calculó C₂₅ = p2.F¹.p5 y se verificó que es
   combinación lineal *exacta* de C₂₃ y C₂₄ (con coeficientes λ₁=λ₂=−1, es decir
   `C₂₅ = −C₂₃ − C₂₄`) — coherente con que solo hay 2 parámetros físicos independientes, ninguna
   contracción adicional aporta información nueva.

## Archivos

- `scripts/formula_inversa_n5.py` — script ejecutable con los 4 chequeos de arriba, todos pasan.

## Qué falta todavía

Encontrar cómo actúa S₅ sobre (α₁, α₁') — ver `14-accion-s5-alpha-n5.md`, donde se documenta un
hallazgo importante: la elección de "pivote" (usar siempre las partículas 2 y 3, 2 y 4) para fijar
(α₁, α₁') NO es una elección covariante bajo permutaciones, a diferencia de n=4.
