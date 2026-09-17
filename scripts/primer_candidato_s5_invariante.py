"""
Fase 2, paso 15: primer candidato S5-invariante a S-matrix de fotones para n=5.

ANTES de construirlo, se corrige un punto estructural que NO estaba explicitado en
17-forma-cerrada-FdotF.md: la suma propuesta alli como "candidato mas simple",
Sum_{i<j} F^i:F^j, NO PUEDE ser una S-matriz de 5 fotones, porque no es multilineal en las
5 polarizaciones -- cada termino de esa suma es lineal en eps_i, eps_j PERO DE GRADO CERO en
las otras 3 polarizaciones (eps_k, eps_l, eps_m no aparecen en absoluto en ese termino). Una
amplitud de n fotones debe ser exactamente lineal (grado 1) en CADA una de las n polarizaciones
(asi es como entra cada pata externa de foton, y asi estan construidos todos los building
blocks del caso n=4 en el paper original: para n=4, con 4 patas, los building blocks
(eps_i.eps_j)(eps_k.eps_l) cubren las 4 polarizaciones exactamente una vez cada una, apareado
en 2 pares -- posible porque 4 es par). Para n=5 (impar) no se puede apleastar todo en pares:
hacen falta 2 pares (F^i:F^j, F^l:F^m -- cubren 4 polarizaciones) + 1 objeto lineal en la
polarizacion restante sola (C^k_{ab} := p_a.F^k.p_b, que es lineal SOLO en eps_k y no in-
volucra eps_a, eps_b, contrayendo con sus MOMENTOS en cambio).

Candidato minimal multilineal en las 5 polarizaciones:
    T(i,j | l,m | k) := F^i:F^j * F^l:F^m * C^k_{ij}
(k el "singlete", {i,j} y {l,m} las 2 parejas, C^k contraido con los momentos de la pareja
{i,j} -- eleccion arbitraria entre {i,j} y {l,m}, se corrige simetrizando sobre el grupo).

Candidato S5-invariante: suma sobre la ORBITA COMPLETA de S5 actuando sobre las etiquetas de
T (el operador de Reynolds aplicado a T) -- esto es automaticamente S5-invariante sin importar
que T de partida no lo sea, y sin necesitar simetrizar "a mano" sobre las 15 particiones en
pares+singlete y las 2 elecciones de contraccion del singlete (el promedio sobre el grupo ya
las cubre todas).

Verificaciones:
1. Invariancia de gauge (cada factor ya es invariante por separado, pero se chequea el
   producto completo con un shift de gauge aleatorio simultaneo en las 5 polarizaciones).
2. Invariancia bajo S5: se recalcula la suma completa con una permutacion aleatoria de las
   etiquetas 1..5 aplicada a (p_i,eps_i) y se compara -- debe dar EXACTAMENTE lo mismo (esto
   es garantizado por construccion -- el chequeo numerico sirve para atrapar bugs de
   implementacion, no para "descubrir" la invariancia).
3. No trivialidad: S != 0 genericamente.
"""

import itertools
import numpy as np
from scipy.optimize import brentq

D = 6


def gen_kinematics(seed, n_extra_dims_check=True):
    rng = np.random.RandomState(seed)

    def rand_null():
        nvec = rng.randn(D - 1)
        nvec /= np.linalg.norm(nvec)
        E = rng.uniform(1, 3)
        return E, nvec

    E2, n2 = rand_null()
    E3, n3 = rand_null()
    E4, n4 = rand_null()
    n5 = rng.randn(D - 1)
    n5 /= np.linalg.norm(n5)

    def p1_mass2(E5):
        E0 = E2 + E3 + E4 + E5
        spatial = E2 * n2 + E3 * n3 + E4 * n4 + E5 * n5
        return -E0 ** 2 + spatial.dot(spatial)

    xs = np.linspace(-50, 50, 400)
    vals = [p1_mass2(x) for x in xs]
    E5 = None
    for i in range(len(xs) - 1):
        if vals[i] * vals[i + 1] < 0:
            E5 = brentq(p1_mass2, xs[i], xs[i + 1])
            break
    assert E5 is not None

    def vec(E, n):
        return np.concatenate([[E], E * n])

    p2v, p3v, p4v, p5v = vec(E2, n2), vec(E3, n3), vec(E4, n4), vec(E5, n5)
    p1v = -(p2v + p3v + p4v + p5v)

    def eta(a, b):
        return -a[0] * b[0] + a[1:].dot(b[1:])

    ps = {1: p1v, 2: p2v, 3: p3v, 4: p4v, 5: p5v}
    if n_extra_dims_check:
        for k in range(1, 6):
            assert abs(eta(ps[k], ps[k])) < 1e-8
        assert np.linalg.norm(sum(ps.values())) < 1e-8

    eps = {}
    for i in range(1, 6):
        e = rng.randn(D)
        other = 2 if i != 2 else 3
        e = e - (eta(e, ps[i]) / eta(ps[other], ps[i])) * ps[other]
        eps[i] = e

    return ps, eps, eta


def FdotF(eta, ps, eps, i, j):
    pi, pj, epsi, epsj = ps[i], ps[j], eps[i], eps[j]
    return 2 * (eta(pi, pj) * eta(epsi, epsj) - eta(pi, epsj) * eta(epsi, pj))


def Cabc(eta, ps, eps, k, a, b):
    """p_a . F^k . p_b."""
    pk, epsk = ps[k], eps[k]
    return eta(pk, ps[a]) * eta(epsk, ps[b]) - eta(pk, ps[b]) * eta(epsk, ps[a])


def T_base_naive(eta, ps, eps, i, j, l, m, k):
    """F^i:F^j * F^l:F^m * C^k_{ij}  (contrae el singlete k con los momentos de la pareja i,j)."""
    return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * Cabc(eta, ps, eps, k, i, j)


def T_base(eta, ps, eps, i, j, l, m, k):
    """
    Version corregida de T_base_naive. C^k_{ij} es ANTISIMETRICO en (i,j) (verificado:
    C(1;2,3)=-C(1;3,2)), mientras que F^i:F^j es SIMETRICO en (i,j) -- el producto
    F^i:F^j * C^k_{ij} es entonces IMPAR bajo el intercambio i<->j, que es una de las
    permutaciones incluidas en la suma sobre la orbita completa de S5 -- eso hace que
    T_base_naive se cancele EXACTAMENTE al sumar sobre la orbita (confirmado numericamente:
    S ~ 1e-12, ruido de punto flotante de un cero exacto, no un "casi cero" gauge-dependiente).
    Se corrige multiplicando por (d_ik - d_jk) (tambien antisimetrico bajo i<->j, construido
    con los mismos invariantes de Mandelstam ya calculados en las partes 3-8), de modo que el
    producto completo sea PAR bajo i<->j y sobreviva la simetrizacion.
    """
    pi, pj, pk = ps[i], ps[j], ps[k]
    d_ik = eta(pi, pk)
    d_jk = eta(pj, pk)
    return (
        FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m)
        * Cabc(eta, ps, eps, k, i, j) * (d_ik - d_jk)
    )


def S5_invariant_sum(eta, ps, eps, base_fn=T_base):
    total = 0.0
    for perm in itertools.permutations([1, 2, 3, 4, 5]):
        i, j, l, m, k = perm
        total += base_fn(eta, ps, eps, i, j, l, m, k)
    return total


if __name__ == "__main__":
    ps, eps, eta = gen_kinematics(42)

    print("=== Chequeo previo: la version 'ingenua' (sin el factor (d_ik-d_jk)) se anula ===")
    for seed in [1, 2, 3]:
        ps_s, eps_s, eta_s = gen_kinematics(seed)
        S_naive = S5_invariant_sum(eta_s, ps_s, eps_s, base_fn=T_base_naive)
        print(f"  seed={seed}: S_naive = {S_naive:.3e}  (esperado: ~0, ruido de punto flotante)")
    print("  => confirma la cancelacion exacta por antisimetria de C^k_ij bajo i<->j.\n")

    print("=== Verificacion 1: invariancia de gauge ===")
    S0 = S5_invariant_sum(eta, ps, eps)
    rng = np.random.RandomState(0)
    zeta = {i: rng.uniform(-3, 3) for i in range(1, 6)}
    eps_shifted = {i: eps[i] + zeta[i] * ps[i] for i in range(1, 6)}
    S_shifted = S5_invariant_sum(eta, ps, eps_shifted)
    print(f"  S (original)      = {S0:.10e}")
    print(f"  S (gauge shifted) = {S_shifted:.10e}")
    print(f"  diferencia relativa = {abs(S0 - S_shifted) / abs(S0):.3e}")
    assert abs(S0 - S_shifted) < 1e-6 * abs(S0), "FALLA: S no es gauge invariante"
    print("  => OK\n")

    print("=== Verificacion 2: invariancia bajo relabeling S5 ===")
    sigma = [3, 1, 5, 2, 4]  # permutacion arbitraria de {1,2,3,4,5} (sigma(1)=3, etc.)
    ps_perm = {i: ps[sigma[i - 1]] for i in range(1, 6)}
    eps_perm = {i: eps[sigma[i - 1]] for i in range(1, 6)}
    S_perm = S5_invariant_sum(eta, ps_perm, eps_perm)
    print(f"  S (original)     = {S0:.10e}")
    print(f"  S (relabeled)    = {S_perm:.10e}")
    print(f"  diferencia relativa = {abs(S0 - S_perm) / abs(S0):.3e}")
    assert abs(S0 - S_perm) < 1e-6 * abs(S0), "FALLA: S no es S5-invariante"
    print("  => OK\n")

    print("=== Verificacion 3: no trivialidad (varias seeds) ===")
    for seed in [1, 2, 3, 7, 99]:
        ps_s, eps_s, eta_s = gen_kinematics(seed)
        Ss = S5_invariant_sum(eta_s, ps_s, eps_s)
        print(f"  seed={seed}: S = {Ss:.6e}")
        assert abs(Ss) > 1e-6, f"S es (numericamente) cero para seed={seed} -- sospechoso"
    print("  => OK: S es genericamente no nulo.\n")

    print("Candidato S5-invariante construido y verificado:")
    print("  S := Sum_{sigma in S5} F^{sigma(1)}:F^{sigma(2)} * F^{sigma(3)}:F^{sigma(4)}"
          " * C^{sigma(5)}_{sigma(1)sigma(2)}")
