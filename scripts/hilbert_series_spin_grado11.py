"""
Fase 2, paso 18: extension de la tabla numerica de la "serie de Hilbert con spin" a grado 11
(grados 8 y 10 se resuelven aparte, sin computo, por un teorema de paridad -- ver
21-teorema-paridad-y-grado11.md).

Familias probadas en grado 11:
  (A) "2 parejas F:F + 1 singlete C" (base grado 7) x un multiplicador de Mandelstam de GRADO 4
      (monomio de 2 factores d_ab, de los 10 productos punto entre los 5 slots -- 55 monomios:
      10 cuadrados + 45 productos cruzados), para las 3 formas de contraccion del singlete
      (P,Q,mixta) -- 165 semillas.
  (B) "1 pareja F:F + 3 singletes C" (grado 2+3+3+3=11, SIN corrector, topologia genuinamente
      distinta) -- variantes segun con que par de momentos contrae cada uno de los 3 singletes.

Mismo metodo que en grado 9: cada semilla se simetriza sumando sobre la orbita completa de S5, se
evalua en varias cinematicas numericas independientes, y se calcula el RANGO NUMERICO de la matriz
resultante.
"""

import itertools
import numpy as np
from scipy.optimize import brentq

D = 6


def gen_kinematics(seed):
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
    pk, epsk = ps[k], eps[k]
    return eta(pk, ps[a]) * eta(epsk, ps[b]) - eta(pk, ps[b]) * eta(epsk, ps[a])


def orbit_sum(eta, ps, eps, base_fn):
    total = 0.0
    for perm in itertools.permutations([1, 2, 3, 4, 5]):
        total += base_fn(eta, ps, eps, *perm)
    return total


SLOT_NAMES = ["i", "j", "l", "m", "k"]
ALL_PAIRS = list(itertools.combinations(SLOT_NAMES, 2))  # 10


# =================================================================================================
# Familia (A): 2 parejas F:F + 1 singlete C, x multiplicador de Mandelstam grado 4 (producto de 2
# de los 10 pares posibles, con repeticion -- 55 monomios).
# =================================================================================================
def make_seed_A(c_contract, mult_two_pairs):
    """mult_two_pairs: tupla de 2 elementos de ALL_PAIRS (cada uno un par de slots) -- su PRODUCTO
    de productos punto es el multiplicador de grado 4."""
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        slot_val = {"i": i, "j": j, "l": l, "m": m, "k": k}
        if c_contract == "P":
            c = Cabc(eta, ps, eps, k, i, j)
        elif c_contract == "Q":
            c = Cabc(eta, ps, eps, k, l, m)
        else:
            c = Cabc(eta, ps, eps, k, i, l)
        mult = 1.0
        for (x, y) in mult_two_pairs:
            a, b = slot_val[x], slot_val[y]
            mult *= eta(ps[a], ps[b])
        return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * c * mult

    return seed_fn


# =================================================================================================
# Familia (B): 1 pareja F:F + 3 singletes C (grado 11 sin corrector). Pareja = (i,j); singletes
# l,m,k, cada uno contraido con la pareja (i,j) -- variante "natural" mas simple.
# =================================================================================================
def make_seed_B_v1():
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        return (
            FdotF(eta, ps, eps, i, j)
            * Cabc(eta, ps, eps, l, i, j)
            * Cabc(eta, ps, eps, m, i, j)
            * Cabc(eta, ps, eps, k, i, j)
        )
    return seed_fn


def make_seed_B_v2():
    """variante: cada singlete contrae con un par distinto de momentos (l con i,j; m con i,k;
    k con j,l) -- para explorar si aporta algo distinto de v1."""
    def seed_fn(eta, ps, eps, i, j, l, m, k):
        return (
            FdotF(eta, ps, eps, i, j)
            * Cabc(eta, ps, eps, l, i, j)
            * Cabc(eta, ps, eps, m, i, k)
            * Cabc(eta, ps, eps, k, j, l)
        )
    return seed_fn


if __name__ == "__main__":
    seeds = {}

    # familia A: 55 monomios de grado 4 x 3 contracciones = 165 semillas
    mono4 = list(itertools.combinations_with_replacement(ALL_PAIRS, 2))
    assert len(mono4) == 55
    for c in ["P", "Q", "M"]:
        for mp in mono4:
            name = f"A_{c}__d{mp[0]}*d{mp[1]}"
            seeds[name] = make_seed_A(c, mp)

    # familia B: 2 semillas representativas
    seeds["B_v1"] = make_seed_B_v1()
    seeds["B_v2"] = make_seed_B_v2()

    print(f"Total de semillas: {len(seeds)} (165 familia A + 2 familia B)\n")

    seed_names = list(seeds.keys())
    n_samples = 25
    M = np.zeros((n_samples, len(seed_names)))

    for row, kseed in enumerate(range(1, n_samples + 1)):
        ps, eps, eta = gen_kinematics(kseed)
        for col, name in enumerate(seed_names):
            M[row, col] = orbit_sum(eta, ps, eps, seeds[name])

    norms = np.linalg.norm(M, axis=0)
    nonzero_cols = norms > 1e-6
    n_null = int(sum(~nonzero_cols))
    print(f"Semillas identicamente nulas (norma < 1e-6 en {n_samples} muestras): "
          f"{n_null} de {len(seed_names)}")

    Mn = M[:, nonzero_cols] / norms[nonzero_cols]
    rank = np.linalg.matrix_rank(Mn, tol=1e-6)
    print(f"\nRango numerico de la matriz (muestras x semillas no nulas) = {rank}")
    print(f"(semillas no nulas: {Mn.shape[1]} de {len(seed_names)}, {n_samples} muestras)")
    print(f"\n=> Dimension real (dentro de esta familia) del espacio de invariantes de grado 11: {rank}")

    sv = np.linalg.svd(Mn, compute_uv=False)
    print("\nPrimeros 15 valores singulares (para ver el salto rango-a-cero):")
    print("  " + ", ".join(f"{v:.3e}" for v in sv[:15]))
