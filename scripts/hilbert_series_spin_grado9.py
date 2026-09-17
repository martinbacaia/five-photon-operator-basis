"""
Fase 2, paso 17: primer paso hacia la "serie de Hilbert con spin" -- en vez de intentar una formula
cerrada tipo Molien (el propio paper de 2017 admite no saber ni si el anillo con spin es
Cohen-Macaulay para n=5; una formula cerrada esta fuera de alcance de una sesion), se calcula
NUMERICAMENTE la dimension real del espacio de invariantes multilineales (S5-invariantes,
gauge-invariantes, Lorentz-invariantes, grado exactamente 1 en cada una de las 5 polarizaciones)
en grado 9 -- el primer grado no trivial encontrado en 19-enumeracion-grado-minimo.md -- para
responder la pregunta que quedo abierta ahi: ¿el espacio de invariantes de grado 9 tiene dimension
1 o mas?

Metodo (generaliza la tecnica ya usada en las partes 4-6 para el anillo de Mandelstam puro, pero
aplicada ahora a construcciones que incluyen polarizaciones): se arma un conjunto de "semillas"
(productos F:F * F:F * C * (un factor de Mandelstam de grado 2), variando (a) que pareja de
momentos contrae el singlete C y (b) cual de los 10 productos punto d_ab (entre los 5 "slots")
se usa como factor corrector), se simetriza CADA semilla sumando sobre la orbita completa de S5
(exactamente como en 18/19), y se evalua el resultado en VARIAS cinematicas numericas aleatorias
independientes. La dimension real del espacio invariante generado por estas semillas es el RANGO
de la matriz (cinematicas x semillas) -- no se asume nada sobre independencia lineal, se calcula.

Nota de honestidad: esto NO es una serie de Hilbert completa (no cubre todas las topologias
posibles ni se ha probado que grado 9 sea el minimo absoluto, ver 19-enumeracion-grado-minimo.md)
-- es un calculo numerico riguroso de la dimension real dentro de la familia de construcciones ya
identificada, en el grado ya confirmado no-trivial.
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


SLOTS = [1, 2, 3, 4, 5]  # nombres de "posiciones" (i,j,l,m,k) -- se instancian via permutacion


def make_seed(c_contract, mult_pair):
    """
    c_contract: 'P' (contrae con la pareja i,j), 'Q' (con l,m), 'M' (mixta, con i,l)
    mult_pair: tupla de 2 nombres de slot in {i,j,l,m,k} cuyo producto punto multiplica
    """
    slot_names = ["i", "j", "l", "m", "k"]

    def seed_fn(eta, ps, eps, i, j, l, m, k):
        slot_val = {"i": i, "j": j, "l": l, "m": m, "k": k}
        if c_contract == "P":
            c = Cabc(eta, ps, eps, k, i, j)
        elif c_contract == "Q":
            c = Cabc(eta, ps, eps, k, l, m)
        else:
            c = Cabc(eta, ps, eps, k, i, l)
        a, b = slot_val[mult_pair[0]], slot_val[mult_pair[1]]
        mult = eta(ps[a], ps[b])
        return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * c * mult

    return seed_fn


def orbit_sum(eta, ps, eps, base_fn):
    total = 0.0
    for perm in itertools.permutations([1, 2, 3, 4, 5]):
        total += base_fn(eta, ps, eps, *perm)
    return total


if __name__ == "__main__":
    slot_names = ["i", "j", "l", "m", "k"]
    mult_pairs = list(itertools.combinations(slot_names, 2))  # 10 pares
    contractions = ["P", "Q", "M"]

    seeds = {}
    for c in contractions:
        for mp in mult_pairs:
            name = f"C_{c}__d({mp[0]}{mp[1]})"
            seeds[name] = make_seed(c, mp)

    print(f"Total de semillas: {len(seeds)} (3 contracciones x 10 pares multiplicadores)\n")

    seed_names = list(seeds.keys())
    n_samples = 25
    M = np.zeros((n_samples, len(seed_names)))

    for row, kseed in enumerate(range(1, n_samples + 1)):
        ps, eps, eta = gen_kinematics(kseed)
        for col, name in enumerate(seed_names):
            M[row, col] = orbit_sum(eta, ps, eps, seeds[name])

    # normalizar columnas (magnitudes muy distintas entre semillas) antes de calcular rango
    norms = np.linalg.norm(M, axis=0)
    nonzero_cols = norms > 1e-6
    print(f"Semillas identicamente nulas (norma < 1e-6 en {n_samples} muestras): "
          f"{sum(~nonzero_cols)} de {len(seed_names)}")
    for name, nz in zip(seed_names, nonzero_cols):
        if not nz:
            print(f"    (nula) {name}")

    Mn = M[:, nonzero_cols] / norms[nonzero_cols]
    rank = np.linalg.matrix_rank(Mn, tol=1e-6)
    print(f"\nRango numerico de la matriz (muestras x semillas no nulas) = {rank}")
    print(f"(semillas no nulas consideradas: {Mn.shape[1]} de {len(seed_names)} totales, "
          f"{n_samples} muestras cinematicas independientes)")
    print(f"\n=> Dimension real (dentro de esta familia de construccion) del espacio de invariantes "
          f"de grado 9: {rank}")

    sv = np.linalg.svd(Mn, compute_uv=False)
    print("\nValores singulares (para ver el salto rango-a-cero):")
    print("  " + ", ".join(f"{v:.3e}" for v in sv))
