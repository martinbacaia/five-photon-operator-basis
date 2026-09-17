"""
Fase 2, paso 16: enumeracion sistematica de candidatos S5-invariantes multilineales de grado
minimo en momentos.

PRIMERO: correccion de un error de conteo de grado cometido en 18-generalizacion-FdotF-y-primer-
candidato.md. Ahi se afirmo "C^k_ab tiene grado 2 en momentos" -- FALSO. Contando potencias de
momento explicitamente en C^k_{ab} = p_a.F^k.p_b = (p_k.p_a)(eps_k.p_b) - (p_k.p_b)(eps_k.p_a):
cada termino tiene 3 momentos (p_k,p_a y luego p_b, o p_k,p_b y p_a) -- GRADO 3, no 2. F^i:F^j
= 2[(p_i.p_j)(eps_i.eps_j) - (p_i.eps_j)(eps_i.p_j)] si tiene grado 2 (dos momentos: p_i,p_j).
Esto cambia la cuenta de grado del candidato ya construido: S = FdotF*FdotF*C*(d_ik-d_jk) tiene
grado 2+2+3+2 = 9 en momentos, NO 7 como se escribio antes.

Con el conteo correcto, la pregunta relevante es: ¿existe un candidato de grado MENOR que 9?
El piso teorico para una estructura "2 parejas F:F + 1 singlete C" (la unica topologia que cubre
las 5 polarizaciones exactamente una vez con los 2 atomos disponibles -- F:F que cubre 2 patas
por objeto de grado 2, y C que cubre 1 pata por objeto de grado 3) es 2+2+3 = 7. El candidato ya
encontrado en la sesion anterior tuvo que agregar un factor extra de grado 2 (d_ik-d_jk) porque la
combinacion "pura" de grado 7 con C^k_{ij} (contrayendo el singlete con AMBOS momentos de la
MISMA pareja que ya aporta el F:F) se cancela identicamente por antisimetria. Este script prueba
la alternativa obvia de grado 7: contraer el singlete con UN momento de CADA pareja (C^k_{il} en
vez de C^k_{ij}), que no tiene el mismo choque de simetria bajo i<->j.

Se prueban tambien variantes de grado 9 adicionales para tener mas de un candidato en la base
(el objetivo final es la BASE completa de grado minimo, no un unico elemento).
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


# -----------------------------------------------------------------------------------------------
# Candidato grado 7, shape (a): C^k_{ij} -- YA SABEMOS que se cancela (chequeo de regresion).
# -----------------------------------------------------------------------------------------------
def T_grado7_shapeA(eta, ps, eps, i, j, l, m, k):
    return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * Cabc(eta, ps, eps, k, i, j)


# -----------------------------------------------------------------------------------------------
# Candidato grado 7, shape (c): C^k_{il} -- singlete contraido con UN momento de CADA pareja.
# -----------------------------------------------------------------------------------------------
def T_grado7_shapeC(eta, ps, eps, i, j, l, m, k):
    return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * Cabc(eta, ps, eps, k, i, l)


# -----------------------------------------------------------------------------------------------
# Candidato ya conocido, grado 9 (parte 18): shape (a) corregida con factor (d_ik - d_jk).
# -----------------------------------------------------------------------------------------------
def T_grado9_v1(eta, ps, eps, i, j, l, m, k):
    d_ik = eta(ps[i], ps[k])
    d_jk = eta(ps[j], ps[k])
    return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * Cabc(eta, ps, eps, k, i, j) * (d_ik - d_jk)


# -----------------------------------------------------------------------------------------------
# Candidato grado 9, variante shape (c) con un factor extra simetrizante -- para ver si ademas
# de shape(c) "pura" (si sobrevive) hay una familia mas rica en grado 9 tambien.
# -----------------------------------------------------------------------------------------------
def T_grado9_v2(eta, ps, eps, i, j, l, m, k):
    d_lk = eta(ps[l], ps[k])
    d_mk = eta(ps[m], ps[k])
    return FdotF(eta, ps, eps, i, j) * FdotF(eta, ps, eps, l, m) * Cabc(eta, ps, eps, k, i, l) * (d_lk - d_mk)


CANDIDATES = {
    "grado7_shapeA (C^k_ij, ya sabido = 0)": T_grado7_shapeA,
    "grado7_shapeC (C^k_il)": T_grado7_shapeC,
    "grado9_v1 (shapeA * (d_ik-d_jk), ya verificado != 0)": T_grado9_v1,
    "grado9_v2 (shapeC * (d_lk-d_mk))": T_grado9_v2,
}


def check_gauge_invariance(name, base_fn, ps, eps, eta, seed_zeta=0):
    S0 = orbit_sum(eta, ps, eps, base_fn)
    rng = np.random.RandomState(seed_zeta)
    zeta = {i: rng.uniform(-3, 3) for i in range(1, 6)}
    eps_shifted = {i: eps[i] + zeta[i] * ps[i] for i in range(1, 6)}
    S1 = orbit_sum(eta, ps, eps_shifted, base_fn)
    return S0, S1


if __name__ == "__main__":
    print("=" * 100)
    print("Enumeracion de candidatos de grado minimo -- gauge invariancia y no trivialidad")
    print("=" * 100)
    for name, fn in CANDIDATES.items():
        print(f"\n--- {name} ---")
        vals = []
        gauge_ok = True
        for seed in [1, 2, 3, 42, 99]:
            ps, eps, eta = gen_kinematics(seed)
            S0, S1 = check_gauge_invariance(name, fn, ps, eps, eta)
            rel_gauge = abs(S0 - S1) / max(1e-12, abs(S0))
            if rel_gauge > 1e-6 and abs(S0) > 1e-6:
                gauge_ok = False
            vals.append(S0)
            print(f"  seed={seed}: S={S0: .6e}   (gauge-shift rel diff = {rel_gauge:.2e})")
        nontrivial = any(abs(v) > 1e-6 for v in vals)
        print(f"  => gauge-invariante: {gauge_ok}   no-trivial (!=0 en alguna semilla): {nontrivial}")
