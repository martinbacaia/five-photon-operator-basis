"""
Fase 2, paso 14: generalizar la forma cerrada de F^i:F^j|_par (derivada en 17-forma-cerrada-FdotF.md
para el par (1,2)) a los 10 pares (i,j) posibles con n=5, SIN asumir el patron -- verificando
explicitamente cada uno con cinematica numerica D=6 explicita (mismo metodo que
kinematics_numeric_check.py).

Convencion de pivote/companeras/leftover para el par (i,j) (misma receta usada en las partes
9-17, aplicada de forma sistematica): sea R = sorted({1..5} - {i,j}) = [A,B,leftover] (orden
ascendente). Para la particula i: pivote=j, companeras=(A,B). Para la particula j: pivote=i,
companeras=(A,B) (mismo orden). Para (i,j)=(1,2): R=[3,4,5] -> A=3,B=4,leftover=5 -- coincide
exactamente con la convencion ya usada y verificada en la parte 17.

Formula conjeturada (generalizacion directa de la ya verificada para (1,2)):
    F^i:F^j|par = -(alpha_i alpha_j + alpha_i' alpha_j')
                  - Q_ij * [ alpha_i alpha_j' / (2 d_iA d_jB) + alpha_i' alpha_j / (2 d_iB d_jA) ]
    Q_ij = d_ij^2 + d_ij d_iA + d_ij d_iB + d_ij d_jA + d_ij d_jB + d_iA d_jB + d_iB d_jA

Metodo de verificacion (igual al de la parte 17, "verificacion numerica de la parte paralela
sola"): se generan eps_i ALEATORIOS (con eps_i.p_i=0, componente perp incluida), se extraen
(alpha_i,alpha_i') via la formula inversa general (gauge-invariante por construccion, no le
importa la parte perp/gauge de eps_i), se reconstruye eps_i^par_puro SOLO con los w's (sin
termino a_i*p_i), y se compara F^i:F^j calculado DIRECTO por vectores sobre esas eps^par_puras
contra la formula cerrada conjeturada.
"""

import itertools
import numpy as np
from scipy.optimize import brentq

D = 6
N = 5


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
    assert E5 is not None, "no se encontro cambio de signo -- reintentar con otra seed"

    def vec(E, n):
        return np.concatenate([[E], E * n])

    p2v, p3v, p4v, p5v = vec(E2, n2), vec(E3, n3), vec(E4, n4), vec(E5, n5)
    p1v = -(p2v + p3v + p4v + p5v)

    def eta(a, b):
        return -a[0] * b[0] + a[1:].dot(b[1:])

    ps = {1: p1v, 2: p2v, 3: p3v, 4: p4v, 5: p5v}
    for k in range(1, 6):
        assert abs(eta(ps[k], ps[k])) < 1e-8, f"p{k}^2 != 0"
    assert np.linalg.norm(sum(ps.values())) < 1e-8

    # eps_i aleatorio con eps_i.p_i = 0, para cada particula (incluye componente perp real)
    eps = {}
    for i in range(1, 6):
        e = rng.randn(D)
        # proyectar contra ALGUN otro momento para fijar eps_i.p_i=0 sin matar la parte perp:
        other = 2 if i != 2 else 3
        e = e - (eta(e, ps[i]) / eta(ps[other], ps[i])) * ps[other]
        assert abs(eta(e, ps[i])) < 1e-8
        eps[i] = e

    return ps, eps, eta


def s_(eta, pi, pm):
    return -2 * eta(pi, pm)


def pFq(eta, pi, epsi, pl, pm):
    """p_l . F^i . p_m = (p_i.p_l)(eps_i.p_m) - (p_i.p_m)(eps_i.p_l)."""
    return eta(pi, pl) * eta(epsi, pm) - eta(pi, pm) * eta(epsi, pl)


def solve_alpha(eta, pi, epsi, pivot, pA, pB):
    """(alpha_i, alpha_i') para particula i con pivote=pivot, companeras=(A,B)."""
    d_i_piv = eta(pi, pivot)
    d_i_A = eta(pi, pA)
    d_i_B = eta(pi, pB)
    d_piv_A = eta(pivot, pA)
    d_piv_B = eta(pivot, pB)
    C_pivA = pFq(eta, pi, epsi, pivot, pA)
    C_pivB = pFq(eta, pi, epsi, pivot, pB)
    Q = (d_i_piv ** 2 + d_i_piv * d_i_A + d_i_piv * d_i_B + d_i_piv * d_piv_A
         + d_i_piv * d_piv_B + d_i_A * d_piv_B + d_i_B * d_piv_A)
    M = np.array([
        [-d_piv_A, -Q / (2 * d_i_B)],
        [-Q / (2 * d_i_A), -d_piv_B],
    ])
    alpha, alphap = np.linalg.solve(M, [C_pivA, C_pivB])
    return alpha, alphap, Q


def w_vec(eta, pi, ps, piv, comp):
    """w_{piv,comp} para particula i: p_piv/s_i(piv) - p_comp/s_i(comp)."""
    return ps[piv] / s_(eta, pi, ps[piv]) - ps[comp] / s_(eta, pi, ps[comp])


def FdotF_direct(eta, pi, epsi, pj, epsj):
    return 2 * (eta(pi, pj) * eta(epsi, epsj) - eta(pi, epsj) * eta(epsi, pj))


def test_pair(ps, eps, eta, i, j, verbose=True):
    remaining = sorted(set(range(1, 6)) - {i, j})
    A, B, leftover = remaining[0], remaining[1], remaining[2]

    pi, pj = ps[i], ps[j]
    epsi, epsj = eps[i], eps[j]

    alpha_i, alphap_i, Q_i = solve_alpha(eta, pi, epsi, ps[j], ps[A], ps[B])
    alpha_j, alphap_j, Q_j = solve_alpha(eta, pj, epsj, ps[i], ps[A], ps[B])

    # Q debe ser simetrica en i<->j (misma cantidad para ambas particulas, ya que solo depende de
    # {i,j,A,B} como conjunto etiquetado -- chequeo explicito, no asumido.
    assert abs(Q_i - Q_j) < 1e-6 * max(1, abs(Q_i)), f"Q_i != Q_j para par ({i},{j}): {Q_i} vs {Q_j}"
    Q = Q_i

    # reconstruccion pura-paralela (sin termino a_i*p_i)
    eps_i_par = alpha_i * w_vec(eta, pi, ps, j, A) + alphap_i * w_vec(eta, pi, ps, j, B)
    eps_j_par = alpha_j * w_vec(eta, pj, ps, i, A) + alphap_j * w_vec(eta, pj, ps, i, B)

    FdotF_pure = FdotF_direct(eta, pi, eps_i_par, pj, eps_j_par)

    d_ij = eta(pi, pj)
    d_iA, d_iB = eta(pi, ps[A]), eta(pi, ps[B])
    d_jA, d_jB = eta(pj, ps[A]), eta(pj, ps[B])

    formula = (
        -(alpha_i * alpha_j + alphap_i * alphap_j)
        - Q * (alpha_i * alphap_j / (2 * d_iA * d_jB) + alphap_i * alpha_j / (2 * d_iB * d_jA))
    )

    diff = FdotF_pure - formula
    rel = abs(diff) / max(1e-12, abs(FdotF_pure))
    if verbose:
        print(f"  par ({i},{j})  A={A},B={B},leftover={leftover}:"
              f"  directo={FdotF_pure: .6e}  formula={formula: .6e}  diff={diff: .2e}  rel={rel:.2e}")
    return rel


if __name__ == "__main__":
    pairs = list(itertools.combinations(range(1, 6), 2))
    print(f"Probando los {len(pairs)} pares (i,j) con la formula cerrada generalizada, 3 seeds distintas:\n")
    worst = 0.0
    for seed in [1, 7, 123]:
        ps, eps, eta = gen_kinematics(seed)
        print(f"--- seed={seed} ---")
        for (i, j) in pairs:
            rel = test_pair(ps, eps, eta, i, j)
            worst = max(worst, rel)
        print()
    print(f"Peor error relativo sobre todos los pares/seeds: {worst:.3e}")
    assert worst < 1e-8, "FALLA: la formula generalizada NO reproduce F^i:F^j para todos los pares"
    print("=> OK: la formula cerrada generalizada de F^i:F^j|_par vale para los 10 pares, "
          "con la convencion de companeras = 2 indices mas chicos entre los 3 restantes, "
          "leftover = el mas grande.")
