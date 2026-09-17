"""
Fase 2, paso 11 (parte A): double-check de las partes 12, 13 y 14(v2) usando cinematica EXPLICITA
(vectores numericos concretos en D=6), en vez de solo algebra simbolica formal sobre productos
punto. Esta es una verificacion mas fuerte porque tambien confirma que resolver la conservacion de
momento simbolicamente (paso usado en TODA la sesion) corresponde a cinematica fisicamente
realizable, no solo a algebra autoconsistente.

Metodo: 4 momentos nulos aleatorios p2,p3,p4 (energia + direccion espacial unitaria), un 5o momento
p5 con direccion aleatoria y energia E5 ajustada (via bisecion) para que p1=-(p2+p3+p4+p5) tambien
sea nulo. Metrica mostly-plus: p.q = -p0 q0 + vec(p).vec(q).
"""

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
    assert E5 is not None, "no se encontro cambio de signo -- reintentar con otra seed"

    def vec(E, n):
        return np.concatenate([[E], E * n])

    p2v, p3v, p4v, p5v = vec(E2, n2), vec(E3, n3), vec(E4, n4), vec(E5, n5)
    p1v = -(p2v + p3v + p4v + p5v)

    def eta(a, b):
        return -a[0] * b[0] + a[1:].dot(b[1:])

    e = rng.randn(D)
    eps1 = e - (eta(e, p1v) / eta(p2v, p1v)) * p2v
    return p1v, p2v, p3v, p4v, p5v, eps1, eta


def aF1b(a, b, eps1_, p1_, eta):
    return eta(a, p1_) * eta(eps1_, b) - eta(a, eps1_) * eta(p1_, b)


def solve_alpha_general(pivot, pA, pB, p1_, eps1_, eta):
    d1piv, d1A, d1B, dpivA, dpivB = (
        eta(p1_, pivot), eta(p1_, pA), eta(p1_, pB), eta(pivot, pA), eta(pivot, pB)
    )
    CpivA = aF1b(pivot, pA, eps1_, p1_, eta)
    CpivB = aF1b(pivot, pB, eps1_, p1_, eta)
    A = -dpivA
    Q = d1piv ** 2 + d1piv * d1A + d1piv * d1B + d1piv * dpivA + d1piv * dpivB + d1A * dpivB + d1B * dpivA
    B = -Q / (2 * d1B)
    Dd = -Q / (2 * d1A)
    E = -dpivB
    M = np.array([[A, B], [Dd, E]])
    return np.linalg.solve(M, [CpivA, CpivB])


if __name__ == "__main__":
    p1v, p2v, p3v, p4v, p5v, eps1, eta = gen_kinematics(42)
    print("chequeo p_i^2=0 y Sum p_i=0:")
    ps = {1: p1v, 2: p2v, 3: p3v, 4: p4v, 5: p5v}
    for i in range(1, 6):
        print(f"  p{i}^2 =", eta(ps[i], ps[i]))
    print("  Sum p_i =", p1v + p2v + p3v + p4v + p5v)

    alpha_old = solve_alpha_general(p2v, p3v, p4v, p1v, eps1, eta)
    print("\n(alpha1,alpha1') en leftover=5:", alpha_old)

    # reconstruccion de eps1_par y comparacion con eps1 real (chequeo de la parte 12/13)
    def s1(pj):
        return -2 * eta(p1v, pj)

    w23 = p2v / s1(p2v) - p3v / s1(p3v)
    w24 = p2v / s1(p2v) - p4v / s1(p4v)
    eps1_par = alpha_old[0] * w23 + alpha_old[1] * w24
    print("\nchequeo (eps1.pm - eps1_par.pm)/d1m debe ser CONSTANTE para m=2,3,4,5:")
    for pm, name in [(p2v, "2"), (p3v, "3"), (p4v, "4"), (p5v, "5")]:
        diff = eta(eps1, pm) - eta(eps1_par, pm)
        print(f"  m={name}: ratio =", diff / eta(p1v, pm))

    # chequeo de la transposicion (2 3) -- parte 14 v2
    alpha_23 = solve_alpha_general(p3v, p2v, p4v, p1v, eps1, eta)
    rho_23 = np.array([[-1, -1], [0, 1]])
    print("\nchequeo (2 3): real=", alpha_23, " predicho=", rho_23.dot(alpha_old))

    # chequeo de la transicion de leftover 5->4 -- parte 14
    alpha_l4 = solve_alpha_general(p2v, p3v, p5v, p1v, eps1, eta)
    s13, s14, s15 = s1(p3v), s1(p4v), s1(p5v)
    T54 = np.array([[1, -s13 / s14], [0, -s15 / s14]])
    print("chequeo leftover 5->4: real=", alpha_l4, " predicho=", T54.dot(alpha_old))
