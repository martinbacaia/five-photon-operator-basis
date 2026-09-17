"""
Fase 2, paso 10 (v2 - CORRECCION de un error de metodo en s5_accion_alpha_n5.py).

Bug encontrado en el primer intento: al calcular rho(sigma), se resolvia el sistema
(C_{sigma(2)sigma(3)}, C_{sigma(2)sigma(4)}) = M * (alpha1_nuevo, alpha1p_nuevo) usando la matriz
M ORIGINAL (construida con el pivote fijo en la particula "2" de la config vieja). Pero la
definicion correcta de "alpha1_nuevo" (el parametro fisico de la particula 1 en la config
RELABELED por sigma) usa como pivote la particula que en la nueva config se llama "2" -- que es
la particula sigma(2) de la config vieja. Por lo tanto hay que usar la matriz M evaluada con el
pivote sigma(2) y companeras sigma(3), sigma(4) (M_sigma), NO la matriz M con pivote 2 fijo.

Este script corrige eso: construye M_for_pivot(pivot,partner1,partner2) de forma generica y
recalcula rho(sigma) = M_sigma^{-1} N(sigma), donde N(sigma) son los coeficientes de
C_{sigma(2)sigma(3)}, C_{sigma(2)sigma(4)} en la base (alpha1,alpha1') VIEJA (pivote 2).

Chequeo clave: es ahora rho(sigma) una matriz CONSTANTE (independiente de la cinematica)?
"""

import itertools
import random

import sympy as sp

n = 5
idx_pairs = list(itertools.combinations(range(1, n + 1), 2))
dot = {p: sp.symbols(f"d{p[0]}{p[1]}") for p in idx_pairs}


def dot_ij(i, j):
    if i == j:
        return sp.Integer(0)
    return dot[(min(i, j), max(i, j))]


eqs = [sum(dot_ij(k, j) for j in range(1, n + 1) if j != k) for k in range(1, n + 1)]
free_vars = [dot[(1, 2)], dot[(1, 3)], dot[(1, 4)], dot[(2, 3)], dot[(2, 4)]]
dependent_vars = [v for v in dot.values() if v not in free_vars]
sol = sp.solve(eqs, dependent_vars, dict=True)
assert len(sol) == 1
sol = sol[0]


def dr(i, j):
    raw = dot_ij(i, j)
    return raw.subs(sol) if raw in sol else raw


def sr(i, j):
    return -2 * dr(i, j)


alpha1, alpha1p = sp.symbols("alpha1 alpha1p")


def w(pivot, partner, m):
    """w_{pivot,partner} . p_m = p_pivot.p_m/s_{1,pivot} - p_partner.p_m/s_{1,partner}."""
    return dr(pivot, m) / sr(1, pivot) - dr(partner, m) / sr(1, partner)


def eps1par_dot(pivot, partner1, partner2, m, a, ap):
    return a * w(pivot, partner1, m) + ap * w(pivot, partner2, m)


def pF1p_generic(l, m, pivot, partner1, partner2, a, ap):
    """p_l.F1.p_m con eps1_par escrita en la base (pivot,partner1,partner2)."""
    return sp.expand(
        dr(1, l) * eps1par_dot(pivot, partner1, partner2, m, a, ap)
        - dr(1, m) * eps1par_dot(pivot, partner1, partner2, l, a, ap)
    )


def coeffs_of(expr, a, ap):
    return expr.coeff(a, 1).coeff(ap, 0), expr.coeff(ap, 1).coeff(a, 0)


def M_for_pivot(pivot, partner1, partner2):
    """Matriz 2x2 tal que (C_{pivot,partner1}, C_{pivot,partner2}) = M (alpha, alpha')
    en la base definida por ESE pivote/companeras."""
    c1 = pF1p_generic(pivot, partner1, pivot, partner1, partner2, alpha1, alpha1p)
    c2 = pF1p_generic(pivot, partner2, pivot, partner1, partner2, alpha1, alpha1p)
    a1c, b1c = coeffs_of(c1, alpha1, alpha1p)
    a2c, b2c = coeffs_of(c2, alpha1, alpha1p)
    return sp.Matrix([[a1c, b1c], [a2c, b2c]])


# Matriz vieja: pivote 2, companeras 3,4
M_old = M_for_pivot(2, 3, 4)
M_old_inv = M_old.inv()
print("M (pivote 2, companeras 3,4):")
sp.pprint(sp.simplify(M_old))
print()


def rho_v2(sigma):
    """sigma: dict {2,3,4,5}->{2,3,4,5}. Devuelve rho(sigma) usando la matriz M_sigma CORRECTA
    (pivote sigma(2), companeras sigma(3), sigma(4)), no la matriz vieja."""
    j2, j3, j4 = sigma[2], sigma[3], sigma[4]
    # N(sigma): C_{j2,j3} y C_{j2,j4} en terminos de alpha1,alpha1' VIEJOS (base pivote=2,3,4)
    new_C23 = pF1p_generic(j2, j3, 2, 3, 4, alpha1, alpha1p)
    new_C24 = pF1p_generic(j2, j4, 2, 3, 4, alpha1, alpha1p)
    a23, b23 = coeffs_of(new_C23, alpha1, alpha1p)
    a24, b24 = coeffs_of(new_C24, alpha1, alpha1p)
    N = sp.Matrix([[a23, b23], [a24, b24]])

    M_sigma = M_for_pivot(j2, j3, j4)
    return sp.simplify(M_sigma.inv() * N)


elements_234 = [2, 3, 4, 5]
sample_perms = [(3, 2, 4, 5), (2, 4, 3, 5), (2, 3, 5, 4), (4, 3, 2, 5), (5, 3, 4, 2), (3, 4, 5, 2)]

print("=== Chequeo: es rho(sigma) (v2, con M_sigma correcta) independiente de la cinematica? ===")
random.seed(1)
free_syms = [dot[(1, 2)], dot[(1, 3)], dot[(1, 4)], dot[(2, 3)], dot[(2, 4)]]
for perm in sample_perms:
    sigma = dict(zip(elements_234, perm))
    R = rho_v2(sigma)
    free = R.free_symbols
    is_const = len(free) == 0
    print(f"  sigma: 2->{sigma[2]},3->{sigma[3]},4->{sigma[4]},5->{sigma[5]}   simbolico const: {is_const}")
    if not is_const:
        # verificacion numerica en 2 puntos, para no confiar solo en "tiene simbolos"
        vals1 = {v: sp.Rational(random.randint(-9, 9), random.choice([1, 2, 3])) for v in free_syms}
        vals2 = {v: sp.Rational(random.randint(-9, 9), random.choice([1, 2, 3])) for v in free_syms}
        R1 = R.subs(vals1)
        R2 = R.subs(vals2)
        print(f"    R en punto 1: {[sp.nsimplify(x) for x in R1]}")
        print(f"    R en punto 2: {[sp.nsimplify(x) for x in R2]}")
    else:
        print(f"    R = {R.tolist()}")
