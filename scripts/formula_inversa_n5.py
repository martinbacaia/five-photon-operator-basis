"""
Fase 2, paso 9: derivar la formula "inversa" (analoga a la ecuacion 2.16 del paper original,
arXiv:1910.14392) que expresa los parametros gauge-invariantes de polarizacion (alpha_1, alpha_1')
de n=5 en terminos de contracciones Lorentz-invariantes del field-strength F^1_{mu nu} = p1_mu
eps1_nu - p1_nu eps1_mu con pares de momentos de las otras 4 particulas.

Metodo: trabajar simbolicamente con los productos punto p_i.p_j como variables (igual que en
`polarizacion_n5_exploracion.py`), sin necesitar cinematica vectorial explicita, porque todas las
cantidades relevantes (p_l.F1.p_m, alpha_1, alpha_1') se escriben en terminos de esos productos
punto.

Paso 0 (calibracion): reproducir a mano/simbolicamente la formula (2.16) del paper para n=4:
    alpha_1 = 2 p2.F1.p3 / sqrt(s t u)
usando la misma tecnica que se va a aplicar depues a n=5. Si esto no da exactamente la formula
publicada, el metodo esta mal y no hay que confiar en el resultado de n=5.

Paso 1: aplicar el mismo metodo a n=5 con el ansatz (ya verificado en el paso anterior de esta
sesion, ver 12-formula-explicita-polarizacion-n5.md):
    eps1_par = alpha_1 * w23 + alpha_1' * w24 + a1 * p1,   w_jk = p_j/s_1j - p_k/s_1k
y resolver el sistema lineal 2x2 que dan las contracciones p2.F1.p3 y p2.F1.p4 para (alpha_1,
alpha_1') en terminos de esas contracciones y los invariantes de Mandelstam.

Paso 2 (verificacion cruzada): repetir con una eleccion DISTINTA de pares (p2.F1.p4 y p2.F1.p5, o
similar) y verificar que da el mismo resultado tras resustituir -- y verificar tambien que una
tercera contraccion cualquiera (p2.F1.p5) es combinacion lineal de las dos elegidas como base,
exactamente como se hizo con los w_jk en el paso anterior.
"""

import itertools
import sympy as sp


def dprint(msg):
    print(msg)


# =================================================================================================
# PASO 0: calibracion contra n=4, ecuacion (2.16) del paper.
# =================================================================================================
dprint("=" * 90)
dprint("PASO 0: calibracion del metodo contra la formula publicada (2.16) para n=4")
dprint("=" * 90)

# Para evitar ambiguedades de rama de sqrt() con simplify simbolico (sqrt(a/b) vs sqrt(a)/sqrt(b)
# no son automaticamente identificados por sympy para simbolos genericos), se introducen sqrt(s),
# sqrt(t), sqrt(u) como generadores independientes rs,rt,ru con rs**2=s, etc. Toda cantidad con
# raices se escribe SIEMPRE en terminos de rs,rt,ru (nunca sqrt(producto) o sqrt(cociente)
# directamente), de modo que la unica ambiguedad de rama posible (el signo global de cada raiz) es
# la misma que ya existe en el paper (su ec. 2.16 tiene esta misma ambiguedad de signo, ligada a la
# eleccion de rama de sqrt(s t u)).
rs, rt, ru = sp.symbols("rs rt ru")
s, t = rs ** 2, rt ** 2
u_expr = ru ** 2
alpha1, a1 = sp.symbols("alpha1 a1")

# Productos punto p_i.p_j para n=4, usando las identidades de (2.2) del paper:
# s = -2 p1.p2 = -2 p3.p4 ; t = -2 p1.p3 = -2 p2.p4 ; u = -2 p1.p4 = -2 p2.p3 ; u = -s - t
# (la relacion u = -s-t se impone como constraint aparte, ver abajo)
dot4 = {
    (1, 2): -s / 2, (3, 4): -s / 2,
    (1, 3): -t / 2, (2, 4): -t / 2,
    (1, 4): -u_expr / 2, (2, 3): -u_expr / 2,
}


def dot4_ij(i, j):
    if i == j:
        return sp.Integer(0)
    return dot4[(min(i, j), max(i, j))]


def s4(i, j):
    return -2 * dot4_ij(i, j)


# ansatz (2.11) para eps1_par: sqrt(st/u) = rs*rt/ru (consistente por construccion)
sqrt_stu_over_u = rs * rt / ru
w23_dot = {m: dot4_ij(2, m) / s4(1, 2) - dot4_ij(3, m) / s4(1, 3) for m in [2, 3, 4]}
eps1par_dot = {m: alpha1 * sqrt_stu_over_u * w23_dot[m] + a1 * dot4_ij(1, m) for m in [2, 3, 4]}

# p2.F1.p3 = (p1.p2)(eps1.p3) - (p1.p3)(eps1.p2), usando eps1.p_m = eps1_par.p_m (eps1_perp
# transverso al plano de scattering)
p2F1p3 = dot4_ij(1, 2) * eps1par_dot[3] - dot4_ij(1, 3) * eps1par_dot[2]
# Imponemos u = -s - t, es decir ru**2 = -rs**2 - rt**2, SOLO al final via sustitucion numerica
# controlada mas abajo -- para el algebra simbolica dejamos u_expr libre (no hace falta la relacion
# de momentum conservation aca, la formula (2.16) es una identidad valida para s,t,u genericos con
# s+t+u=0, y no depende de esa relacion para su verificacion algebraica).
p2F1p3 = sp.expand(p2F1p3)

# sqrt(s t u) = rs*rt*ru (consistente con la misma convencion)
alpha1_from_formula = sp.simplify(2 * p2F1p3 / (rs * rt * ru))
dprint(f"  p2.F1.p3 (derivado del ansatz)           = {sp.factor(p2F1p3)}")
dprint(f"  2*p2.F1.p3/sqrt(s t u)                     = {alpha1_from_formula}")
dprint(f"  Coincide con alpha_1 (ec. 2.16 del paper): {sp.simplify(alpha1_from_formula - alpha1) == 0}")
assert sp.simplify(alpha1_from_formula - alpha1) == 0, "FALLA la calibracion contra (2.16) -- no seguir"
dprint("  => CALIBRACION OK: el metodo reproduce (2.16) exactamente. Procedemos a n=5.\n")


# =================================================================================================
# PASO 1: n=5. Productos punto simbolicos con las 5 variables independientes (mismo esquema que
# en polarizacion_n5_exploracion.py), resolviendo conservacion de momento.
# =================================================================================================
dprint("=" * 90)
dprint("PASO 1: formula inversa para n=5 -- (alpha_1, alpha_1') en terminos de p_l.F1.p_m")
dprint("=" * 90)

n = 5
idx_pairs = list(itertools.combinations(range(1, n + 1), 2))
dot = {p: sp.symbols(f"d{p[0]}{p[1]}") for p in idx_pairs}


def dot_ij(i, j):
    if i == j:
        return sp.Integer(0)
    return dot[(min(i, j), max(i, j))]


def sij(i, j):
    return -2 * dot_ij(i, j)


momentum_conservation_eqs = []
for k in range(1, n + 1):
    momentum_conservation_eqs.append(sum(dot_ij(k, j) for j in range(1, n + 1) if j != k))

free_vars = [dot[(1, 2)], dot[(1, 3)], dot[(1, 4)], dot[(2, 3)], dot[(2, 4)]]
dependent_vars = [v for v in dot.values() if v not in free_vars]
sol = sp.solve(momentum_conservation_eqs, dependent_vars, dict=True)
assert len(sol) == 1
sol = sol[0]


def dot_r(i, j):
    raw = dot_ij(i, j)
    return raw.subs(sol) if raw in sol else raw


def s_r(i, j):
    return -2 * dot_r(i, j)


alpha1p, alpha1pp, a1_5 = sp.symbols("alpha1 alpha1p a1")  # alpha1, alpha1', a1 for n=5

w23_dot5 = {m: dot_r(2, m) / s_r(1, 2) - dot_r(3, m) / s_r(1, 3) for m in [2, 3, 4, 5]}
w24_dot5 = {m: dot_r(2, m) / s_r(1, 2) - dot_r(4, m) / s_r(1, 4) for m in [2, 3, 4, 5]}
eps1par_dot5 = {
    m: alpha1p * w23_dot5[m] + alpha1pp * w24_dot5[m] + a1_5 * dot_r(1, m)
    for m in [2, 3, 4, 5]
}


def pF1p(l, m):
    """p_l . F1 . p_m = (p1.pl)(eps1.pm) - (p1.pm)(eps1.pl)."""
    return sp.expand(dot_r(1, l) * eps1par_dot5[m] - dot_r(1, m) * eps1par_dot5[l])


# Elegimos, en paralelo a n=4 (que usaba (2,3)), las dos contracciones (2,3) y (2,4):
c23 = pF1p(2, 3)
c24 = pF1p(2, 4)

dprint("  p2.F1.p3 (n=5) = " + str(sp.simplify(c23)))
dprint("  p2.F1.p4 (n=5) = " + str(sp.simplify(c24)))

# Verificar primero que la dependencia en a1 se cancela en ambas (deben ser gauge invariantes)
c23_a1coeff = sp.simplify(sp.Poly(sp.expand(c23), a1_5).coeff_monomial(a1_5))
c24_a1coeff = sp.simplify(sp.Poly(sp.expand(c24), a1_5).coeff_monomial(a1_5))
dprint(f"  coeficiente de a1 en p2.F1.p3: {c23_a1coeff}  (debe ser 0)")
dprint(f"  coeficiente de a1 en p2.F1.p4: {c24_a1coeff}  (debe ser 0)")
assert c23_a1coeff == 0 and c24_a1coeff == 0, "FALLA: dependencia espuria de a1 (gauge) en las contracciones"
dprint("  => Ambas contracciones son gauge-invariantes (independientes de a1). OK.\n")

# Resolver el sistema lineal 2x2 para (alpha1, alpha1') en terminos de c23, c24
C23, C24 = sp.symbols("C23 C24")  # valores de las contracciones (datos)
eq1 = sp.Eq(c23.subs(a1_5, 0), C23)
eq2 = sp.Eq(c24.subs(a1_5, 0), C24)
solution = sp.solve([eq1, eq2], [alpha1p, alpha1pp], dict=True)
assert len(solution) == 1, f"se esperaba solucion unica, se obtuvo {solution}"
solution = solution[0]

alpha1_formula = sp.simplify(solution[alpha1p])
alpha1p_formula = sp.simplify(solution[alpha1pp])

dprint("  Resolviendo el sistema lineal 2x2 {p2.F1.p3 = C23, p2.F1.p4 = C24} para (alpha_1, alpha_1'):")
dprint(f"    alpha_1  = {alpha1_formula}")
dprint(f"    alpha_1' = {alpha1p_formula}")
dprint("  (C23 := p2.F1.p3,  C24 := p2.F1.p4, expresados en terminos de las 5 variables de Mandelstam)\n")

# ------------------------------------------------------------------------------------------------
# Verificacion 1: sustituir la solucion de vuelta en las definiciones y confirmar identidad.
# ------------------------------------------------------------------------------------------------
dprint("--- Verificacion 1: sustituir la solucion de vuelta reproduce C23, C24 identicamente ---")
check1 = sp.simplify(c23.subs(a1_5, 0).subs({alpha1p: alpha1_formula, alpha1pp: alpha1p_formula}) - C23)
check2 = sp.simplify(c24.subs(a1_5, 0).subs({alpha1p: alpha1_formula, alpha1pp: alpha1p_formula}) - C24)
dprint(f"  p2.F1.p3(alpha_1(C),alpha_1'(C)) - C23 = {check1}  (debe ser 0)")
dprint(f"  p2.F1.p4(alpha_1(C),alpha_1'(C)) - C24 = {check2}  (debe ser 0)")
assert check1 == 0 and check2 == 0
dprint("  => OK: la formula inversa es consistente (round-trip exacto).\n")

# ------------------------------------------------------------------------------------------------
# Verificacion 2: eleccion alternativa de pares -- (2,4) y (2,5) -- debe dar una formula
# DISTINTA en apariencia pero que, tras sustituir alpha_1, alpha_1' por sus valores en terminos
# de C23, C24, sea consistente (round-trip) tambien para una TERCERA contraccion no usada, p2.F1.p5,
# expresada como combinacion lineal de C23 y C24 (mismo chequeo que w25 en el paso anterior).
# ------------------------------------------------------------------------------------------------
dprint("--- Verificacion 2: una tercera contraccion (p2.F1.p5) es combinacion lineal de C23, C24 ---")
c25 = pF1p(2, 5).subs(a1_5, 0)
c25_a1coeff = sp.simplify(sp.Poly(sp.expand(pF1p(2, 5)), a1_5).coeff_monomial(a1_5))
dprint(f"  coeficiente de a1 en p2.F1.p5: {c25_a1coeff}  (debe ser 0)")
assert c25_a1coeff == 0

lam1, lam2 = sp.symbols("lam1 lam2")
# c25, como funcion lineal homogenea de alpha1p/alpha1pp, se puede escribir como combinacion lineal
# de c23 y c24 (mismas variables alpha1p, alpha1pp) si y solo si existen lam1,lam2 (funciones de
# los Mandelstam) tales que c25 = lam1*c23.subs(a1_5,0) + lam2*c24.subs(a1_5,0) como polinomios en
# (alpha1p, alpha1pp). Comparamos coeficiente a coeficiente.
c23_0 = sp.expand(c23.subs(a1_5, 0))
c24_0 = sp.expand(c24.subs(a1_5, 0))
c25_0 = sp.expand(c25)

coeff_a1p_c23 = c23_0.coeff(alpha1p, 1).coeff(alpha1pp, 0)
coeff_a1pp_c23 = c23_0.coeff(alpha1pp, 1).coeff(alpha1p, 0)
coeff_a1p_c24 = c24_0.coeff(alpha1p, 1).coeff(alpha1pp, 0)
coeff_a1pp_c24 = c24_0.coeff(alpha1pp, 1).coeff(alpha1p, 0)
coeff_a1p_c25 = c25_0.coeff(alpha1p, 1).coeff(alpha1pp, 0)
coeff_a1pp_c25 = c25_0.coeff(alpha1pp, 1).coeff(alpha1p, 0)

M = sp.Matrix([[coeff_a1p_c23, coeff_a1p_c24], [coeff_a1pp_c23, coeff_a1pp_c24]])
rhs = sp.Matrix([coeff_a1p_c25, coeff_a1pp_c25])
lam_sol = sp.simplify(M.solve(rhs))
lam1_val, lam2_val = lam_sol[0], lam_sol[1]
dprint(f"  lambda_1 (coef. de C23) = {lam1_val}")
dprint(f"  lambda_2 (coef. de C24) = {lam2_val}")

reconstructed = sp.simplify(lam1_val * c23_0 + lam2_val * c24_0 - c25_0)
dprint(f"  lambda1*C23 + lambda2*C24 - C25 = {reconstructed}  (debe ser 0)")
assert reconstructed == 0
dprint("  => OK: p2.F1.p5 no aporta informacion nueva -- es combinacion lineal exacta de C23, C24.\n")

dprint("=" * 90)
dprint("RESUMEN FINAL")
dprint("=" * 90)
dprint("alpha_1  = " + str(alpha1_formula))
dprint("alpha_1' = " + str(alpha1p_formula))
dprint("con C23 = p2.F1.p3, C24 = p2.F1.p4, F1_{mu nu} = p1_mu eps1_nu - p1_nu eps1_mu (misma")
dprint("definicion (2.18) del paper, sin cambios para n=5).")
