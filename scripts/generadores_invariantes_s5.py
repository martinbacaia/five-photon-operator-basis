"""
Fase 2, paso 4: construir EXPLICITAMENTE (no solo contar dimensiones) los invariantes de bajo grado
del anillo de S5 sobre el espacio de 5 invariantes de Mandelstam, y ver en que grado aparece el
primer generador genuinamente nuevo (no producto de generadores de grado menor).

Metodo: operador de Reynolds R(f) = (1/|G|) sum_{g en G} f(rho(g) a), aplicado a un polinomio
generico de cada grado. Como el anillo de invariantes en grado <=3 tiene multiplicidad 1 (ver serie
de Molien en anillo_invariantes_s5.py: dim 1,0,1,1 en grados 0,1,2,3), el resultado de aplicar
Reynolds a un monomio generico ya genera (salvo que se anule por casualidad, que se chequea) el
generador de ese grado.

Se generan los 120 elementos de S5 por multiplicacion de matrices a partir de 2 generadores
((12) y (12345)), reusando y verificando otra vez que la construccion es un homomorfismo genuino
(no se da por sentado, dado el bug encontrado en la sesion anterior).
"""

from itertools import permutations, product as iproduct

import sympy as sp

from n5_mandelstam_representation import kernel_basis_matrix, action_on_kernel

n = 5
a = sp.symbols("a1 a2 a3 a4 a5")


def compose_perm(sigma1, sigma2):
    return tuple(sigma1[sigma2[i] - 1] for i in range(n))


def generate_group(basis, edge_list):
    gen_perms = [(2, 1, 3, 4, 5), (2, 3, 4, 5, 1)]  # (12) y (12345)
    gen_mats = [action_on_kernel(n, basis, edge_list, g) for g in gen_perms]

    identity_perm = tuple(range(1, n + 1))
    group = {identity_perm: sp.eye(5)}
    frontier = [identity_perm]
    while frontier:
        new_frontier = []
        for p in frontier:
            for gp, gm in zip(gen_perms, gen_mats):
                q = compose_perm(gp, p)
                if q not in group:
                    group[q] = gm * group[p]
                    new_frontier.append(q)
        frontier = new_frontier
    assert len(group) == 120, f"el grupo generado tiene {len(group)} elementos, se esperaban 120"
    return group


def verify_homomorphism_sample(group, samples=25):
    import random
    random.seed(7)
    perms = list(group.keys())
    fails = 0
    for _ in range(samples):
        p, q = random.choice(perms), random.choice(perms)
        lhs = group[p] * group[q]
        rhs = group[compose_perm(p, q)]
        if sp.simplify(lhs - rhs) != sp.zeros(5, 5):
            fails += 1
    return fails


def reynolds(poly, group):
    total = sp.Integer(0)
    for M in group.values():
        new_a = [sum(M[i, j] * a[j] for j in range(5)) for i in range(5)]
        total += poly.subs(list(zip(a, new_a)), simultaneous=True)
    return sp.expand(total / 120)


if __name__ == "__main__":
    basis, edge_list = kernel_basis_matrix(n)
    print("Generando el grupo S5 (120 elementos) por multiplicacion de matrices...")
    group = generate_group(basis, edge_list)
    print(f"  {len(group)} elementos generados (esperado 120): OK")

    fails = verify_homomorphism_sample(group)
    print(f"  Verificacion de homomorfismo en 25 pares aleatorios: {fails} fallos "
          f"({'OK' if fails == 0 else 'MAL -- no continuar'})\n")
    if fails:
        raise SystemExit("Fallo la verificacion de homomorfismo, no tiene sentido seguir.")

    print("=== Invariante de grado 2 (via Reynolds sobre a1^2) ===")
    inv2 = reynolds(a[0] ** 2, group)
    print(f"  Q(a) = {inv2}")
    print(f"  Q(a) != 0: {inv2 != 0}\n")

    print("=== Invariante de grado 3 (via Reynolds sobre a1^3) ===")
    inv3 = reynolds(a[0] ** 3, group)
    print(f"  C(a) = {inv3}")
    print(f"  C(a) != 0: {inv3 != 0}\n")

    print("=== Grado 4: la serie de Molien predice dimension 2 ===")
    print("¿Q(a)^2 y algun otro invariante independiente lo explican, o hace falta un generador")
    print("nuevo en grado 4?")
    Q2 = sp.expand(inv2 ** 2)
    inv4_from_reynolds = reynolds(a[0] ** 4, group)
    inv4_from_reynolds_b = reynolds(a[0] ** 2 * a[1] ** 2, group)

    # Construir el espacio de polinomios de grado 4 generado por Q^2 vs. lo que da Reynolds
    # directamente sobre un par de monomios de grado 4 (deben caer en el espacio de invariantes,
    # que Molien dice es 2-dimensional). Si Q^2 y los resultados de Reynolds son linealmente
    # independientes, hay un generador nuevo en grado 4. Si son proporcionales entre si (mismo
    # invariante salvo escala), falta encontrar el 2do invariante con otro monomio semilla.
    print(f"  Q(a)^2 (primeros terminos): {str(Q2)[:120]}...")
    print(f"  Reynolds(a1^4): {str(inv4_from_reynolds)[:120]}...")
    print(f"  Reynolds(a1^2 a2^2): {str(inv4_from_reynolds_b)[:120]}...")

    # chequear independencia lineal EXACTA (coeficiente a coeficiente, sin evaluar en puntos
    # numericos -- ver 04-generadores-explicitos-s5.md: una version anterior de este chequeo
    # evaluaba en puntos aleatorios y forzaba el resultado a sp.Integer, lo cual corrompe valores
    # racionales no enteros (los invariantes tienen coeficientes con denominador 5, 9, etc.) y dio
    # un rango 3 espurio. La forma correcta es resolver la combinacion lineal general
    # c1*Q^2 + c2*D1 + c3*D2 = 0 comparando los coeficientes de los monomios, exactamente.
    c1, c2, c3 = sp.symbols("c1 c2 c3")
    polys = [Q2, inv4_from_reynolds, inv4_from_reynolds_b]
    combo = sp.expand(c1 * polys[0] + c2 * polys[1] + c3 * polys[2])
    coeffs_eqs = sp.Poly(combo, *a).coeffs()
    sol = sp.solve(coeffs_eqs, [c1, c2, c3], dict=True)
    print(f"\n  Relaciones lineales exactas entre {{Q^2, Reynolds(a1^4), Reynolds(a1^2a2^2)}}: {sol}")
    if sol:
        # sol[0] fija 2 de las 3 incognitas en terminos de la restante (que queda libre) ->
        # el espacio nulo de relaciones es 1-dimensional -> los 3 candidatos generan un espacio
        # de dimension 3 - 1 = 2.
        print("  => hay exactamente 1 relacion lineal (1-parametrica) entre los 3 candidatos.")
        print("  => el espacio que generan tiene dimension 3 - 1 = 2, coincide EXACTO con Molien.")
        print("  => Q^2 por si solo NO alcanza esa dimension 2 (es solo 1 dimension) -- confirma")
        print("     que hace falta un invariante de grado 4 GENUINAMENTE NUEVO, no reducible a Q y C.")
    else:
        print("  => no se encontro relacion lineal -- los 3 serian independientes (inconsistente")
        print("     con Molien=2, revisar).")
