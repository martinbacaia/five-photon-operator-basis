"""
Fase 2, paso 5: revisar grados 5 y 6 del anillo de invariantes de S5 -- ver si los PRODUCTOS de
los generadores ya encontrados (Q grado2, C grado3, E4 grado4) alcanzan las dimensiones que predice
Molien (2 en grado 5, 5 en grado 6), o si hace falta un generador nuevo (o aparece la primera
relacion).

Reutiliza el grupo y los invariantes Q, C, E4 ya re-verificados en reverificacion_generadores.py.
"""

import sympy as sp

from n5_mandelstam_representation import kernel_basis_matrix, action_on_kernel

n = 5
a = sp.symbols("a1 a2 a3 a4 a5")


def compose_perm(sigma1, sigma2):
    return tuple(sigma1[sigma2[i] - 1] for i in range(n))


def generate_group():
    basis, edge_list = kernel_basis_matrix(n)
    gen_perms = [(2, 1, 3, 4, 5), (2, 3, 4, 5, 1)]
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
    assert len(group) == 120
    return group


def reynolds(poly, group):
    total = sp.Integer(0)
    for M in group.values():
        new_a = [sum(M[i, j] * a[j] for j in range(5)) for i in range(5)]
        total += poly.subs(list(zip(a, new_a)), simultaneous=True)
    return sp.expand(total / len(group))


def rank_of_span(polys):
    """
    Rango EXACTO (no numerico) del espacio generado por una lista de polinomios homogeneos del
    mismo grado, via reduccion por filas sobre sus vectores de coeficientes (monomios como base).

    NOTA DE UN BUG REAL ENCONTRADO Y CORREGIDO: una semilla mal escrita a mano (grado 5 en vez de
    6) mezclada con genuinos invariantes de grado 6 daba un rango espurio +1 (vive en un subespacio
    de monomios totalmente disjunto, luego "parece" independiente sin serlo de verdad frente a la
    pregunta que interesa). Por eso esta funcion ahora EXIGE que todos los polinomios pasados sean
    homogeneos del mismo grado, y frena con un error explicito si no lo son -- en vez de devolver
    silenciosamente un numero que mezcla grados distintos.
    """
    all_monoms = set()
    poly_dicts = []
    degrees_seen = set()
    for p in polys:
        pd = sp.Poly(p, *a).as_dict()
        degs = {sum(k) for k in pd.keys()}
        degrees_seen |= degs
        poly_dicts.append(pd)
        all_monoms |= set(pd.keys())
    assert len(degrees_seen) <= 1, (
        f"rank_of_span recibio polinomios de grados MEZCLADOS {degrees_seen} -- "
        f"esto es casi siempre un error en la semilla (ver nota en el docstring)."
    )
    all_monoms = sorted(all_monoms)
    rows = []
    for pd in poly_dicts:
        rows.append([pd.get(m, 0) for m in all_monoms])
    M = sp.Matrix(rows)
    return M.rank(), M


if __name__ == "__main__":
    print("Generando el grupo S5 (120 elementos)...")
    group = generate_group()
    print("  listo.\n")

    Q = reynolds(a[0] ** 2, group)
    C = reynolds(a[0] ** 3, group)
    D2 = reynolds(a[0] ** 2 * a[1] ** 2, group)  # usado como representante del generador de grado 4
    Q2 = sp.expand(Q ** 2)

    # Aislar E4: un invariante de grado 4 independiente de Q^2 (ya demostrado). Se usa D2
    # directamente (no hace falta "ortogonalizar" para que sirva como generador).
    E4 = D2

    print("=" * 90)
    print("GRADO 5 (Molien predice dimension 2)")
    print("=" * 90)
    QC = sp.expand(Q * C)
    seed_candidates_5 = [
        reynolds(a[0] ** 5, group),
        reynolds(a[0] ** 3 * a[1] ** 2, group),
        reynolds(a[0] * a[1] * a[2] * a[3] * a[4], group),
    ]
    all_5 = [QC] + seed_candidates_5
    rank5, _ = rank_of_span(all_5)
    print(f"Rango del espacio generado por {{Q*C}} + 3 semillas de Reynolds nuevas: {rank5}")
    print(f"(Molien dice que la dimension real de grado 5 es 2)")
    if rank5 == 1:
        print("=> Q*C por si solo (dimension 1) NO alcanza -- pero las otras semillas tampoco dan")
        print("   nada nuevo (coinciden con Q*C) -- revisar con mas semillas.")
    elif rank5 >= 2:
        print("=> Se confirma un invariante de grado 5 ADEMAS de Q*C: si es proporcional a Q*C,")
        print("   productos de generadores ya conocidos alcanzan la dimension 2 (nada nuevo).")
        print("   Si es independiente de Q*C, hace falta un generador nuevo en grado 5.")
        # chequear especificamente si algo es independiente de QC:
        rank_QC_alone, _ = rank_of_span([QC])
        for i, cand in enumerate(seed_candidates_5):
            combo_rank, _ = rank_of_span([QC, cand])
            print(f"   semilla {i}: junto con Q*C, rango = {combo_rank} "
                  f"({'independiente de Q*C' if combo_rank == 2 else 'proporcional a Q*C'})")

    print()
    print("=" * 90)
    print("GRADO 6 (Molien predice dimension 5)")
    print("=" * 90)
    Q3 = sp.expand(Q ** 3)
    C2 = sp.expand(C ** 2)
    E4Q = sp.expand(E4 * Q)
    known_products_6 = [Q3, C2, E4Q]
    rank6_known, _ = rank_of_span(known_products_6)
    print(f"Rango de {{Q^3, C^2, E4*Q}} (productos de generadores ya conocidos): {rank6_known}")
    print(f"(Molien dice que la dimension real de grado 6 es 5)")

    # buscar invariantes adicionales de grado 6 via Reynolds con varias semillas
    extra_seeds_6 = [
        a[0] ** 6, a[0] ** 4 * a[1] ** 2, a[0] ** 3 * a[1] ** 3,
        a[0] ** 2 * a[1] ** 2 * a[2] ** 2, a[0] ** 2 * a[1] * a[2] * a[3] * a[4],  # grado 6 (antes: bug, grado 5)
    ]
    extra_invariants_6 = [reynolds(seed, group) for seed in extra_seeds_6]
    all_6 = known_products_6 + extra_invariants_6
    rank6_all, _ = rank_of_span(all_6)
    print(f"Rango de {{Q^3, C^2, E4*Q}} + 5 semillas nuevas de Reynolds: {rank6_all}")
    if rank6_all > rank6_known:
        print(f"=> Aparecen {rank6_all - rank6_known} dimension(es) NUEVAS no cubiertas por productos")
        print("   de Q, C, E4 -- evidencia de que hace falta AL MENOS un generador nuevo en grado 6")
        print("   (o en un grado intermedio no explorado).")
    if rank6_all < 5:
        print(f"   (todavia no se alcanzo la dimension 5 predicha por Molien con estas semillas;")
        print(f"    puede hacer falta explorar mas semillas o hay generadores en grados intermedios)")
    elif rank6_all == 5:
        print("   => Se alcanzo EXACTAMENTE la dimension 5 predicha por Molien.")
