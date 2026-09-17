"""
Fase 2, paso 3: determinar si el anillo de polinomios S5-invariantes sobre la representacion de
5 dimensiones de los invariantes de Mandelstam (ver n5_mandelstam_representation.py) es un anillo
de polinomios LIBRE (como el de S3 sobre s,t,u -- ver s3_partition_functions.py) o si tiene
relaciones no triviales.

Criterio usado: el teorema de Chevalley-Shephard-Todd dice que el anillo de invariantes de un grupo
finito G actuando en un espacio vectorial V es un anillo de polinomios (libremente generado) SI Y
SOLO SI G esta generado, EN ESA REPRESENTACION ESPECIFICA, por pseudo-reflexiones (elementos que
fijan un hiperplano, es decir, cuyo espacio propio de autovalor 1 tiene dimension = dim(V)-1).

Para S3 actuando sobre (s,t) (la "2M"), las 3 transposiciones SI son reflexiones genuinas (fijan una
recta, autovalores (1,-1)) -- por eso el anillo de invariantes de esa representacion es libre
(generado por 2 invariantes, de grado 2 y 3: s^2+t^2+u^2 y stu).

Para S5 actuando sobre el espacio de 5 invariantes de Mandelstam independientes, se verifica aqui
si ALGUN elemento de S5 actua como reflexion en esa representacion especifica (no en la
representacion estandar de permutacion, que es una cosa distinta).
"""

import sympy as sp

from n5_mandelstam_representation import kernel_basis_matrix, action_on_kernel, class_representatives


def eigenvalue_multiplicities(A):
    """Diccionario {autovalor: multiplicidad algebraica} de una matriz sympy exacta."""
    return A.eigenvals()


def is_reflection(A, dim):
    """
    Una pseudo-reflexion (en una representacion REAL) tiene autovalor 1 con multiplicidad
    exactamente dim-1, y un unico autovalor restante que (por ser matriz real de orden finito)
    debe ser -1 (los autovalores complejos vienen en pares conjugados; si solo uno no es 1, tiene
    que ser real, y como no es 1, tiene que ser -1).
    """
    ev = eigenvalue_multiplicities(A)
    mult_one = ev.get(sp.Integer(1), 0)
    return mult_one == dim - 1


if __name__ == "__main__":
    n = 5
    basis, edge_list = kernel_basis_matrix(n)
    dim = basis.shape[1]
    reps, class_sizes = class_representatives(n)

    print(f"=== Buscando reflexiones de S{n} en la representacion de {dim} dimensiones ===\n")
    print(f"{'tipo de ciclo':>18} | {'tamano clase':>12} | {'traza':>6} | {'autovalores (mult)':>40} | es reflexion?")

    any_reflection = False
    for ct, sigma in sorted(reps.items(), key=lambda kv: (-len(kv[0]), kv[0])):
        A = action_on_kernel(n, basis, edge_list, sigma)
        ev = eigenvalue_multiplicities(A)
        ev_str = ", ".join(f"{val}(x{mult})" for val, mult in ev.items())
        refl = is_reflection(A, dim)
        any_reflection = any_reflection or (refl and ct != tuple([1] * n))
        print(f"{str(ct):>18} | {class_sizes[ct]:>12} | {A.trace()!s:>6} | {ev_str:>40} | "
              f"{'SI' if refl else 'no'}")

    print()
    if any_reflection:
        print("Se encontro al menos una reflexion (fuera de la identidad) -> el anillo de "
              "invariantes PODRIA ser libre (Chevalley-Shephard-Todd no lo garantiza sin generar "
              "el grupo completo con reflexiones, pero es un indicio a favor).")
    else:
        print("NINGUN elemento no-trivial de S5 actua como reflexion (autovalor 1 con "
              "multiplicidad dim-1=4) en esta representacion de 5 dimensiones.")
        print("=> Por Chevalley-Shephard-Todd (direccion 'solo si'), el anillo de polinomios "
              "S5-invariantes sobre este espacio de 5 invariantes de Mandelstam NO es un anillo "
              "de polinomios libre -- tiene relaciones no triviales entre sus generadores.")
        print("Esto es CUALITATIVAMENTE distinto del caso n=4 (donde S3 SI actua por reflexiones")
        print("sobre (s,t), y el anillo de invariantes es libre, generado en grado 2 y 3).")

    # -------------------------------------------------------------------------------------------
    # Serie de Molien: el Hilbert series del anillo de invariantes es calculable SIEMPRE (sea o no
    # libre), via la formula de Molien: H(x) = (1/|G|) sum_g 1/det(I - x*rho(g))
    # -------------------------------------------------------------------------------------------
    print("\n=== Serie de Molien del anillo de invariantes de S5 sobre esta representacion ===")
    x = sp.symbols("x")
    I = sp.eye(dim)
    total = sp.Integer(0)
    from math import factorial
    for ct, sigma in reps.items():
        A = action_on_kernel(n, basis, edge_list, sigma)
        detmat = (I - x * A).det()
        term = 1 / detmat
        total += class_sizes[ct] * term
    molien = sp.simplify(total / factorial(n))
    molien_series = sp.series(molien, x, 0, 12).removeO()
    print("H(x) truncada a orden 11:")
    sp.pprint(sp.expand(molien_series))

    coeffs = sp.Poly(molien_series, x).all_coeffs()[::-1]
    print("\nDimension del espacio de invariantes de grado d (d=0,1,2,...):")
    for d, c in enumerate(coeffs):
        print(f"  d={d}: {c}")
