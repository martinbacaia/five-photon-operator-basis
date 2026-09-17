"""
Cross-verificacion (version corregida) entre:
  (A) la representacion de S4 sobre invariantes de Mandelstam de n=4, calculada como el nucleo del
      mapa de incidencia de K4 (n5_mandelstam_representation.py), restringida al subgrupo S3 que
      fija la particula 4;
  (B) la representacion "2M" de S3 sobre (s,t) ya validada en s3_partition_functions.py /
      verificacion_rigurosa.py.

Primer intento (ver historial): se comparo usando una correspondencia ADIVINADA entre las
permutaciones fisicas de {1,2,3,4} y las etiquetas abstractas "(12)","(13)","(23)","(123)","(132)"
del script anterior. Las trazas y determinantes coincidian para los 6 elementos (lo cual NO prueba
nada por si solo, ya que traza/det son constantes en cada clase de conjugacion y las clases de S3
solo tienen 3 tipos), pero al buscar una matriz de cambio de base P que conjugue TODOS los
generadores simultaneamente, fallo para (13) y (23) -- señal de que la correspondencia entre
etiquetas abstractas y permutaciones fisicas estaba mal, no que la representacion estuviera mal.

Esta version deriva la accion fisica DESDE CERO, sin adivinar nombres: usando la conservacion de
momento (s12=s34, s13=s24, s14=s23=-s12-s13, exactamente como se derivo a mano en esta sesion) se
calcula como cada permutacion fisica de {1,2,3,4} que fija la particula 4 transforma el par
(s12, s13). Esto se compara despues, con una matriz de cambio de base P resuelta explicitamente y
verificada en los 6 elementos, contra la representacion nueva basada en el nucleo de incidencia.
"""

import sympy as sp

from n5_mandelstam_representation import kernel_basis_matrix, action_on_kernel

s, t = sp.symbols("s t")
u = -s - t

# pares (i,j) -> expresion en (s,t) via la conservacion de momento para n=4:
# s12=s34=s ; s13=s24=t ; s14=s23=u=-s-t   (derivado a mano y verificado en esta sesion)
raw_sij = {
    (1, 2): s, (2, 1): s,
    (3, 4): s, (4, 3): s,
    (1, 3): t, (3, 1): t,
    (2, 4): t, (4, 2): t,
    (1, 4): u, (4, 1): u,
    (2, 3): u, (3, 2): u,
}


def physical_action(sigma):
    """
    sigma: tupla de longitud 4, sigma[i-1] = imagen de la particula i.
    Devuelve (new_s, new_t) = valor de (s12, s13) en la configuracion permutada, en terminos de
    las variables (s,t) de la configuracion original.

    NOTA DE UN BUG REAL ENCONTRADO Y CORREGIDO (ver 02-correccion-convencion-homomorfismo.md):
    la primera version de esta funcion usaba new_s_ij = old_s_{sigma^-1(i),sigma^-1(j)}, que
    resulto ser un ANTI-homomorfismo (A(g1)A(g2) = A(g2 o g1), no A(g1 g2)) -- verificado
    explicitamente por computo, no supuesto. La convencion correcta, consistente con
    action_on_kernel (que es P(sigma) e_i = e_{sigma(i)}, un homomorfismo genuino por
    construccion), es new_s_ij = old_s_{sigma(i),sigma(j)} (sigma, no su inversa).
    """
    new_s = raw_sij[(sigma[0], sigma[1])]   # s_{sigma(1), sigma(2)}
    new_t = raw_sij[(sigma[0], sigma[2])]   # s_{sigma(1), sigma(3)}
    return new_s, new_t


def physical_matrix(sigma):
    im_s, im_t = physical_action(sigma)
    col_s = sp.Matrix([[im_s.coeff(s)], [im_s.coeff(t)]])
    col_t = sp.Matrix([[im_t.coeff(s)], [im_t.coeff(t)]])
    return sp.Matrix.hstack(col_s, col_t)


perms_fixing_4 = {
    "e": (1, 2, 3, 4),
    "(12)": (2, 1, 3, 4),
    "(13)": (3, 2, 1, 4),
    "(23)": (1, 3, 2, 4),
    "(123)": (2, 3, 1, 4),
    "(132)": (3, 1, 2, 4),
}

if __name__ == "__main__":
    basis, edge_list = kernel_basis_matrix(4)

    print("Accion fisica derivada de la conservacion de momento (s,t)=(s12,s13):")
    matrices_physical = {}
    for name, sigma in perms_fixing_4.items():
        Mphys = physical_matrix(sigma)
        matrices_physical[name] = Mphys
        print(f"  {name}: (s,t) -> {physical_action(sigma)}   matriz={Mphys.tolist()}")

    print("\nComparacion traza/det (necesario pero no suficiente):")
    matrices_new = {}
    for name, sigma in perms_fixing_4.items():
        A_new = action_on_kernel(4, basis, edge_list, sigma)
        matrices_new[name] = A_new
        A_phys = matrices_physical[name]
        print(f"  {name}: traza(nucleo)={A_new.trace()} traza(fisica)={A_phys.trace()}  "
              f"det(nucleo)={A_new.det()} det(fisica)={A_phys.det()}")

    # Resolver P tal que A_phys(g) = P^-1 A_new(g) P para TODOS los generadores simultaneamente
    p11, p12, p21, p22 = sp.symbols("p11 p12 p21 p22")
    P = sp.Matrix([[p11, p12], [p21, p22]])
    eqs = []
    for name in ["(12)", "(123)"]:  # (12) y (123) generan S3
        diff = matrices_physical[name] * P - P * matrices_new[name]
        eqs.extend([diff[0, 0], diff[0, 1], diff[1, 0], diff[1, 1]])
    sol = sp.solve(eqs, [p11, p12, p21, p22], dict=True)
    print(f"\nSoluciones para P (resuelta con (12) y (123)): {sol}")

    if sol:
        Pmat = sp.Matrix([[sol[0].get(p11, p11), sol[0].get(p12, p12)],
                           [sol[0].get(p21, p21), sol[0].get(p22, p22)]])
        free_syms = list(Pmat.free_symbols)
        Pmat_num = Pmat.subs({sym: 1 for sym in free_syms}) if free_syms else Pmat
        print(f"P concreta: {Pmat_num.tolist()}  det(P)={Pmat_num.det()}")

        full_check = True
        for name in perms_fixing_4:
            lhs = matrices_physical[name] * Pmat_num
            rhs = Pmat_num * matrices_new[name]
            ok = sp.simplify(lhs - rhs) == sp.zeros(2, 2)
            full_check = full_check and ok
            print(f"  chequeo conjugacion para TODOS los 6 elementos -- {name}: "
                  f"{'OK' if ok else 'FALLA'}")
        print(f"\n=> Representaciones isomorfas (P valida para los 6 elementos, det(P)!=0): "
              f"{full_check and Pmat_num.det() != 0}")
