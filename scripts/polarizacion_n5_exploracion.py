"""
Fase 2, paso 8 (continuacion): derivar explicitamente la parametrizacion de los 2 parametros
gauge-invariantes por particula para fotones de 5 puntos, generalizando la ecuacion (2.11) del
paper original (arXiv:1910.14392, que lo hace para n=4 con UN parametro por particula).

Estrategia: construir cinematica explicita de 5 momentos sin masa en D dimensiones (D=6, generico,
evita degeneraciones de bajo D), y verificar numericamente:
  1. Que los vectores w_jk = p_j/s_ij - p_k/s_ik (para pares j,k entre las "otras" 4 particulas)
     son automaticamente ortogonales a p_i (misma identidad que en n=4).
  2. Cual es el rango real del espacio que generan estos w_jk (se espera n-3=2 para n=5, no mas).
  3. Encontrar una base explicita de 2 vectores independientes (analogos a la formula 2.11).
"""

from fractions import Fraction
import itertools
import random

import sympy as sp


def build_5pt_kinematics(D=6, seed=0):
    """
    Construye 5 momentos nulos en D dimensiones con conservacion de momento, usando el metodo
    estandar: elegir 4 momentos nulos "independientes" al azar (con componentes racionales
    pequenas para poder verificar todo exactamente) y definir el 5to como menos la suma -- pero
    hay que corregir para que el 5to tambien sea nulo. Se usa un metodo iterativo simple:
    parametrizar cada momento nulo como (E, vec n * E) con n unitario, y ajustar numericamente.

    Para simplificar y mantener todo EXACTO (racional), se usa el truco de "espacio de momento nulo"
    vía parametros de Gross-Mende/BCFW no hace falta aca: alcanza con verificar las identidades de
    ortogonalidad ALGEBRAICAMENTE (sin necesitar cinematica fisica concreta), ya que la identidad
    p_i . (p_j/s_ij - p_k/s_ik) = 0 se sigue PURAMENTE de las definiciones s_ij = -2 p_i.p_j y no
    depende de que la cinematica sea fisica. Por eso este script trabaja simbolicamente con los
    p_i.p_j como variables libres (sujetas solo a las relaciones de conservacion de momento), en
    vez de construir vectores concretos.
    """
    pass


# ---------------------------------------------------------------------------------------------
# Enfoque simbolico: tratar los productos p_i.p_j (i<j) como las UNICAS variables (ya sabemos,
# de la Fase 2 anterior, que hay exactamente 5 independientes para n=5). No hace falta cinematica
# vectorial explicita para verificar las identidades de ortogonalidad de los w_jk, porque esas
# identidades solo usan productos punto p_i.p_j, nunca componentes de los vectores en si.
# ---------------------------------------------------------------------------------------------

# Notacion: dot[i,j] = p_i . p_j para i<j (i,j en 1..5). s_ij = -2*dot[i,j].
n = 5
idx_pairs = list(itertools.combinations(range(1, n + 1), 2))
dot = {p: sp.symbols(f"d{p[0]}{p[1]}") for p in idx_pairs}


def dot_ij(i, j):
    if i == j:
        return sp.Integer(0)  # p_i^2 = 0 (sin masa)
    key = (min(i, j), max(i, j))
    return dot[key]


def s(i, j):
    return -2 * dot_ij(i, j)


# Conservacion de momento: para cada k, suma_{j != k} p_k.p_j = 0 (ya usado en Fase 2 para derivar
# la representacion de S5). Aplicamos estas 5 relaciones para reducir a 5 variables independientes.
momentum_conservation_eqs = []
for k in range(1, n + 1):
    momentum_conservation_eqs.append(sum(dot_ij(k, j) for j in range(1, n + 1) if j != k))

# Resolver: expresar 5 de los 10 dot[i,j] en terminos de los otros 5 (elegimos resolver por
# d15,d25,d35,d45 y uno mas, dejando como libres d12,d13,d14,d23,d24 -- una eleccion valida entre
# varias).
free_vars = [dot[(1, 2)], dot[(1, 3)], dot[(1, 4)], dot[(2, 3)], dot[(2, 4)]]
dependent_vars = [v for v in dot.values() if v not in free_vars]

sol = sp.solve(momentum_conservation_eqs, dependent_vars, dict=True)
assert len(sol) == 1, f"se esperaba una solucion unica, se obtuvieron {len(sol)}"
sol = sol[0]

print("Relaciones de conservacion de momento resueltas (variables dependientes en terminos de las 5 libres):")
for k, v in sol.items():
    print(f"  {k} = {v}")
print()


def dot_ij_reduced(i, j):
    """p_i.p_j ya expresado solo en terminos de las 5 variables libres."""
    raw = dot_ij(i, j)
    if raw == 0:
        return sp.Integer(0)
    return raw.subs(sol) if raw in sol else raw


def s_reduced(i, j):
    return -2 * dot_ij_reduced(i, j)


# ---------------------------------------------------------------------------------------------
# Verificar la identidad de ortogonalidad p_1 . w_jk = 0 donde w_jk "=" p_j/s_1j - p_k/s_1k
# Como no tenemos componentes vectoriales explicitas, verificamos la identidad ESCALAR:
# p_1.p_j / s_1j - p_1.p_k / s_1k = 0  (que es lo que significa "p_1 . w_jk = 0" cuando w_jk se
# define como esa combinacion lineal de p_j y p_k).
# ---------------------------------------------------------------------------------------------
print("=== Verificar ortogonalidad p_1 . (p_j/s_1j - p_k/s_1k) = 0 para pares j,k en {2,3,4,5} ===")
others = [2, 3, 4, 5]
all_ok = True
for j, k in itertools.combinations(others, 2):
    lhs = sp.simplify(dot_ij_reduced(1, j) / s_reduced(1, j) - dot_ij_reduced(1, k) / s_reduced(1, k))
    ok = (lhs == 0)
    all_ok = all_ok and ok
    print(f"  j={j}, k={k}: p1.(pj/s1j - pk/s1k) = {lhs}  {'OK (=0)' if ok else 'FALLA'}")
print(f"=> Identidad de ortogonalidad se cumple para TODOS los pares: {all_ok}\n")
print("(Nota: esta identidad es puramente algebraica -- p1.pj/s1j = p1.pj/(-2 p1.pj) = -1/2 "
      "siempre, no depende de la conservacion de momento. La pregunta real es el RANGO del "
      "espacio que generan los w_jk como VECTORES, no esta identidad escalar trivial.)\n")

# ---------------------------------------------------------------------------------------------
# La pregunta real: rango del espacio generado por los w_jk = p_j/s_1j - p_k/s_1k, EXPRESADOS
# como combinacion lineal de {p_2,p_3,p_4,p_5} (que son genericamente linealmente independientes,
# ya que junto con p_1 abarcan el plano de 4 dimensiones). El rango de los COEFICIENTES en esa
# base es exactamente el rango real de los w_jk como vectores.
# ---------------------------------------------------------------------------------------------
print("=== Rango real del espacio de los w_jk (como combinaciones de p2,p3,p4,p5) ===")
basis_order = [2, 3, 4, 5]  # basis of p2,p3,p4,p5
pairs = list(itertools.combinations(others, 2))
rows = []
for (j, k) in pairs:
    row = [sp.Integer(0)] * 4
    row[basis_order.index(j)] = 1 / s_reduced(1, j)
    row[basis_order.index(k)] = -1 / s_reduced(1, k)
    rows.append(row)

M = sp.Matrix(rows)
M = sp.simplify(M)
rank = M.rank()
print(f"Numero de vectores w_jk construidos (pares entre {{2,3,4,5}}): {len(pairs)}")
print(f"Rango real (numero de w_jk independientes), como combinaciones de {{p2,p3,p4,p5}}: {rank}")
print(f"(Nota historica: se esperaba ingenuamente n-3=2; el resultado real es 3 -- ver "
      f"12-formula-explicita-polarizacion-n5.md, error de conteo detectado y corregido.)\n")

# p1 = -(p2+p3+p4+p5) por conservacion de momento -- su fila en la base {p2,p3,p4,p5} es (-1,-1,-1,-1).
p1_row = [sp.Integer(-1)] * 4
M_with_p1 = sp.Matrix(rows + [p1_row])
rank_with_p1 = sp.simplify(M_with_p1).rank()
print(f"Rango incluyendo p1 en el span: {rank_with_p1} "
      f"({'p1 YA esta contenido en el span de los w_jk' if rank_with_p1 == rank else 'p1 agrega una dimension nueva'})")
print(f"Dimension fisica (tras cocientar por la direccion de gauge p1): {rank} - 1 = {rank - 1}\n")

assert rank == 3 and rank_with_p1 == 3, "rango inesperado -- revisar antes de confiar en la base explicita de abajo"

# Base explicita de dimension 3: {w_23, w_24, p_1} (la misma que usa la formula final del paper,
# fase2-objetivo-tecnico/12-formula-explicita-polarizacion-n5.md).
print("Base explicita (3 vectores, coeficientes sobre p2,p3,p4,p5): {w_23, w_24, p_1}")
w23 = rows[pairs.index((2, 3))]
w24 = rows[pairs.index((2, 4))]
base_candidate = sp.Matrix([w23, w24, p1_row])
print(f"  w_23 = {w23}")
print(f"  w_24 = {w24}")
print(f"  p_1  = {p1_row}")
base_rank = base_candidate.rank()
print(f"  rango de {{w_23, w_24, p_1}}: {base_rank} (debe ser 3 para que sea una base completa)\n")

assert base_rank == 3, "{w_23, w_24, p_1} no es una base completa -- no confiar en la descomposicion de abajo"

print("Expresando los demas w_jk como combinacion lineal de w_23, w_24, p_1:")
c1, c2, c3 = sp.symbols("c1 c2 c3")
for (j, k), row in zip(pairs, rows):
    if (j, k) in [(2, 3), (2, 4)]:
        continue
    eqs = [sp.Eq(c1 * w23[m] + c2 * w24[m] + c3 * p1_row[m], row[m]) for m in range(4)]
    sol_c = sp.solve(eqs, [c1, c2, c3], dict=True)
    assert len(sol_c) == 1, f"w_{j}{k}: se esperaba una unica solucion, se obtuvo {sol_c}"
    print(f"    w_{j}{k} = {sol_c[0]}")
print("\n=> Todos los w_jk restantes se expresan exactamente como combinacion de {w_23, w_24, p_1} "
      "-- confirma que esa base captura toda la informacion, sin perder ni sobrar nada "
      "(reproduce la fórmula final citada en el paper).")
