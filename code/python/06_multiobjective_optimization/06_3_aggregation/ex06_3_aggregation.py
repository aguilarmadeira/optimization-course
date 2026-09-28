"""Reproduz os exemplos do deck 6.3 (agregação de objetivos).

Exemplo convexo da UC (6.1):  f1 = x1^2 + x2^2, f2 = (x1 - 1)^2 + x2^2,
x em [0,1]^2, x1^2 + x2^2 <= 1.
  1. «O que fazem os pesos?» e «Pesos uniformes ≠ pontos uniformes»: soma
     ponderada com 11 pesos w1 = 0, 0.1, ..., 1, x0 = (0.5, 0.3), arranque a
     quente, SLSQP: x1* = w2, f1 = w2^2, f2 = (1 - w2)^2 (tabela); n_f = 75.
Exemplo não convexo da UC (6.2):  f1 = x1, f2 = 1 - x1^2 + x2,
x em [0,1] x [0, 0.6].
  2. «A limitação fundamental»: com w1, w2 > 0 a soma ponderada só dá os dois
     extremos (supported): 11, 100 e 10 000 pesos -> 2 pontos.
  3. «Tchebycheff ponderado»: min theta s.a. w_i (f_i - z_i^id) <= theta,
     z^id = (0, 0): encontra também os pontos non-supported.
4. «A escala importa»: projetos A, B, C; w = (0.9, 0.1) sem normalizar dá
   81.4 / 85.8 / 100.5 (ganha A); normalizado com ideal e nadir, C / B / A.

Usa aggregation.py (nesta pasta). O n_f = 75 é o do SLSQP do SciPy (como no
script das figuras) e verifica-se aqui, não no MATLAB/Octave (outro solver).
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex06_3_aggregation.py  (numpy e scipy)

Otimização — deck 6.3.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from aggregation import weighted_sum, weighted_tchebycheff


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def perto(v, s, tol):
    """v está a menos de tol do valor exato s?"""
    return bool(np.all(np.abs(np.ravel(np.asarray(v, float))
                              - np.ravel(np.asarray(s, float))) <= tol))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


limpa = lambda v: np.where(np.abs(v) >= 5e-5, v, 0.0) + 0.0     # não imprimir -0.0000
ok = True
print("Otimização — deck 6.3: métodos de agregação de objetivos")

# ------------------------------------------------ 1. Exemplo convexo
F = lambda x: np.array([x[0]**2 + x[1]**2, (x[0] - 1)**2 + x[1]**2])
cons = [{"type": "ineq", "fun": lambda x: 1 - x[0]**2 - x[1]**2}]   # forma da UC: g >= 0
b = [(0, 1), (0, 1)]
x0 = [0.5, 0.3]
w1 = np.linspace(0, 1, 11); W = np.c_[w1, 1 - w1]
print("\n1. Exemplo convexo: min w1 f1 + w2 f2, 11 pesos, x0 = (%g, %g), arranque a quente" % tuple(x0))
r = weighted_sum(F, x0, b, W, cons, verbose=True)
print("   solver: %s" % r.message)
print("   n_f = %d (soma de r.nfev; inclui as avaliações dos gradientes por diferenças finitas)" % r.nfev)
T = [[0.0, 1.0, 1.000, 1.000, 0.000], [0.1, 0.9, 0.900, 0.810, 0.010], [0.2, 0.8, 0.800, 0.640, 0.040],
     [0.3, 0.7, 0.700, 0.490, 0.090], [0.4, 0.6, 0.600, 0.360, 0.160], [0.5, 0.5, 0.500, 0.250, 0.250],
     [0.6, 0.4, 0.400, 0.160, 0.360], [0.7, 0.3, 0.300, 0.090, 0.490], [0.8, 0.2, 0.200, 0.040, 0.640],
     [0.9, 0.1, 0.100, 0.010, 0.810], [1.0, 0.0, 0.000, 0.000, 1.000]]
c = [confere(np.c_[W, r.x[:, 0], r.fx], T, 3), perto(r.x[:, 0], W[:, 1], 1e-5), perto(r.x[:, 1], 0, 1e-5),
     confere(r.x[[10, 5, 0]], [[0, 0], [0.5, 0], [1, 0]], 3), r.flag == 0, r.nfev == 75]
print("  tabela (11 pesos): %s | x1* = w2: %s | x2* = 0: %s | w = (1;0), (0.5;0.5), (0;1) -> (0;0), (0.5;0), (1;0): %s"
      " | sem falhas: %s | n_f = 75: %s" % simnao(c))
ok = ok and all(c)
d = np.diff(r.fx[:, 0])
print("   pesos uniformes, pontos não uniformes: passos em f1 de %.2f (junto a f1 = 1) a %.2f (junto a f1 = 0)"
      % (-d[0], -d[-1]))

# ------------------------------------------------ 2. Exemplo não convexo
F = lambda x: np.array([x[0], 1 - x[0]**2 + x[1]])
b = [(0, 1), (0, 0.6)]
X0 = [[0.1, 0.1], [0.5, 0.1], [0.9, 0.1]]          # três pontos iniciais (solver local)
print("\n2. Exemplo não convexo: soma ponderada, 3 pontos iniciais por peso, fica o melhor")
w1 = np.array([0.2, 0.5, 0.8])
r = weighted_sum(F, X0, b, np.c_[w1, 1 - w1])
for wk, xk, Fk in zip(w1, r.x, r.fx):
    print("   w1 = %.1f: x* = (%.4f; %.4f), (f1, f2) = (%.4f; %.4f)" % (wk, *limpa(np.r_[xk, Fk])))
ra = weighted_sum(F, [0.5, 0.1], b, [[0.5, 0.5]])
xa = ra.x[0]
print("   w1 = 0.5 só a partir de (0.5; 0.1): x = (%.4f; %.4f), f~ = %.4f"
      % (*limpa(xa), 0.5 * ra.fx[0].sum()), end="")
if 0.01 < xa[0] < 0.99:
    print("\n   (cuidado: ponto estacionário de f~, côncava em x1 — é um máximo, não um mínimo)")
else:
    print(" (um extremo)")
c = []
for N in (11, 100, 10000):
    w1 = np.linspace(0, 1, N)
    pos = (w1 > 0) & (w1 < 1)                       # w1, w2 > 0
    if N < 10000:
        Fx = weighted_sum(F, X0, b, np.c_[w1, 1 - w1]).fx
        como = "solver"
    else:                                           # 10 000 pesos: pesquisa exaustiva em x1
        t = np.linspace(0, 1, 1001)                 # (com w2 > 0 o ótimo tem x2 = 0)
        wp = w1[pos][:, None]
        j = np.argmin(wp * t + (1 - wp) * (1 - t**2), axis=1)
        Fx = np.zeros((N, 2)); Fx[pos] = np.c_[t[j], 1 - t[j]**2]
        como = "grelha"
    P = np.unique(np.round(Fx[pos] * 1e6) / 1e6, axis=0) + 0.0
    print("   %5d pesos (%s), %5d com w1, w2 > 0 -> %d pontos distintos: %s"
          % (N, como, pos.sum(), len(P), " ".join("(%g; %g)" % tuple(p) for p in P)))
    c.append(len(P) == 2 and perto(P[np.lexsort(P.T[::-1])], [[0, 1], [1, 0]], 1e-6))
print("   (com w1 = 1, w2 = 0, x2 fica livre: o solver devolve (0; 1.1), só fracamente eficiente)")
print("  11 pesos -> 2 pontos: %s | 100 -> 2: %s | 10 000 -> 2: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 3. Tchebycheff ponderado
zid = np.array([0.0, 0.0])
w1 = np.linspace(0, 1, 11); W = np.maximum(np.c_[w1, 1 - w1], 1e-3)   # como no script das figuras
print("\n3. Tchebycheff ponderado (min theta), z^id = (0, 0), 11 pesos (w_i >= 1e-3), a quente")
r = weighted_tchebycheff(F, [0.5, 0.3], b, W, zid, verbose=True)
print("   solver: %s; n_f = %d chamadas a F" % (r.message, r.nfev))
a, bb = W[:, 0], W[:, 1]                             # exato: w1 f1 = w2 (1 - f1^2)
f1ex = (-a + np.sqrt(a**2 + 4 * bb**2)) / (2 * bb)
nint = int(np.sum((r.fx[:, 0] > 0.01) & (r.fx[:, 0] < 0.99)))
print("   f1 exato (w1 f1 = w2 (1 - f1^2)):" + "".join(" %.4f" % v for v in f1ex))
print("   %d dos 11 pontos são interiores (non-supported)" % nint)
c = [perto(r.fx[:, 0], f1ex, 1e-3), perto(r.fx[:, 1], 1 - f1ex**2, 2e-3), nint == 9, r.flag == 0]
print("  pontos na frente (cruzamento w1 f1 = w2 f2): %s | f2 = 1 - f1^2: %s | 9 non-supported: %s"
      " | sem falhas: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 4. A escala importa
print("\n4. A escala importa: f1 = massa (kg), f2 = custo (€)")
Y = np.array([[1.5, 800], [0.9, 850], [0.5, 1000]]); nomes = "ABC"
s = Y @ [0.9, 0.1]; j = int(np.argmin(s))
print("   w = (0.9; 0.1), sem normalizar: A = %.2f, B = %.2f, C = %.2f -> ganha %s" % (*s, nomes[j]))
zi, zn = Y.min(0), Y.max(0)                          # os três são eficientes
Yn = (Y - zi) / (zn - zi)
Wn = np.array([[0.9, 0.1], [0.5, 0.5], [0.1, 0.9]])
print("   normalizado com z^id = (%g; %g), z^nad = (%g; %g):" % (*zi, *zn))
for w in Wn:
    print("     w = (%.1f; %.1f): A = %.3f, B = %.3f, C = %.3f" % (*w, *(Yn @ w)))
esc = "".join(nomes[i] for i in np.argmin(Yn @ Wn.T, axis=0))
print("     escolhe %s" % " / ".join(esc))
c = [confere(s, [81.4, 85.8, 100.5], 1), j == 0, esc == "CBA"]
print("  81.4 / 85.8 / 100.5: %s | ganha A: %s | normalizado C / B / A: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
