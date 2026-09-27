"""Reproduz os exemplos do deck 6.4 (epsilon-constrangimento).

Exemplo não convexo da UC (6.2):  min (f1, f2) = (x1, 1 - x1^2 + x2),
x em [0,1] x [0, 0.6].  Subproblema:  min f2  s.a.  f1 <= e.
  1. «O mesmo exemplo do 6.3, à mão»: e = 0.25, 0.50, 0.75.
  2. «Em MATLAB: varrer e com fmincon» / «O mesmo em Python»: 10 valores
     e = 0.05, 0.15, ..., 0.95, x0 = (0.5, 0.3), arranque a quente, SLSQP:
     10 pontos, 60 avaliações de f2, 90 de f1 (150 chamadas), u = 2e.
  3. «O multiplicador: a taxa de troca»: e = 0.5 -> u = 1; f2(0.51) = 0.7399;
     e = 0.9 -> u = 1.8, a mesma folga compra ~0.018.
  4. «Restrição inativa»: e = 1.05 e 1.2 -> a mesma solução (1, 0), u = 0.
  5. «Nantes–Lille»: o mais barato com duração <= e (e = 4, 5, 10 h).

As contagens 60/90/150 são as do SLSQP do SciPy (como no script das figuras,
com a restrição na forma da UC, {'type': 'ineq', 'fun': e - f1}); verificam-se
aqui e não no MATLAB/Octave (outro solver).
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex06_4_epsilon_constraint.py  (numpy e scipy)

Otimização — deck 6.4.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import warnings

import numpy as np
from scipy.optimize import NonlinearConstraint, minimize  # NonlinearConstraint: trust-constr (r.v)

from epsilon_constraint import epsilon_constraint

# o quasi-Newton do trust-constr avisa que f1 é linear («delta_grad == 0.0»); é esperado
warnings.filterwarnings("ignore", message="delta_grad == 0.0")


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


ok = True
print("Otimização — deck 6.4: método do epsilon-constrangimento")

f1 = lambda x: x[0]
f2 = lambda x: 1 - x[0]**2 + x[1]
b = [(0, 1), (0, 0.6)]
x0 = [0.5, 0.3]

# ------------------------------------------------ 1. À mão: três valores
print("\n1. min f2 = 1 - x1^2 + x2  s.a.  f1 = x1 <= e  (x* = (e, 0))")
r = epsilon_constraint(f1, f2, x0, b, [0.25, 0.50, 0.75])
for Fk in r.fx:
    print("   e = %.2f  =>  (f1, f2) = (%.4f; %.4f)" % (Fk[0], Fk[0], Fk[1]))
c = confere(r.fx, [[0.25, 0.9375], [0.50, 0.75], [0.75, 0.4375]], 4)
print("  (0.25; 0.9375), (0.50; 0.75), (0.75; 0.4375): %s" % simnao([c]))
ok = ok and c

# ------------------------------------------------ 2. Varrimento do slide
print("\n2. Varrimento e = 0.05, 0.15, ..., 0.95 a partir de x0 = (%g, %g), arranque a quente" % tuple(x0))
e = np.linspace(0.05, 0.95, 10)
r = epsilon_constraint(f1, f2, x0, b, e, verbose=True)
print("   solver: %s" % r.message)
nd = len(np.unique(np.round(r.fx * 1e6), axis=0))
print("   %d valores de e -> %d pontos distintos, todos na frente f2 = 1 - f1^2" % (len(e), nd))
print("   contagens (SLSQP): %d avaliações de f2 (soma de r.nfev), %d de f1 (%d chamadas); %d iterações"
      % (r.nfev, r.nfev1, r.ncalls, r.nit))
c = [len(e) == 10, nd == 10, perto(r.x, np.c_[e, np.zeros(10)], 1e-6),
     perto(r.fx[:, 1], 1 - e**2, 1e-6), perto(r.u, 2 * e, 1e-4), r.flag == 0,
     r.nfev == 60, r.nfev1 == 90, r.ncalls == 150]
print("  10 valores: %s | 10 pontos: %s | x* = (e, 0): %s | f2 = 1 - e^2: %s | u = 2e: %s | sem falhas: %s"
      " | 60 de f2: %s | 90 de f1: %s | 150 chamadas: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 3. O multiplicador
print("\n3. Multiplicador de f1 <= e (g = e - f1 >= 0): u = 2e = -d f2*/de")
r = epsilon_constraint(f1, f2, x0, b, [0.50, 0.51, 0.90, 0.91])
d1 = r.fx[1, 1] - r.fx[0, 1]; d2 = r.fx[2, 1] - r.fx[3, 1]
print("   e = 0.50: x* = (%.4f; %.4f), f2* = %.4f, u = %.4f" % (*r.x[0], r.fx[0, 1], r.u[0]))
print("   e = 0.51: f2* = %.4f; diferença %.4f  (aproximação -u*0.01 = %.4f)"
      % (r.fx[1, 1], d1, -r.u[0] * 0.01))
print("   e = 0.90: u = %.4f; a folga 0.90 -> 0.91 compra %.4f em f2" % (r.u[2], d2))
rt = minimize(f2, x0, bounds=b, constraints=NonlinearConstraint(f1, -np.inf, 0.5), method="trust-constr")
print("   e = 0.50 com trust-constr (como no slide): u = r.v[0] = %.4f (x* = (%.4f; %.4f))"
      % (rt.v[0][0], *rt.x))
c = [confere(r.x[0, 0], 0.5, 1), confere(r.fx[0, 1], 0.75, 2), confere(r.u[0], 1, 3),
     confere(r.fx[1, 1], 0.7399, 4), confere(d1, -0.0101, 4), confere(r.u[2], 1.8, 3),
     confere(d2, 0.018, 3), confere(rt.v[0][0], 1, 2)]
print("  x* = 0.5: %s | f2* = 0.75: %s | u = 1: %s | f2(0.51) = 0.7399: %s | -0.0101: %s"
      " | u(0.9) = 1.8: %s | ~0.018: %s | trust-constr u ~ 1: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 4. Restrição inativa
print("\n4. Restrição inativa: e > 1")
r = epsilon_constraint(f1, f2, x0, b, [1.05, 1.2])
for h, xk, Fk, uk in zip(r.history, r.x, r.fx, r.u):
    print("   e = %.2f: x* = (%.4f; %.4f), f1 = %.4f < e, u = %.4f" % (h[0], *xk, Fk[0], abs(uk)))
c = perto(r.x, [[1, 0], [1, 0]], 1e-6) and perto(r.u, [0, 0], 1e-6)
print("  mesma solução (1, 0) com u = 0: %s" % simnao([c]))
ok = ok and c

# ------------------------------------------------ 5. Nantes–Lille (6.1)
print("\n5. Nantes–Lille: o mais barato com duração <= e")
nomes = ["Carro via Paris", "Carro via Rouen", "Carro sem portagens", "Autocarro 1",
         "Autocarro 2", "TGV", "TGV low-cost", "Avião"]
dur = np.array([5 + 33 / 60, 5 + 49 / 60, 7 + 18 / 60, 9.5, 9.5, 4 + 2 / 60, 4 + 20 / 60, 1 + 5 / 60])
custo = np.array([230, 251, 204, 38, 28, 84, 40, 134])
esc = []
for E in (4, 5, 10):
    adm = np.flatnonzero(dur <= E)
    j = adm[np.argmin(custo[adm])]              # empate: fica o primeiro
    esc.append(j)
    print("   e = %2d h: %s (%d €)" % (E, nomes[j], custo[j]))
c = esc == [7, 6, 4] and list(custo[esc]) == [134, 40, 28]
print("  avião (134), TGV low-cost (40), autocarro 2 (28): %s" % simnao([c]))
ok = ok and c

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
