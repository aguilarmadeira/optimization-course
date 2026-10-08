"""Reproduz os exemplos do deck 4.2 (condições KKT).

(a) Resolve os exemplos dos slides pelas condições KKT, com as contas do
    slide, e verifica-as numericamente com kkt_check (estacionariedade,
    admissibilidade, complementaridade u_j g_j = 0 e sinal u_j >= 0).
    Convenção da UC: g_j(x) >= 0, h_l(x) = 0, L = f - u'g - lambda'h.
    Nas aulas: L = f - lambda'g - beta'h; no código, u = lambda_j e lambda = beta_l das aulas.
      1. De 4.1 para 4.2 (círculo largo / apertado) e Geometria (2) (u < 0).
      2. Exemplo 1: os 2^J = 4 casos da receita; sobrevive o caso 3.
      3. «Cuidado: ativa não implica u > 0»: min x^2 s.a. x >= 0.
      4. Exemplo 2: verificar o candidato (2,1): u1 = 1/3, u2 = 2/3.
      5. «KKT dá candidatos»: min -x^2 s.a. x + 1 >= 0, 2 - x >= 0.
      6. O que medem os multiplicadores: exemplo 1 com b = 3.1 e a produção
         do 1.1 (preços-sombra 10, 20, 0).
      7. Quando KKT falha: LICQ.
      8. Exercícios 1, 2 e 4 (as respostas a cinzento).
(b) Resolve o exemplo 1, o exemplo 2 e os exercícios 2 e 4 com
    scipy.optimize.minimize(method='trust-constr'), como no slide
    «No computador»; os multiplicadores vêm em r.v, um vetor por restrição.
    Sinal: o SciPy escreve L = f + v'c. Com NonlinearConstraint(g, 0, np.inf)
    (g >= 0) e NonlinearConstraint(h, 0, 0):  u_UC = -r.v[0],  lambda_UC = -r.v[1].
    A produção com scipy.optimize.linprog (como no 1.1): os preços-sombra são
    -res.ineqlin.marginals (o SciPy minimiza -L e dá derivadas em ordem a b).

Usa kkt_check.py (nesta pasta). Os slides usam vírgula decimal; aqui o ponto.
Não há função da UC neste deck: é um exemplo. Implementação didática.
Correr (de qualquer pasta):  python ex04_2_kkt.py  (numpy e scipy)

Otimização — deck 4.2.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import warnings

import numpy as np
from scipy.optimize import NonlinearConstraint, linprog, minimize

from kkt_check import kkt_check

# Com restrições lineares, o quasi-Newton do trust-constr avisa que a Hessiana
# da restrição é nula («delta_grad == 0.0»); é esperado e não afeta o resultado.
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


sn = lambda b: "sim" if b else "não"
vec = lambda v: ", ".join("%.4f" % (t + 0.0) for t in np.ravel(v))
ok = True
tolA = 1e-10        # contas exatas (a)
tolB = 1e-4         # solver (b)

print("Otimização — deck 4.2: condições de Karush-Kuhn-Tucker")
print("Convenção da UC: g_j >= 0, h_l = 0, L = f - u'g - lambda'h  (nas aulas: u = lambda_j, lambda = beta_l)")
print("\n=============== (a) pelas condições KKT ===============")

# ------------------------------------------------ 1. De 4.1 para 4.2
df = lambda x: 2 * (np.asarray(x) - 1)
dg = lambda x: np.array([[-2 * x[0], -2 * x[1]]])        # grad g' para g = R2 - x1^2 - x2^2
ucalc = lambda x: np.linalg.lstsq(dg(x).T, df(x), rcond=None)[0]   # grad f = u grad g
print("\n1. min (x1-1)^2 + (x2-1)^2  s.a.  g = R^2 - x1^2 - x2^2 >= 0")
x = np.array([1.0, 1.0]); g = 4 - x @ x
c1, r = kkt_check(df(x), [g], dg(x), [0.0], tol=tolA)
print("   largo (R^2 = 4):    x* = (1, 1), g = %g > 0 (inativa), u* = 0, grad f = 0;  KKT: %s"
      % (g, sn(c1)))
x = np.array([1.0, 1.0]) / np.sqrt(2); g = 1 - x @ x; uap = ucalc(x)
c2, r = kkt_check(df(x), [g], dg(x), uap, tol=tolA)
print("   apertado (R^2 = 1): x* = (%.4f, %.4f), g = %.1e (ativa), u* = %.4f = sqrt(2) - 1;  KKT: %s"
      % (*x, g, uap[0], sn(c2)))
x = np.array([1.0, 1.0]) * np.sqrt(2); g = 4 - x @ x; u = ucalc(x)
c3, r = kkt_check(df(x), [g], dg(x), u, tol=tolA)
print("   Geometria (2): x = (sqrt 2, sqrt 2) no largo: u = %.4f = 1/sqrt(2) - 1 < 0;  KKT: %s "
      "(falha o sinal: %.4f)" % (u[0], sn(c3), r.sinal))
c = [c1, c2, perto(uap, np.sqrt(2) - 1, tolA), perto(u, 1 / np.sqrt(2) - 1, tolA), not c3]
print("   largo KKT: %s | apertado KKT: %s | u* = sqrt2 - 1: %s | u = 1/sqrt2 - 1: %s"
      " | (sqrt2, sqrt2) não é KKT: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 2. Exemplo 1: 2^J casos
f1 = lambda x: (x[0] - 2)**2 + 2 * (x[1] - 1)**2
df1 = lambda x: np.array([2 * (x[0] - 2), 4 * (x[1] - 1)])
g1 = lambda x: np.array([3 - x[0] - 4 * x[1], x[0] - x[1]])
Jg1 = np.array([[-1.0, -4], [1, -1]])
print("\n2. Exemplo 1: min (x1-2)^2 + 2(x2-1)^2  s.a.  g1 = 3 - x1 - 4x2 >= 0,  g2 = x1 - x2 >= 0")
print("   2(x1-2) + u1 - u2 = 0,  4(x2-1) + 4u1 + u2 = 0,  u1 g1 = 0,  u2 g2 = 0,  u >= 0")
casos = [(0, 0), (0, 1), (1, 0), (1, 1)]      # 1 = restrição ativa (g_j = 0); 0 = u_j = 0
hip = ["u1 = u2 = 0     ", "u1 = 0, g2 = 0  ", "g1 = 0, u2 = 0  ", "g1 = g2 = 0     "]
Xs = [(2, 1), (4/3, 4/3), (5/3, 1/3), (3/5, 3/5)]          # slide
Us = [(0, 0), (0, -4/3), (2/3, 0), (22/25, -48/25)]


def caso_ex1(a1, a2, b=3.0):
    """Sistema linear do caso (a1, a2) em (x1, x2, u1, u2)."""
    M = [[2, 0, 1, -1], [0, 4, 4, 1]]; rhs = [4, 4]        # estacionariedade
    M.append([1, 4, 0, 0] if a1 else [0, 0, 1, 0]); rhs.append(b if a1 else 0)
    M.append([1, -1, 0, 0] if a2 else [0, 0, 0, 1]); rhs.append(0)
    z = np.linalg.solve(np.array(M, float), np.array(rhs, float))
    return z[:2], z[2:]


print("   %-5s %-17s %-19s %-19s %-16s %s" % ("caso", "hipótese", "x", "u", "g", "verificação"))
c = []
for k, (a1, a2) in enumerate(casos):
    x, u = caso_ex1(a1, a2); g = g1(x)
    ver = ["g%d < 0 (inadmissível)" % (j + 1) for j in range(2) if g[j] < -tolA]
    ver += ["u%d < 0" % (j + 1) for j in range(2) if u[j] < -tolA]
    ver = "; ".join(ver) + "  x" if ver else "ponto KKT  V"
    print("   %-5d %s (%7.4f,%7.4f) (%7.4f,%7.4f) (%6.3f,%6.3f) %s"
          % (k + 1, hip[k], *x, *(u + 0.0), *g, ver))
    c += [perto(x, Xs[k], tolA), perto(u, Us[k], tolA)]
    if k == 2:
        x3, u3 = x, u
cK, r = kkt_check(df1(x3), g1(x3), Jg1, u3, tol=tolA)
gx3 = g1(x3)
print("   Sobrevive o caso 3: x* = (5/3, 1/3), f* = %.4f; ativas I = {%s}; u1 g1 = %.1e, u2 g2 = %.1e"
      % (f1(x3), ", ".join(map(str, r.ativas)), *(u3 * gx3)))
print("   Geometria (2): grad f(x*) = (%.4f, %.4f) = u1 grad g1 = (2/3)(-1,-4)" % tuple(df1(x3)))
c = [all(c), cK, perto(f1(x3), 1, tolA), perto(gx3[1], 4/3, tolA),
     perto(df1(x3), [-2/3, -8/3], tolA)]
print("   tabela dos 4 casos: %s | caso 3 é KKT: %s | f* = 1: %s | g2 = 4/3: %s"
      " | grad f = (2/3)(-1,-4): %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 3. Ativa com u = 0
print("\n3. min x^2  s.a.  g = x >= 0:  x* = 0, g(x*) = 0 (ativa), 2x - u = 0 => u* = 0")
cK, r = kkt_check([0.0], [0.0], [[1.0]], [0.0], tol=tolA)
print("   KKT: %s; ativas I = {%s} com u = 0 (ativa não implica u > 0)"
      % (sn(cK), ", ".join(map(str, r.ativas))))
ok = ok and cK and r.ativas == [1]

# ------------------------------------------------ 4. Exemplo 2
f2 = lambda x: (x[0] - 3)**2 + (x[1] - 2)**2
df2 = lambda x: np.array([2 * (x[0] - 3), 2 * (x[1] - 2)])
g2 = lambda x: np.array([5 - x[0]**2 - x[1]**2, 4 - x[0] - 2 * x[1], x[0], x[1]])
Jg2 = lambda x: np.array([[-2 * x[0], -2 * x[1]], [-1, -2], [1, 0], [0, 1]], float)
print("\n4. Exemplo 2: min (x1-3)^2 + (x2-2)^2  s.a.  5 - x1^2 - x2^2 >= 0, 4 - x1 - 2x2 >= 0,"
      " x1 >= 0, x2 >= 0")
x = np.array([2.0, 1.0]); g = g2(x); J = Jg2(x)
I = np.flatnonzero(np.abs(g) <= tolA)                       # ativas
u = np.zeros(4)
u[I] = np.linalg.solve(J[I].T, df2(x))                       # grad f = sum_{j em I} u_j grad g_j
cK, r = kkt_check(df2(x), g, J, u, tol=tolA)
rk = np.linalg.matrix_rank(J[I])
print("   x = (2, 1): g = (%s), ativas I = {%s}; grad f = (%g, %g), grad g1 = (%g, %g), grad g2 = (%g, %g)"
      % (" ".join("%g" % t for t in g), ", ".join(str(j + 1) for j in I), *df2(x), *J[0], *J[1]))
print("   u = (%s)  (slide: 1/3, 2/3, 0, 0);  LICQ: característica de [grad g1 grad g2] = %d;"
      "  f = %g;  KKT: %s" % (vec(u), rk, f2(x), sn(cK)))
c = [perto(u, [1/3, 2/3, 0, 0], tolA), rk == 2, perto(f2(x), 2, tolA), cK]
print("   u1 = 1/3, u2 = 2/3: %s | LICQ: %s | f = 2: %s | ponto KKT: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 5. KKT dá candidatos
print("\n5. min -x^2  s.a.  g1 = x + 1 >= 0,  g2 = 2 - x >= 0  (os 4 casos; «ambas ativas» é impossível)")
fm = lambda x: -x**2
gm = lambda x: np.array([x + 1, 2 - x])
Jgm = np.array([[1.0], [-1.0]])
c = []
for x, u, nome in [(0.0, (0, 0), "máximo"), (-1.0, (2, 0), "mínimo local"),
                   (2.0, (0, 4), "mínimo global")]:
    cK, r = kkt_check([-2 * x], gm(x), Jgm, u, tol=tolA)
    print("   x = %2g: u = (%g, %g), f = %2g, f'' = -2;  KKT: %s -> %s" % (x, *u, fm(x) + 0.0, sn(cK), nome))
    c.append(cK)
cand = [0.0, -1.0, 2.0]
c.append(cand[int(np.argmin([fm(t) for t in cand]))] == 2)
print("   os três são pontos KKT: %s | o menor f (global) é x = 2: %s" % (sn(all(c[:3])), sn(c[3])))
ok = ok and all(c)

# ------------------------------------------------ 6. O que medem os multiplicadores
print("\n6. Sensibilidade no exemplo 1: g1 = b - x1 - 4x2 >= 0, b = 3 -> 3.1 (caso 3: g1 ativa, u2 = 0)")
fb = lambda b: f1(caso_ex1(1, 0, b)[0])
f31 = fb(3.1)
print("   f*(3) = %.4f, f*(3.1) = %.5f (exato: (6 - b)^2/9); aproximação 1 - u1*0.1 = %.4f"
      % (fb(3.0), f31, 1 - 2 / 3 * 0.1))
c = [perto(fb(3.0), 1, tolA), confere(f31, 0.934, 3), confere(1 - 2 / 3 * 0.1, 0.933, 3),
     perto(f31, (6 - 3.1)**2 / 9, tolA)]
print("   f*(3) = 1: %s | f*(3.1) = 0.934: %s | aproximação 0.933: %s | (6-b)^2/9: %s" % simnao(c))
ok = ok and all(c)

# produção do 1.1: max 40x1 + 30x2  <=>  min -(40x1 + 30x2)
print("\n   Produção do 1.1: max 40x1 + 30x2 s.a. 2x1 + x2 <= 100 (máquina), x1 + x2 <= 80 (mão de obra),")
print("   x1 <= 40 (procura), x >= 0.  Em g_j >= 0: g = b - A x >= 0 e x >= 0;  f = -(40x1 + 30x2)")
Ap = np.array([[2.0, 1], [1, 1], [1, 0]]); bp = np.array([100.0, 80, 40]); cp = np.array([40.0, 30])
x = np.array([20.0, 60.0])
gp = np.concatenate([bp - Ap @ x, x]); Jgp = np.vstack([-Ap, np.eye(2)])
I = np.flatnonzero(np.abs(gp) <= tolA)
u = np.zeros(5)
u[I] = np.linalg.solve(Jgp[I].T, -cp)                        # grad f = -c = sum u_j grad g_j
cK, r = kkt_check(-cp, gp, Jgp, u, tol=tolA)
print("   x* = (20, 60), L* = %g;  ativas I = {%s} (máquina, mão de obra): 2u1 + u2 = 40, u1 + u2 = 30"
      % (cp @ x, ", ".join(str(j + 1) for j in I)))
print("   u = (%s)  -> preços-sombra 10 EUR/h, 20 EUR/h, 0;  KKT: %s" % (vec(u), sn(cK)))
c = [perto(cp @ x, 2600, tolA), perto(u, [10, 20, 0, 0, 0], tolA), cK]
print("   L* = 2600: %s | u = (10, 20, 0): %s | ponto KKT: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 7. Quando KKT falha: LICQ
print("\n7. min -x1  s.a.  g1 = (1 - x1)^3 - x2 >= 0,  g2 = x2 >= 0;  x* = (1, 0)")
x = np.array([1.0, 0.0]); gf = np.array([-1.0, 0.0])
J = np.array([[-3 * (1 - x[0])**2 + 0.0, -1], [0, 1]])
u = np.linalg.lstsq(J.T, gf, rcond=None)[0]                  # melhor u (mínimos quadrados)
resmin = np.linalg.norm(gf - J.T @ u)
print("   grad g1 = (%g, %g), grad g2 = (%g, %g): característica %d (dependentes)"
      % (*J[0], *J[1], np.linalg.matrix_rank(J)))
print("   min_u |grad f - J'u| = %.4f: a componente x1 de grad f = (-1, 0) não se obtém -> não há u" % resmin)
c = [np.linalg.matrix_rank(J) == 1, perto(resmin, 1, tolA)]
print("   gradientes dependentes: %s | sem multiplicadores (resíduo 1): %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 8. Exercícios
print("\n8. Exercícios")
# Ex. 1: min -ln(x1+1) - x2 s.a. 3 - 2x1 - x2 >= 0, x1 >= 0, x2 >= 0
fe1 = lambda x: -np.log(x[0] + 1) - x[1]
dfe1 = lambda x: np.array([-1 / (x[0] + 1), -1.0])
ge1 = lambda x: np.array([3 - 2 * x[0] - x[1], x[0], x[1]])
Jge1 = np.array([[-2.0, -1], [1, 0], [0, 1]])
x = np.array([0.0, 3.0]); g = ge1(x); I = np.flatnonzero(np.abs(g) <= tolA)
u = np.zeros(3); u[I] = np.linalg.solve(Jge1[I].T, dfe1(x))
cK1, r = kkt_check(dfe1(x), g, Jge1, u, tol=tolA)
print("   Ex. 1: x* = (0, 3), ativas I = {%s}, u = (%s), f* = %.4f;  KKT: %s"
      % (", ".join(str(j + 1) for j in I), vec(u), fe1(x), sn(cK1)))
c = [list(I + 1) == [1, 2], cK1]
# Ex. 2: min x1^2 + x2^2 s.a. 5 - x1^2 - x2^2 >= 0, x1 + 2x2 = 4, x >= 0
fe2 = lambda x: x[0]**2 + x[1]**2
dfe2 = lambda x: 2 * np.asarray(x, float)
ge2 = lambda x: np.array([5 - x[0]**2 - x[1]**2, x[0], x[1]])
Jge2 = lambda x: np.array([[-2 * x[0], -2 * x[1]], [1, 0], [0, 1]], float)
he2 = lambda x: np.array([x[0] + 2 * x[1] - 4])
Jhe2 = np.array([[1.0, 2.0]])
x = np.array([0.8, 1.6]); g = ge2(x)
lam = np.linalg.lstsq(Jhe2.T, dfe2(x), rcond=None)[0]        # só a igualdade é ativa: u = 0
cK2, r = kkt_check(dfe2(x), g, Jge2(x), np.zeros(3), he2(x), Jhe2, lam, tol=tolA)
print("   Ex. 2: x* = (4/5, 8/5), g = (%s) > 0 (inativas, u = 0), lambda = %.4f (slide: 8/5);  KKT: %s"
      % (vec(g), lam[0], sn(cK2)))
c += [r.ativas == [], perto(lam, 8 / 5, tolA), cK2]
# Ex. 4: a lata com h <= 6 (g = 6 - h >= 0), z = (r, h)
V = 330.0
A = lambda z: 2 * np.pi * z[0]**2 + 2 * np.pi * z[0] * z[1]
dA = lambda z: np.array([4 * np.pi * z[0] + 2 * np.pi * z[1], 2 * np.pi * z[0]])
hv = lambda z: np.array([np.pi * z[0]**2 * z[1] - V])
Jhv = lambda z: np.array([[2 * np.pi * z[0] * z[1], np.pi * z[0]**2]])
gh = lambda z: np.array([6 - z[1]])
Jgh = np.array([[0.0, -1.0]])
hh = 6.0; rr = np.sqrt(V / (6 * np.pi)); z = np.array([rr, hh])
lam = (2 * rr + hh) / (rr * hh)                              # de dL/dr = 0
u = lam * np.pi * rr**2 - 2 * np.pi * rr                     # de dL/dh = 0: 2 pi r - lambda pi r^2 + u = 0
cK4, r = kkt_check(dA(z), gh(z), Jgh, [u], hv(z), Jhv(z), [lam], tol=1e-9)
print("   Ex. 4: h* = 6, r* = %.4f cm, A* = %.4f cm^2, lambda = %.4f, u* = %.4f cm^2/cm;  KKT: %s"
      % (rr, A(z), lam, u, sn(cK4)))
c += [confere(rr, 4.18, 2), confere(A(z), 267.7, 1), confere(u, 5.19, 2), cK4]
print("   Ex. 1 ativas {1,2}: %s | KKT: %s | Ex. 2 inativas: %s | lambda = 8/5: %s | KKT: %s"
      " | Ex. 4 r* = 4.18: %s | A* = 267.7: %s | u* = 5.19: %s | KKT: %s" % simnao(c))
ok = ok and all(c)

# =================================================== (b) com um solver
print("\n=============== (b) com um solver ===============")
print("Solver: scipy.optimize.minimize(method='trust-constr'), multiplicadores em r.v.")
print("O SciPy escreve L = f + v'c: com NonlinearConstraint(g, 0, np.inf) e NonlinearConstraint(h, 0, 0),")
print("u_UC = -r.v[0] e lambda_UC = -r.v[1]  (fmincon: u_UC = lambda.ineqnonlin, lambda_UC = -lambda.eqnonlin).")
opts = dict(gtol=1e-10, xtol=1e-12, maxiter=3000)
# f, grad f, g (>= 0), jac g, h (= 0) ou None, jac h, x0, x* exato, u* exato, lambda* exato, nome
probs = [(f1, df1, g1, lambda x: Jg1, None, None, [0, 0], [5/3, 1/3], [2/3, 0], None, "Exemplo 1"),
         (f2, df2, g2, Jg2, None, None, [0.5, 0.5], [2, 1], [1/3, 2/3, 0, 0], None, "Exemplo 2"),
         (fe2, dfe2, ge2, Jge2, he2, lambda x: Jhe2, [1, 1], [0.8, 1.6], [0, 0, 0], [1.6],
          "Exercício 2"),
         (A, dA, gh, lambda z: Jgh, hv, Jhv, [3, 6], None, None, None, "Exercício 4 (lata, h <= 6)")]
for fo, dfo, gc, Jg, hc, Jh, x0, xs, us, ls, nome in probs:
    cons = [NonlinearConstraint(gc, 0, np.inf, jac=Jg)]
    if hc is not None:
        cons.append(NonlinearConstraint(hc, 0, 0, jac=Jh))
    r = minimize(fo, x0, jac=dfo, method="trust-constr", constraints=cons, options=opts)
    u = -r.v[0]
    lam = -r.v[1] if hc is not None else None
    print("   %s: x = (%s), f = %.4f, u = -r.v[0] = (%s)%s  [status %d, %d it., %d aval. de f]"
          % (nome, vec(r.x), r.fun, vec(u),
             "" if lam is None else ", lambda = -r.v[1] = (%s)" % vec(lam), r.status, r.nit, r.nfev))
    if xs is None:      # exercício 4: comparar com os valores do slide
        c = [confere(r.x[1], 6, 2), confere(r.x[0], 4.18, 2), confere(r.fun, 267.7, 1),
             confere(u, 5.19, 2)]
        print("      h* = 6: %s | r* = 4.18: %s | A* = 267.7: %s | u* = 5.19: %s" % simnao(c))
    else:
        c = [perto(r.x, xs, tolB), perto(u, us, tolB)]
        if ls is not None:
            c.append(perto(lam, ls, tolB))
        print("      x*: %s | u*: %s" % simnao(c[:2]) + ("" if ls is None else " | lambda*: %s" % sn(c[2])))
    ok = ok and all(c)

# produção (PL) com linprog, como no 1.1
res = linprog(c=-cp, A_ub=Ap, b_ub=bp, bounds=[(0, None)] * 2)
u = -res.ineqlin.marginals
print("   Produção com linprog (min -L): x = (%s), L = %.4f, res.ineqlin.marginals = (%s)"
      % (vec(res.x), -res.fun, vec(res.ineqlin.marginals)))
print("      preços-sombra = -marginals = (%s)" % vec(u))
c = [perto(res.x, [20, 60], tolB), perto(-res.fun, 2600, tolB), perto(u, [10, 20, 0], tolB)]
print("      x* = (20, 60): %s | L* = 2600: %s | preços-sombra (10, 20, 0): %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
