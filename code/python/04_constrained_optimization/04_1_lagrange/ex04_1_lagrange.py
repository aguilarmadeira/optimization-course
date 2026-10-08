"""Reproduz os exemplos do deck 4.1 (multiplicadores de Lagrange).

(a) Resolve os exemplos dos slides pelas condições de Lagrange, com as contas
    do slide, e verifica numericamente as condições
    grad f(x*) + sum_l lambda_l grad h_l(x*) = 0  e  h_l(x*) = 0.
    Convenção das aulas:  L(x, lambda) = f(x) + sum_l lambda_l h_l(x).
    (Com L = f - sum lambda_l h_l, como em Nocedal e Wright, só muda o sinal de lambda.)
      1. O exemplo resolvido: min 2x1^2 + x2^2 s.a. x1 + x2 = 1
         (e «Geometria (1)» e «O que significa lambda?»).
      2. Candidatos: min x1 + x2 s.a. x1^2 + x2^2 = 2, classificados também
         com a Hessiana orlada (n = 2, uma restrição: det < 0 mínimo, det > 0 máximo).
      3. Exemplo 2, a lata: min 2 pi r^2 + 2 pi r h s.a. pi r^2 h = V.
      4. Várias restrições: min x1^2+x2^2+x3^2 s.a. x1+x2+x3 = 3, x1-x2 = 1.
      5. Quando a condição de Lagrange falha (gradientes dependentes).
(b) Resolve 1, 3 e 4 com scipy.optimize.minimize(method='trust-constr'),
    como no deck 4.2 («No computador»); os multiplicadores vêm em r.v.
    Sinal: o SciPy escreve L = f + v'c, com c = h e limites 0 <= c <= 0, logo
    lambda = r.v, sem troca de sinal (como no fmincon: lambda = lambda.eqnonlin).

Os slides usam vírgula decimal; aqui usa-se o ponto.
Não há função da UC neste deck: é um exemplo. Implementação didática.
Correr (de qualquer pasta):  python ex04_1_lagrange.py  (numpy e scipy)

Otimização — deck 4.1.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import warnings

import numpy as np
from scipy.optimize import NonlinearConstraint, minimize

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


def orlada(gh, HL):
    """Hessiana orlada (n = 2, uma restrição): [[0, gh'], [gh, HL]]."""
    gh = np.ravel(gh)
    return np.block([[np.zeros((1, 1)), gh[None, :]], [gh[:, None], np.atleast_2d(HL)]])


def lagrange_lin(Q, c, A, b):
    """min 1/2 x'Qx + c'x s.a. Ax = b: grad_x L = Qx + c + A'lambda = 0, Ax = b."""
    Q = np.atleast_2d(Q); A = np.atleast_2d(A); m = A.shape[0]
    M = np.block([[Q, A.T], [A, np.zeros((m, m))]])
    z = np.linalg.solve(M, np.concatenate([-np.ravel(c), np.ravel(b)]))
    return z[:Q.shape[0]], z[Q.shape[0]:]


fmt = lambda v: ", ".join("%.4f" % t for t in np.ravel(v))
ok = True
tolA = 1e-10        # contas exatas (a)
tolB = 1e-5         # solver (b)

print("Otimização — deck 4.1: multiplicadores de Lagrange")
print("Convenção das aulas: L = f + sum lambda_l h_l")
print("\n=============== (a) pelas condições de Lagrange ===============")

# ------------------------------------------------ 1. O exemplo resolvido
f1 = lambda x: 2 * x[0]**2 + x[1]**2
df1 = lambda x: np.array([4 * x[0], 2 * x[1]])
h1 = lambda x: x[0] + x[1] - 1
dh1 = np.array([1.0, 1.0])
print("\n1. min 2x1^2 + x2^2  s.a.  h = x1 + x2 - 1 = 0")
xg = np.array([0.8, 0.2]); d = np.array([1.0, -1.0])
print("   Geometria (1): em x = (0.8, 0.2), grad f = (%.1f, %.1f), grad f'd = %.1f (d = (1,-1))"
      % (*df1(xg), df1(xg) @ d))
c = [confere(df1(xg), [3.2, 0.4], 1), confere(df1(xg) @ d, 2.8, 1)]
# 4x1 + lambda = 0, 2x2 + lambda = 0, x1 + x2 = 1
x, lam = lagrange_lin(np.diag([4.0, 2.0]), [0, 0], [[1, 1]], [1])
res = np.linalg.norm(df1(x) + lam[0] * dh1, np.inf)
print("   Sistema: 4x1 + lambda = 0, 2x2 + lambda = 0, x1 + x2 = 1")
print("   x* = (%.4f, %.4f), lambda* = %.4f, f* = %.4f  (slide: (1/3, 2/3), -4/3, 2/3)"
      % (*x, lam[0], f1(x)))
print("   grad f(x*) = (%.4f, %.4f) = -lambda*(1,1);  |grad_x L| = %.1e, |h| = %.1e"
      % (*df1(x), res, abs(h1(x))))
detH1 = np.linalg.det(orlada(dh1, np.diag([4.0, 2.0])))
print("   Hessiana orlada [[0,1,1],[1,4,0],[1,0,2]]: det = %.4f < 0 -> mínimo" % detH1)
# O que significa lambda?  x1 + x2 = 1 + delta  ->  df*/d(delta) = -lambda*
fd = lambda delta: f1(lagrange_lin(np.diag([4.0, 2.0]), [0, 0], [[1, 1]], [1 + delta])[0])
dl = 1e-4
decl = (fd(dl) - fd(-dl)) / (2 * dl)
print("   Sensibilidade: f*(delta) = (2/3)(1+delta)^2, declive em 0 (dif. centrais) = %.6f = -lambda*"
      % decl)
c += [perto(x, [1/3, 2/3], tolA), perto(lam, -4/3, tolA), perto(f1(x), 2/3, tolA),
      res <= tolA, abs(h1(x)) <= tolA, perto(decl, 4/3, 1e-8), perto(detH1, -6, 1e-9)]
print("   Geometria (1): %s %s | x*: %s | lambda*: %s | f*: %s | grad_x L = 0: %s | h = 0: %s"
      " | declive = -lambda*: %s | orlada -6: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 2. Candidatos
f2 = lambda x: x[0] + x[1]
df2 = lambda x: np.array([1.0, 1.0])
h2 = lambda x: x[0]**2 + x[1]**2 - 2
dh2 = lambda x: np.array([2 * x[0], 2 * x[1]])
print("\n2. min x1 + x2  s.a.  h = x1^2 + x2^2 - 2 = 0")
print("   (1,1) + lambda (2x1, 2x2) = 0  =>  x1 = x2 = -1/(2 lambda);  na restrição, x1 = x2 = +-1")
dT = np.array([1.0, -1.0]) / np.sqrt(2)          # direção tangente (a mesma nos dois)
c = []
for x, slam, sf, sdet, nome in [((-1.0, -1.0), 0.5, -2, -8, "mínimo"), ((1.0, 1.0), -0.5, 2, 8, "máximo")]:
    x = np.array(x)
    lam = -1 / (2 * x[0])
    res = np.linalg.norm(df2(x) + lam * dh2(x), np.inf)
    H = 2 * lam * np.eye(2)                      # Hessiana de L em x
    q = dT @ H @ dT                              # 2.ª ordem em T
    dH = np.linalg.det(orlada(dh2(x), H))        # Hessiana orlada: det = -d'H d, d = (h_x2, -h_x1)
    print("   x = (%2d, %2d): lambda = %5.2f, f = %2d, |grad_x L| = %.1e, |h| = %.1e, "
          "d'H_L d = %5.2f, det orlada = %5.2f -> %s" % (*x, lam, f2(x), res, abs(h2(x)), q, dH, nome))
    c += [perto(lam, slam, tolA), perto(f2(x), sf, tolA), res <= tolA, abs(h2(x)) <= tolA,
          np.sign(q) == np.sign(lam), perto(dH, sdet, 1e-9), np.sign(dH) == -np.sign(q)]
print("   lambda, f, grad_x L = 0, h = 0 e classificação (2.ª ordem e Hessiana orlada) nos dois candidatos: %s"
      % simnao([all(c)]))
ok = ok and all(c)

# ------------------------------------------------ 3. A lata
V = 330.0
A3 = lambda z: 2 * np.pi * z[0]**2 + 2 * np.pi * z[0] * z[1]      # z = (r, h), h = altura
dA3 = lambda z: np.array([4 * np.pi * z[0] + 2 * np.pi * z[1], 2 * np.pi * z[0]])
h3 = lambda z: np.pi * z[0]**2 * z[1] - V
dh3 = lambda z: np.array([2 * np.pi * z[0] * z[1], np.pi * z[0]**2])
print("\n3. Lata: min A = 2 pi r^2 + 2 pi r h  s.a.  pi r^2 h = V = %g cm^3" % V)
print("   dL/dh = 2 pi r + lambda pi r^2 = 0 => lambda = -2/r;  dL/dr = 0 => h = 2r;  pi r^2 (2r) = V")
r = (V / (2 * np.pi))**(1 / 3); h = 2 * r; lam = -2 / r; z = np.array([r, h])
res = np.linalg.norm(dA3(z) + lam * dh3(z), np.inf)
HL = np.pi * np.array([[4 + 2 * lam * h, 2 + 2 * lam * r],
                       [2 + 2 * lam * r, 0.0]])  # Hessiana de L = pi [-4 -2; -2 0]
dT = np.array([1.0, -4.0])                       # dh3(z)'dT = 0
print("   r* = %.4f cm, h* = %.4f cm, A* = %.4f cm^2, lambda* = -2/r* = %.4f cm^2/cm^3 (dA*/dV = -lambda*)"
      % (r, h, A3(z), lam))
print("   |grad_x L| = %.1e, |h| = %.1e;  2.ª ordem: d = (1,-4) em T, d'H_L d = %.4f = 12 pi > 0 (mínimo)"
      % (res, abs(h3(z)), dT @ HL @ dT))
c = [confere(r, 3.745, 3), confere(r, 3.74, 2), confere(h, 7.49, 2), confere(A3(z), 264.4, 1),
     confere(lam, -0.534, 3), res <= 1e-9 * np.linalg.norm(dA3(z)), abs(h3(z)) <= 1e-9 * V,
     perto(dT @ HL @ dT, 12 * np.pi, 1e-9)]
print("   r* = 3.745 (3.74): %s %s | h* = 7.49: %s | A* = 264.4: %s | lambda* = -0.534: %s"
      " | grad_x L = 0: %s | h = 0: %s | 12 pi: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 4. Várias restrições
f4 = lambda x: x[0]**2 + x[1]**2 + x[2]**2
df4 = lambda x: 2 * np.asarray(x)
A4 = np.array([[1.0, 1, 1], [1, -1, 0]]); b4 = np.array([3.0, 1])   # h = A4 x - b4
print("\n4. min x1^2 + x2^2 + x3^2  s.a.  x1 + x2 + x3 = 3,  x1 - x2 = 1")
x, lam = lagrange_lin(2 * np.eye(3), np.zeros(3), A4, b4)
res = np.linalg.norm(df4(x) + A4.T @ lam, np.inf)
d = np.array([1.0, 1, -2])
print("   x* = (%s), lambda* = (%s), f* = %.4f  (slide: (3/2, 1/2, 1), (-2, -1), 7/2)"
      % (fmt(x), fmt(lam), f4(x)))
print("   grad f(x*) = (%g, %g, %g) = 2 (1,1,1) + 1 (1,-1,0) = -lambda1 grad h1 - lambda2 grad h2;"
      "  grad f'd = %g (d = (1,1,-2));"
      "  |grad_x L| = %.1e, |h| = %.1e"
      % (*df4(x), df4(x) @ d + 0.0, res, np.linalg.norm(A4 @ x - b4, np.inf)))
c = [perto(x, [1.5, 0.5, 1], tolA), perto(lam, [-2, -1], tolA), perto(f4(x), 3.5, tolA),
     perto(df4(x) @ d, 0, tolA), res <= tolA, np.linalg.norm(A4 @ x - b4, np.inf) <= tolA]
print("   x*: %s | lambda*: %s | f*: %s | grad f'd = 0: %s | grad_x L = 0: %s | h = 0: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 5. Quando Lagrange falha
print("\n5. min x1 + x2 + x3^2  s.a.  x1 - 1 = 0,  x1^2 + x2^2 - 1 = 0;  x* = (1, 0, 0)")
x = np.array([1.0, 0, 0])
gf = np.array([1.0, 1, 2 * x[2]])
J = np.array([[1.0, 0, 0], [2 * x[0], 2 * x[1], 0]])   # linhas: grad h1', grad h2'
lam = np.linalg.lstsq(J.T, -gf, rcond=None)[0]         # melhor lambda (mínimos quadrados)
resmin = np.linalg.norm(gf + J.T @ lam)
print("   grad h1 = (%g,%g,%g), grad h2 = (%g,%g,%g): característica %d (dependentes)"
      % (*J[0], *J[1], np.linalg.matrix_rank(J)))
print("   min_lambda |grad f + J'lambda| = %.4f: dL/dx2 = 1 para todos os lambda -> sistema sem solução"
      % resmin)
c = [np.linalg.matrix_rank(J) == 1, perto(resmin, 1, tolA)]
print("   gradientes dependentes: %s | sem multiplicadores (resíduo 1): %s" % simnao(c))
ok = ok and all(c)

# =================================================== (b) com um solver
print("\n=============== (b) com um solver ===============")
print("Solver: scipy.optimize.minimize(method='trust-constr'), multiplicadores em r.v.")
print("O SciPy escreve L = f + v'c (aqui c = h, com 0 <= c <= 0), como nas aulas: lambda = r.v.")
opts = dict(gtol=1e-10, xtol=1e-12, maxiter=2000)
probs = [(f1, df1, lambda x: [h1(x)], lambda x: [dh1], [0, 0], [1/3, 2/3], [-4/3],
          "1. exemplo resolvido"),
         (A3, dA3, lambda z: [h3(z)], lambda z: [dh3(z)], [3, 6], None, None,
          "3. lata (x0 = (3, 6))"),
         (f4, df4, lambda x: A4 @ x - b4, lambda x: A4, [0, 0, 0], [1.5, 0.5, 1], [-2, -1],
          "4. várias restrições")]
for fo, dfo, hc, Jh, x0, xs, ls, nome in probs:
    r = minimize(fo, x0, jac=dfo, method="trust-constr", options=opts,
                 constraints=[NonlinearConstraint(hc, 0, 0, jac=Jh)])
    lam = r.v[0]
    print("   %s: x = (%s), f = %.4f, lambda = r.v = (%s)  [status %d, %d it., %d aval. de f]"
          % (nome, fmt(r.x), r.fun, fmt(lam), r.status, r.nit, r.nfev))
    if xs is None:      # a lata: comparar com os valores do slide
        c = [confere(r.x[0], 3.745, 3), confere(r.x[1], 7.49, 2), confere(r.fun, 264.4, 1),
             confere(lam, -0.534, 3)]
        print("      r*: %s | h*: %s | A*: %s | lambda*: %s" % simnao(c))
    else:
        c = [perto(r.x, xs, tolB), perto(lam, ls, tolB)]
        print("      x*: %s | lambda*: %s" % simnao(c))
    ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
