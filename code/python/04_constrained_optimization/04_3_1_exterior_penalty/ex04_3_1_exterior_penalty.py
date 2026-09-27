"""Reproduz os exemplos do deck 4.3.1 (penalização exterior).

1. «Ideia»: min x^2 s.a. x - 1 >= 0; x_R = R/(1+R): 1/2, 0.909 (R = 1, 10).
2. «O exemplo nas contas» e «As KKT reaparecem»: min x1^2 + x2^2 s.a.
   x1 + x2 - 1 >= 0, R = 1, 10, 100, 1000: x_R, g(x_R), ||x_R - x*|| = O(1/R), u_R.
3. Restrição inativa: min (x1-1)^2 + (x2-1)^2 s.a. 9 - x1^2 - x2^2 >= 0:
   x_R = (1,1) e u_R = 0 para todo o R.
4. «A escala das restrições importa»: g1 = 1 - x1, g2 = 100(1 - x2), R = 1, 10, 100.
5. «E com igualdades?»: min 2x1^2 + x2^2 s.a. x1 + x2 - 1 = 0: x_R, h(x_R), lambda_R.
6. A chamada do slide «No computador» (tolviol = tolx = 1e-6, tmax = 12).
7. A comparação do deck 4.3.2 (lado exterior): de (0,0), 6 ciclos, n_f = 844;
   de (1,1), 6 ciclos, n_f = 864 (até ||x^(t) - x*|| <= 1e-4, sem a avaliação final).

Minimizador interno: o Nelder–Mead da UC na versão do script da comparação
(nelder_mead_comp, tolx = 1e-10, tolf = 1e-12, kmax = 2000; ver o cabeçalho
de nelder_mead_comp.py e o README). As tabelas 2-5 são uma única corrida
com arranque a quente (R^(0) = 1, c = 10), lida em history.
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex04_3_1_exterior_penalty.py

Otimização — deck 4.3.1.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from nelder_mead_comp import nelder_mead_comp
from penalty_exterior import penalty_exterior



def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


NM = lambda P, x: nelder_mead_comp(P, x, 1e-10, 1e-12, 2000)   # interno dos exemplos
sem = lambda x: []                                              # sem restrições deste tipo
ok = True
print("Otimização — deck 4.3.1: penalização exterior")

# ------------------------------------------------------------ 1. Ideia (1D)
r = penalty_exterior(lambda x: x[0]**2, lambda x: x[0] - 1, sem, [0.0], 1, 10, 0, 0, 2, NM)
xR = r.history[1:, 5]
print("\n1. min x^2 s.a. x - 1 >= 0: x_R = R/(1+R) (R = 0: x_R = 0, o mínimo livre)")
for Ri, xi in zip(r.history[1:, 1], xR):
    print("   R = %4g: x_R = %.4f  (exato %.4f)" % (Ri, xi, Ri / (1 + Ri)))
c = [confere(xR, [0.5, 0.909], 3)]
print("  x_R = 1/2, 0.909: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 2. O exemplo nas contas
f = lambda x: x[0]**2 + x[1]**2
g = lambda x: x[0] + x[1] - 1
xs = np.array([0.5, 0.5])
r = penalty_exterior(f, g, sem, [0, 0], 1, 10, 0, 0, 4, NM)
H = r.history[1:]
Rv = H[:, 1]; X = H[:, 5:7]; u = H[:, 7]
gv = X.sum(axis=1) - 1; err = np.linalg.norm(X - xs, axis=1)
print("\n2. min x1^2 + x2^2 s.a. x1 + x2 - 1 >= 0, de (0,0), R^(0) = 1, c = 10 (x* = (1/2,1/2), u* = 1)")
print("  %6s %20s %10s %12s %12s %8s %8s" % ("R", "x_R", "g(x_R)", "||x_R-x*||", "R*||x_R-x*||", "u_R", "n_f"))
for i in range(4):
    print("  %6g   (%.4f; %.4f) %10.4f %12.5f %12.4f %8.4f %8d"
          % (Rv[i], X[i, 0], X[i, 1], gv[i], err[i], Rv[i] * err[i], u[i], H[i, 2]))
print("  erro O(1/R): R*||x_R - x*|| -> 1/(2 sqrt 2) = %.4f;  u_R = 2R/(1+2R) -> 1" % (1 / (2 * np.sqrt(2))))
c = [confere(X, np.repeat([0.333, 0.476, 0.4975, 0.4998], 2).reshape(4, 2),
             np.repeat([3, 3, 4, 4], 2).reshape(4, 2)),
     confere(gv, [-0.333, -0.0476, -0.0050, -0.0005], [3, 4, 4, 4]),
     confere(err, [0.236, 0.0337, 0.0035, 0.00035], [3, 4, 4, 5]),
     confere(u, [0.667, 0.952, 0.995, 0.9995], [3, 3, 3, 4])]
print("  x_R: %s | g(x_R): %s | ||x_R - x*||: %s | u_R: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 3. Restrição inativa
r = penalty_exterior(lambda x: (x[0] - 1)**2 + (x[1] - 1)**2, lambda x: 9 - x[0]**2 - x[1]**2,
                     sem, [0, 0], 1, 10, 0, 0, 4, NM)
H = r.history[1:]
print("\n3. Restrição inativa: min (x1-1)^2 + (x2-1)^2 s.a. 9 - x1^2 - x2^2 >= 0")
for row in H:
    print("   R = %4g: x_R = (%.4f; %.4f), g = %.4f, u_R = %.4f"
          % (row[1], row[5], row[6], 9 - row[5]**2 - row[6]**2, row[7] + 0.0))
print("   (%s: x não mudou e a violação é 0)" % r.message)
c = [confere(H[:, 5:7], np.ones_like(H[:, 5:7]), 4), confere(H[:, 7], 0 * H[:, 7], 4)]
print("  x_R = (1,1): %s | u_R = 0: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 4. Escala das restrições
fe = lambda x: (x[0] - 2)**2 + (x[1] - 2)**2
r = penalty_exterior(fe, lambda x: [1 - x[0], 100 * (1 - x[1])], sem, [0, 0], 1, 10, 0, 0, 3, NM)
H = r.history[1:]
kap = (2 + 2e4 * H[:, 1]) / (2 + 2 * H[:, 1])     # Hessiana de P em x_R: diag(2+2R, 2+2e4 R)
print("\n4. Escala: min (x1-2)^2 + (x2-2)^2 s.a. 1 - x1 >= 0 e 100(1 - x2) >= 0 (x* = (1,1))")
print("  %6s %10s %12s %10s" % ("R", "x1_R", "x2_R", "kappa"))
for i in range(3):
    print("  %6g %10.4f %12.7f %10.0f" % (H[i, 1], H[i, 5], H[i, 6], kap[i]))
r2 = penalty_exterior(fe, lambda x: [1 - x[0], 1 - x[1]], sem, [0, 0], 1, 10, 0, 0, 1, NM)
print("  com g2 = 1 - x2 e R = 1: x_R = (%.4f; %.4f) = ((2+R)/(1+R), (2+R)/(1+R)), kappa = 1"
      % tuple(r2.x))
c = [confere(H[:, 5], [1.50, 1.091, 1.0099], [2, 3, 4]),
     confere(H[:, 6], [1.0001, 1.00001, 1.000001], [4, 5, 6]),
     confere(kap, [5000, 9091, 9901], 0), confere(r2.x, [1.5, 1.5], 4)]
print("  x1_R: %s | x2_R: %s | kappa: %s | com g2: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 5. Igualdades
r = penalty_exterior(lambda x: 2 * x[0]**2 + x[1]**2, sem, lambda x: x[0] + x[1] - 1,
                     [0, 0], 1, 10, 0, 0, 4, NM)
H = r.history[1:]
X = H[:, 5:7]; hv = X.sum(axis=1) - 1; lam = H[:, 7]
print("\n5. Igualdade: min 2x1^2 + x2^2 s.a. x1 + x2 - 1 = 0 (x* = (1/3, 2/3), lambda* = 4/3)")
for i in range(4):
    print("   R = %4g: x_R = (%.4f; %.4f), h = %.4f, lambda_R = -2Rh = %.3f"
          % (H[i, 1], X[i, 0], X[i, 1], hv[i], lam[i]))
c = [confere(X, [[0.2000, 0.4000], [0.3125, 0.6250], [0.3311, 0.6623], [0.3331, 0.6662]], 4),
     confere(hv, [-0.4000, -0.0625, -0.0066, -0.0007], 4),
     confere(lam, [0.800, 1.250, 1.325, 1.332], 3)]
print("  x_R: %s | h(x_R): %s | lambda_R: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 6. Chamada do slide
print("\n6. Chamada do slide: penalty_exterior(f, g, h, [0, 0], 0.1, 10, 1e-6, 1e-6, 12)")
print("   (interno por omissão: nelder_mead_comp; o slide usa o fminsearch do MATLAB)")
r = penalty_exterior(f, g, sem, [0, 0], 0.1, 10, 1e-6, 1e-6, 12, verbose=True)
nfc = r.history[1:, 2]
print("   nit = %d (R = %g), u = %.4f, n_f por ciclo %d-%d, nfev = %d (inclui f(x) final)"
      % (r.nit, r.R, r.u[0], nfc.min(), nfc.max(), r.nfev))
print("   %s" % r.message)
H = r.history
print("   nota: com x_R exatos o critério valeria no ciclo 9 (viol = 5e-8, passo = %.1e); aqui o passo"
      % (np.sqrt(2) * (1e7 / (1 + 2e7) - 1e6 / (1 + 2e6))))
print("         no ciclo 9 é %.1e > 1e-6 (erro do interno ao longo de g = 0, P mal condicionada)"
      % (np.linalg.norm(H[9, 5:7] - H[8, 5:7]) / max(1, np.linalg.norm(H[8, 5:7]))))
print("   slide (fminsearch do MATLAB): nit = 10 (R = 1e8); u ~ 1; nfev ~ 130-210 por ciclo")
c = [abs(r.u[0] - 1) <= 0.01]
print("  u ~ 1 (a menos de 1 %%): %s | nit = 10: %s (informativo: depende do minimizador interno)"
      % (simnao(c)[0], "sim" if r.nit == 10 else "não"))
ok = ok and all(c)

# ------------------------------------------------------------ 7. Comparação (deck 4.3.2)
print("\n7. Comparação exterior vs. interior (deck 4.3.2), lado exterior: R^(0) = 0.1, c = 10,")
print("   tolerâncias internas 1e-10; para no 1.º ciclo com ||x^(t) - x*|| <= 1e-4 (lido em history)")


def ciclos_ate(r, tol=1e-4):
    """1.º ciclo t com ||x^(t) - x*|| <= tol e n_f acumulado até aí (sem a avaliação final)."""
    H = r.history
    e = np.linalg.norm(H[1:, 5:7] - xs, axis=1)
    t = int(np.argmax(e <= tol)) + 1
    return t, int(H[t, 3]), H[1:t + 1, 2]


res = {}
for x0 in ([0.0, 0.0], [1.0, 1.0]):
    r = penalty_exterior(f, g, sem, x0, 0.1, 10, 0, 0, 12, NM)
    t, nf, nfc = ciclos_ate(r)
    res[tuple(x0)] = (t, nf, int(np.median(nfc)))
    print("   de (%g,%g):" % tuple(x0))
    for row in r.history[1:t + 1]:
        print("     t = %d  R = %8.1e  x = (%.5f, %.5f)  ||x - x*|| = %.1e  n_f = %3d  acum. = %d"
              % (row[0], row[1], row[5], row[6], np.linalg.norm(row[5:7] - xs), row[2], row[3]))
    print("     -> %d ciclos, n_f = %d, mediana por ciclo = %d" % res[tuple(x0)])
c = [res[(0.0, 0.0)][:2] == (6, 844), res[(0.0, 0.0)][2] == 140, res[(1.0, 1.0)][:2] == (6, 864)]
print("  (0,0): 6 ciclos, n_f = 844: %s | mediana 140 (~140 por ciclo): %s | (1,1): 6 ciclos, n_f = 864: %s"
      % simnao(c))
ok = ok and all(c)

# informativo: o mesmo com o nelder_mead do deck 3.2.3 (contração exterior aceite com <=)
from nelder_mead import nelder_mead  # noqa: E402  (deck 3.2.3)
NM3 = lambda P, x: nelder_mead(P, x, 1e-10, 1e-12, 2000)
v = [ciclos_ate(penalty_exterior(f, g, sem, x0, 0.1, 10, 0, 0, 12, NM3))[:2]
     for x0 in ([0, 0], [1, 1])]
print("   informativo — com nelder_mead do deck 3.2.3 (aceita a contração exterior com <=):")
print("   (0,0): %d ciclos, n_f = %d;  (1,1): %d ciclos, n_f = %d" % (*v[0], *v[1]))

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
