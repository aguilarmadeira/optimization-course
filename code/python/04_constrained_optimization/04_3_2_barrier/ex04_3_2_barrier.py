"""Reproduz os exemplos do deck 4.3.2 (barreira logarítmica).

1. «Ideia»: min x^2 s.a. x - 1 >= 0; x_R = (1 + sqrt(1+2R))/2: 1.366, 1.048, 1.005
   (R = 1, 0.1, 0.01).
2. «O exemplo nas contas» e «As KKT reaparecem --- perturbadas»: min x1^2 + x2^2
   s.a. x1 + x2 - 1 >= 0, R = 1, 0.1, 0.01, 0.001: x_R, g(x_R), ||x_R - x*|| = O(R), u_R.
3. Barreira inversa R/g, R = 0.01: x_R = (0.548; 0.548), contra (0.5050; 0.5050) da log.
4. Restrição inativa: min (x1-1)^2 + (x2-1)^2 s.a. 9 - x1^2 - x2^2 >= 0:
   R = 1 -> (0.882; 0.882), R = 0.1 -> (0.986; 0.986); u_R = R/g -> 0.
5. A chamada do slide «No computador» (tolcomp = 1e-8, tolx = 1e-6, tmax = 12).
6. A comparação exterior vs. interior (slide «Exterior vs. interior»): barreira de
   (1,1) em 5 ciclos, 743 chamadas a P, 26 fora de X, n_f = 717; exterior de (0,0)
   em 6 ciclos, n_f = 844, e de (1,1), n_f = 864 (até ||x^(t) - x*|| <= 1e-4, sem a
   avaliação final). Usa penalty_exterior da pasta 04_3_1_exterior_penalty.

Minimizador interno: o Nelder–Mead da UC na versão do script da comparação
(nelder_mead_comp, tolx = 1e-10, tolf = 1e-12, kmax = 2000; ver o cabeçalho
de nelder_mead_comp.py e o README). As tabelas são uma única corrida com
arranque a quente (c = 0.1), lida em history.
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex04_3_2_barrier.py

Otimização — deck 4.3.2.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from barrier_log import barrier_log
from nelder_mead_comp import nelder_mead_comp
from penalty_exterior import penalty_exterior  # (4.3.1)


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


NM = lambda P, x: nelder_mead_comp(P, x, 1e-10, 1e-12, 2000)   # interno dos exemplos
ok = True
print("Otimização — deck 4.3.2: barreira logarítmica")

# ------------------------------------------------------------ 1. Ideia (1D)
r = barrier_log(lambda x: x[0]**2, lambda x: x[0] - 1, [2.0], 1, 0.1, 0, 0, 3, NM)
xR = r.history[1:, 5]
print("\n1. min x^2 s.a. x - 1 >= 0, de x = 2: x_R = (1 + sqrt(1+2R))/2")
for Ri, xi in zip(r.history[1:, 1], xR):
    print("   R = %5g: x_R = %.4f  (exato %.4f)" % (Ri, xi, (1 + np.sqrt(1 + 2 * Ri)) / 2))
c = [confere(xR, [1.366, 1.048, 1.005], 3)]
print("  x_R = 1.366, 1.048, 1.005: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 2. O exemplo nas contas
f = lambda x: x[0]**2 + x[1]**2
g = lambda x: x[0] + x[1] - 1
xs = np.array([0.5, 0.5])
r = barrier_log(f, g, [1, 1], 1, 0.1, 0, 0, 4, NM)
H = r.history[1:]
Rv = H[:, 1]; X = H[:, 5:7]; u = H[:, 7]
gv = X.sum(axis=1) - 1; err = np.linalg.norm(X - xs, axis=1)
print("\n2. min x1^2 + x2^2 s.a. x1 + x2 - 1 >= 0, de (1,1), R^(0) = 1, c = 0.1 (x* = (1/2,1/2), u* = 1)")
print("  %6s %20s %10s %12s %12s %8s %6s %6s"
      % ("R", "x_R", "g(x_R)", "||x_R-x*||", "||x_R-x*||/R", "u_R", "P", "fora"))
for i in range(4):
    print("  %6g   (%.4f; %.4f) %10.4f %12.5f %12.4f %8.4f %6d %6d"
          % (Rv[i], X[i, 0], X[i, 1], gv[i], err[i], err[i] / Rv[i], u[i], H[i, 2], H[i, 3]))
print("  erro O(R): ||x_R - x*||/R -> 1/sqrt 2 = %.4f;  u_R = R/g = (1 + sqrt(1+4R))/2 -> 1;"
      "  u_R g(x_R) = R" % (1 / np.sqrt(2)))
c = [confere(X, np.repeat([0.809, 0.5458, 0.5050, 0.5005], 2).reshape(4, 2),
             np.repeat([3, 4, 4, 4], 2).reshape(4, 2)),
     confere(gv, [0.618, 0.0916, 0.0099, 0.0010], [3, 4, 4, 4]),
     confere(err, [0.437, 0.0648, 0.0070, 0.00071], [3, 4, 4, 5]),
     confere(u, [1.618, 1.092, 1.010, 1.001], 3),
     all(gv > 0)]
print("  x_R: %s | g(x_R): %s | ||x_R - x*||: %s | u_R: %s | todos admissíveis: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 3. Barreira inversa
R = 0.01
Pinv = lambda x: np.inf if g(x) <= 0 else f(x) + R / g(x)
ri = nelder_mead_comp(Pinv, [1, 1], 1e-10, 1e-12, 2000)
print("\n3. Barreira inversa P = f + R/g, R = 0.01, de (1,1): x_R = (%.4f; %.4f), ||x_R - x*|| = %.4f"
      % (ri.x[0], ri.x[1], np.linalg.norm(ri.x - xs)))
print("   log (tabela acima), R = 0.01: x_R = (%.4f; %.4f), ||x_R - x*|| = %.4f  -> inversa O(sqrt R), log O(R)"
      % (X[2, 0], X[2, 1], err[2]))
c = [confere(ri.x, [0.548, 0.548], 3), confere(X[2], [0.505, 0.505], 3)]
print("  inversa (0.548; 0.548): %s | log (0.505; 0.505): %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 4. Restrição inativa
r = barrier_log(lambda x: (x[0] - 1)**2 + (x[1] - 1)**2, lambda x: 9 - x[0]**2 - x[1]**2,
                [0, 0], 1, 0.1, 0, 0, 2, NM)
H = r.history[1:]
print("\n4. Restrição inativa: min (x1-1)^2 + (x2-1)^2 s.a. 9 - x1^2 - x2^2 >= 0, de (0,0) (x* = (1,1), u* = 0)")
for row in H:
    print("   R = %4g: x_R = (%.4f; %.4f), g = %.4f, u_R = R/g = %.4f"
          % (row[1], row[5], row[6], 9 - row[5]**2 - row[6]**2, row[7]))
c = [confere(H[:, 5:7], [[0.882, 0.882], [0.986, 0.986]], 3), bool(H[1, 7] < H[0, 7] < 1)]
print("  x_R = (0.882; 0.882), (0.986; 0.986): %s | u_R a descer para 0: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 5. Chamada do slide
print("\n5. Chamada do slide: barrier_log(f, g, [1, 1], 1, 0.1, 1e-8, 1e-6, 12)")
print("   (interno por omissão: nelder_mead_comp; o slide usa o fminsearch do MATLAB)")
r = barrier_log(f, g, [1, 1], 1, 0.1, 1e-8, 1e-6, 12, verbose=True)
nfc = r.history[1:, 2]
print("   nit = %d (R = %.3g), u = %.4f, chamadas a P por ciclo %d-%d, nfev = chamadas a P + 1 = %d"
      % (r.nit, r.R, r.u[0], nfc.min(), nfc.max(), r.nfev))
print("   fora de X: %d; n_f efetivo = %d (nfev é majorante de n_f); %s"
      % (r.nfora, r.nfev_efetivo, r.message))
print("   nota: R no ciclo 9 = 0.1^8 calculado como %.17g > 1e-8 = tolcomp: J*R <= tolcomp"
      % r.history[9, 1])
print("         só pode valer a partir do ciclo 10 (R ~ 1e-9), qualquer que seja o minimizador interno")
print("   slide (fminsearch do MATLAB): nit = 10 (R = 1e-9), u ~ 1; nfev ~ 140-190 por ciclo")
c = [abs(r.u[0] - 1) <= 0.01]
print("  u ~ 1 (a menos de 1 %%): %s | nit = 10: %s (informativo: depende do minimizador interno)"
      % (simnao(c)[0], "sim" if r.nit == 10 else "não"))
ok = ok and all(c)

# ------------------------------------------------------------ 6. Comparação
print("\n6. Comparação exterior vs. interior: o mesmo Nelder–Mead interno (tolerâncias 1e-10),")
print("   arranque a quente; para no 1.º ciclo com ||x^(t) - x*|| <= 1e-4 (lido em history)")


def ate(H, tol=1e-4):
    """1.º ciclo t (linha de history) com ||x^(t) - x*|| <= tol."""
    e = np.linalg.norm(H[1:, 5:7] - xs, axis=1)
    return int(np.argmax(e <= tol)) + 1


rb = barrier_log(f, g, [1, 1], 1, 0.1, 0, 0, 12, NM)
tb = ate(rb.history); Hb = rb.history
nPb = int(Hb[tb, 4]); forab = int(Hb[1:tb + 1, 3].sum())
print("   barreira de (1,1), R^(0) = 1, c = 0.1:")
for row in Hb[1:tb + 1]:
    print("     t = %d  R = %8.1e  x = (%.5f, %.5f)  ||x - x*|| = %.1e  P = %3d (fora %2d)  acum. = %d"
          % (row[0], row[1], row[5], row[6], np.linalg.norm(row[5:7] - xs), row[2], row[3], row[4]))
print("     -> %d ciclos, %d chamadas a P, %d fora de X, n_f = %d (sem a avaliação final)"
      % (tb, nPb, forab, nPb - forab))
ext = {}
for x0 in ((0.0, 0.0), (1.0, 1.0)):
    re = penalty_exterior(f, g, lambda x: [], list(x0), 0.1, 10, 0, 0, 12, NM)
    te = ate(re.history)
    ext[x0] = (te, int(re.history[te, 3]), int(np.median(re.history[1:te + 1, 2])))
print("   exterior de (0,0), R^(0) = 0.1, c = 10: %d ciclos, n_f = %d" % ext[(0.0, 0.0)][:2])
print("   exterior de (1,1), R^(0) = 0.1, c = 10: %d ciclos, n_f = %d" % ext[(1.0, 1.0)][:2])
medb = int(np.median(Hb[1:tb + 1, 2]))
print("\n                                                 exterior  barreira")
print("   ciclos exteriores t                          %9d %9d" % (ext[(0.0, 0.0)][0], tb))
print("   por ciclo (mediana; n_f | chamadas a P)      %9d %9d" % (ext[(0.0, 0.0)][2], medb))
print("   total até ||x^(t) - x*|| <= 1e-4             %9d %9d" % (ext[(0.0, 0.0)][1], nPb))
c = [(tb, nPb, forab, nPb - forab) == (5, 743, 26, 717),
     ext[(0.0, 0.0)][:2] == (6, 844), ext[(1.0, 1.0)][:2] == (6, 864),
     (ext[(0.0, 0.0)][2], medb) == (140, 146)]
print("  barreira 5 ciclos, 743 P, 26 fora, n_f = 717: %s | exterior (0,0) 6 ciclos, 844: %s |"
      " exterior (1,1) 6 ciclos, 864: %s | medianas 140 e 146: %s" % simnao(c))
ok = ok and all(c)

# informativo: o mesmo com o nelder_mead do deck 3.2.3 (contração exterior aceite com <=)
from nelder_mead import nelder_mead  # noqa: E402  (deck 3.2.3)
rb3 = barrier_log(f, g, [1, 1], 1, 0.1, 0, 0, 12, lambda P, x: nelder_mead(P, x, 1e-10, 1e-12, 2000))
t3 = ate(rb3.history); H3 = rb3.history
print("   informativo — com nelder_mead do deck 3.2.3 (aceita a contração exterior com <=):")
print("   barreira de (1,1): %d ciclos, %d chamadas a P, %d fora de X, n_f = %d"
      % (t3, H3[t3, 4], H3[1:t3 + 1, 3].sum(), H3[t3, 4] - H3[1:t3 + 1, 3].sum()))

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
