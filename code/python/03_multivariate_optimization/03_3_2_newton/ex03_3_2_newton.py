"""Reproduz os exemplos do deck 3.3.2 (método de Newton).

«Podemos usar a curvatura?»: f(x) = (x1^2 + 9 x2^2)/2, x0 = (9, 1):
    H d = -g dá d0 = (-9, -1) e x1 = (0, 0) -- uma iteração.
«Rosenbrock: muito rápido, mas não monótono»: Newton puro de (-1.5; 2),
    tabela k = 0..7 (f sobe para 751 e para 40), 7 iterações, f < 1e-16.
«Newton procura g = 0»: Himmelblau de (0, 0), lambda(H0) = (-42; -26);
    o puro vai em 4 iterações para o máximo local (-0.27; -0.92), f = 181.6;
    o amortecido deteta g^T d >= 0, usa -g e chega a (3, 2) em 5 iterações.
«Newton amortecido»: Rosenbrock, 14 iterações (monótono); os últimos
    minimizantes da linha são 0.91; 0.98; 1.00.
Deck 3.3.3 (quadro de consulta): amortecido com ||g|| < 1e-4: 14 it.,
    n_g = 15, n_H = 14.

Paragem ||g|| <= tolg = 1e-8 (a das figuras), salvo indicação.
Imprime as tabelas, os valores finais e as contagens, e no fim compara com
os slides (vírgula decimal nos slides, ponto aqui).
Correr (de qualquer pasta):  python ex03_3_2_newton.py

Otimização — deck 3.3.2.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from newton_nd import newton_nd, DIRECOES


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def confere_sig(v, s, m=3):
    """v confere com s mostrado com m algarismos significativos?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    e = np.floor(np.log10(np.abs(s)))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (e - m + 1) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 3.3.2: método de Newton")

# ------------------------------------------------ quadrática
A = np.diag([1.0, 9.0])
quad = lambda x: 0.5 * (x[0]**2 + 9 * x[1]**2)     # = x^T A x / 2
gq = lambda x: np.array([x[0], 9 * x[1]])
Hq = lambda x: A
x0 = np.array([9.0, 1.0])
d0 = -np.linalg.solve(A, gq(x0))
r = newton_nd(quad, gq, Hq, x0)
print("\nQuadrática (x1^2 + 9 x2^2)/2 de (9, 1): gradiente -g0 = (%g, %g); Newton H d = -g: d0 = (%g, %g)"
      % (*(-gq(x0)), *d0))
print("x1 = (%g, %g): %d iteração; n_g = %d, n_H = %d" % (*r.x, r.nit, r.ngev, r.nhev))
c = [np.allclose(d0, [-9, -1]), np.allclose(r.x, 0), r.nit == 1]
print("  d0 = (-9, -1): %s | x1 = (0, 0): %s | 1 it.: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ Rosenbrock, puro
rosen = lambda x: (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2
grosen = lambda x: np.array([-2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0]**2), 200 * (x[1] - x[0]**2)])
Hrosen = lambda x: np.array([[2 - 400 * (x[1] - 3 * x[0]**2), -400 * x[0]], [-400 * x[0], 200]])
xr = np.array([-1.5, 2.0])
print("\nRosenbrock de (-1.5; 2), Newton puro (tolg = 1e-8):")
r = newton_nd(rosen, grosen, Hrosen, xr, 1e-8, 50, verbose=True)
H = r.history
print("%s; n_g = %d, n_H = %d, n_f = %d (f só nas iteradas)" % (r.message, r.ngev, r.nhev, r.nfev))
subidas = [(int(H[k, 0]), H[k, 3]) for k in range(1, H.shape[0]) if H[k, 3] > H[k - 1, 3]]
print("f sobe em: " + ", ".join("k = %d (f = %.2f)" % s for s in subidas))
Tx = [[-1.5000, 2.0000], [-1.4510, 2.1029], [0.2044, -2.6986], [0.2059, 0.0424],
      [0.9997, 0.3692], [0.9997, 0.9993], [1.0000, 1.0000], [1.0000, 1.0000]]
Tf = [12.500000, 6.007882, 751.610005, 0.630622, 39.701728]
Tf56 = [1.09e-7, 1.20e-12]
Tg = [1.63e2, 6.31e0, 5.92e2, 1.59e0, 2.82e2, 6.61e-4, 4.89e-5]
c = [H.shape[0] == 8 and confere(H[:, 1:3], Tx, 4) and confere(H[:5, 3], Tf, 6)
     and confere_sig(H[5:7, 3], Tf56) and H[7, 3] == 0 and confere_sig(H[:7, 4], Tg) and H[7, 4] == 0,
     r.nit == 7, r.fx < 1e-16, [s[0] for s in subidas] == [2, 4]]
print("  tabela k = 0..7: %s | 7 it.: %s | f < 1e-16: %s | não monótono (sobe em k = 2 e 4): %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ Rosenbrock, amortecido
print("\nRosenbrock de (-1.5; 2), Newton amortecido (salvaguarda g^T d >= 0 -> -g; pesquisa em linha de Brent):")
r = newton_nd(rosen, grosen, Hrosen, xr, 1e-8, 50, damped=True, verbose=True)
H = r.history
mono = bool(np.all(np.diff(H[:, 3]) <= 0))
print("%s; n_g = %d, n_H = %d, n_f = %d (%d nas pesquisas em linha)" % (r.message, r.ngev, r.nhev, r.nfev, r.nfev_ls))
print("monótono: %s; últimos minimizantes da linha: %.2f; %.2f; %.2f" % ("sim" if mono else "não", *H[-3:, 5]))
c = [r.nit == 14, mono, confere(H[-3:, 5], [0.91, 0.98, 1.00], 2)]
print("  14 it.: %s | monótono: %s | alpha = 0.91; 0.98; 1.00: %s" % simnao(c))
ok = ok and all(c)
r4 = newton_nd(rosen, grosen, Hrosen, xr, 1e-4, 50, damped=True)
print("Com tolg = 1e-4 (quadro de consulta do 3.3.3): %d it., n_g = %d, n_H = %d" % (r4.nit, r4.ngev, r4.nhev))
c = [r4.nit == 14, r4.ngev == 15, r4.nhev == 14]
print("  14 it.: %s | n_g = 15: %s | n_H = 14: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ Rosenbrock, amortecido com Armijo
print("\nRosenbrock de (-1.5; 2), amortecido com recuo de Armijo (alpha = 1, 1/2, ...; c1 = 1e-4):")
ra = newton_nd(rosen, grosen, Hrosen, xr, 1e-8, 50, damped=True, ls="armijo")
Ha = ra.history
x1 = Ha[1, 1:3]; g1 = grosen(x1); d1 = -np.linalg.solve(Hrosen(x1), g1); f1 = rosen(x1)
print("k = 1: f(x1) = %.3f, g^T d = %.3f" % (f1, g1 @ d1))
recuo = [(a, rosen(x1 + a * d1), f1 + 1e-4 * a * (g1 @ d1)) for a in (1, 0.5, 0.25, 0.125)]
for a, fa, lim in recuo:
    print("  alpha = %-6g f = %10.3f  limite de Armijo = %.3f  %s" % (a, fa, lim, "aceite" if fa <= lim else "rejeitado"))
monoa = bool(np.all(np.diff(Ha[:, 3]) <= 0))
ult = int(np.argmax(Ha[1:, 5][::-1] != 1)) if np.any(Ha[1:, 5] != 1) else ra.nit
print("%s; n_g = %d, n_H = %d, n_f = %d; monótono: %s; alpha = 1 nas últimas %d iterações"
      % (ra.message, ra.ngev, ra.nhev, ra.nfev, "sim" if monoa else "não", ult))
c = [confere([r_[1] for r_ in recuo], [751.610, 49.736, 7.145, 5.238], 3), confere(f1, 6.008, 3),
     [r_[1] <= r_[2] for r_ in recuo] == [False, False, False, True],
     ra.nit == 22, ra.nfev == 29, monoa, ult == 7]
print("  751,6 / 49,7 / 7,1 / 5,24: %s | f(x1) = 6,008: %s | aceita 1/8: %s | 22 it.: %s | n_f = 29: %s"
      " | monótono: %s | últimas 7 com alpha = 1: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ Himmelblau
himmel = lambda x: (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2
ghimmel = lambda x: np.array([4 * x[0] * (x[0]**2 + x[1] - 11) + 2 * (x[0] + x[1]**2 - 7),
                              2 * (x[0]**2 + x[1] - 11) + 4 * x[1] * (x[0] + x[1]**2 - 7)])
Hhimmel = lambda x: np.array([[12 * x[0]**2 + 4 * x[1] - 42, 4 * x[0] + 4 * x[1]],
                              [4 * x[0] + 4 * x[1], 4 * x[0] + 12 * x[1]**2 - 26]])
xh = np.array([0.0, 0.0])
lam0 = np.linalg.eigvalsh(Hhimmel(xh))
print("\nHimmelblau de (0, 0): lambda(H0) = (%g; %g) -- definida negativa" % tuple(lam0))
rp = newton_nd(himmel, ghimmel, Hhimmel, xh, 1e-8, 50)
lamx = np.linalg.eigvalsh(Hhimmel(rp.x))
print("Puro: %d it. -> x = (%.4f; %.4f), f = %.4f, lambda(H) = (%.1f; %.1f) < 0: máximo local"
      % (rp.nit, *rp.x, rp.fx, *lamx))
rd = newton_nd(himmel, ghimmel, Hhimmel, xh, 1e-8, 50, damped=True)
print("Amortecido: direção na 1.ª iteração: %s; %d it. -> x = (%.4f; %.4f), f = %.1e"
      % (DIRECOES[int(rd.history[1, 7])], rd.nit, *rd.x, rd.fx))
c = [np.allclose(lam0, [-42, -26]), rp.nit == 4, confere(rp.x, [-0.27, -0.92], 2), confere(rp.fx, 181.6, 1),
     bool(np.all(lamx < 0)), rd.history[1, 7] == 2, rd.nit == 5, np.allclose(rd.x, [3, 2])]
print("  lambda(H0) = (-42; -26): %s | puro 4 it.: %s | (-0.27; -0.92): %s | f = 181.6: %s | máximo: %s"
      " | amortecido usa -g: %s | 5 it.: %s | (3, 2): %s" % simnao(c))
ok = ok and all(c)

ra = newton_nd(himmel, ghimmel, Hhimmel, xh, 1e-8, 50, damped=True, ls="armijo")
rm = newton_nd(himmel, ghimmel, Hhimmel, xh, 1e-8, 50, damped=True, modify=True, ls="armijo")
print("Com Armijo: -g na 1.ª iteração: %d it., n_f = %d; H + mu I: %d it., n_f = %d; ambos em (%.0f, %.0f)"
      % (ra.nit, ra.nfev, rm.nit, rm.nfev, *rm.x))
c = [ra.history[1, 7] == 2, ra.nit == 6, ra.nfev == 14, rm.history[1, 7] == 1, rm.nit == 8, rm.nfev == 10,
     np.allclose(ra.x, [3, 2]), np.allclose(rm.x, [3, 2])]
print("  -g: %s | 6 it.: %s | n_f = 14: %s | H + mu I: %s | 8 it.: %s | n_f = 10: %s | (3, 2): %s, %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
