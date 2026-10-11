"""Reproduz os exemplos do deck 2.3 (interpolação quadrática sucessiva).

Exemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x, x1 = 1, Delta = 0.5,
         tolx = tolf = 1e-3: 4 iterações (k = 0..3), sem reflexões,
         n_f = 3 + 4 = 7, x* ~ 1.6094 (exato 1.6093811).
«Quando a parábola não serve»: f(x) = x^4/4 - x^2/2 com os trios
         (0.3, 0.45, 0.6) -> a2 < 0, xb = -0.46;
         (0.55, 0.6, 0.65) -> a2 ~ 0+, xb = 5.31, f(xb) = 184;
         e a salvaguarda (reflexão) a partir do primeiro trio.

Imprime a tabela das iterações como nos slides, os valores finais e as
contagens, e no fim compara com os slides (vírgula decimal nos slides,
ponto aqui).  Correr (de qualquer pasta):  python ex02_3_powell.py

Otimização — deck 2.3.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from powell_quadratic import powell_quadratic


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 2.3: interpolação quadrática sucessiva (Powell)")

# ------------------------------------------------------------ exemplo
f = lambda x: 8 * np.exp(1 - x) + 7 * np.log(x)
xs = 1.6093810677230782                     # minimizante (Newton com muitas casas)
x1, Delta, tolx, tolf, kmax = 1, 0.5, 1e-3, 1e-3, 50
print("\nExemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x, x1 = %g, Delta = %g, tolx = tolf = %g"
      % (x1, Delta, tolx))
r = powell_quadratic(f, x1, Delta, tolx, tolf, kmax, verbose=True)
print("Fim: x* = %.4f, f(x*) = %.4f; %d iterações, %d reflexões, n_f = 3 + %d = %d"
      % (r.x, r.fx, r.nit, r.nrefl, r.nit + r.nrefl, r.nfev))
print("Erros |xb - x*|:" + "".join(" %.1e" % e for e in np.abs(r.history[:, 9] - xs)))

# tabela do slide: x1 x2 x3 (4 casas) f1 f2 f3 (4) a1 a2 (3) xb f(xb) (4)
T = np.array([[1.0000, 1.5000, 2.0000, 8.0000, 7.6905, 7.7951, -0.619, 0.828, 1.6237, 7.6806],
              [1.5000, 1.6237, 2.0000, 7.6905, 7.6806, 7.7951, -0.080, 0.768, 1.6139, 7.6805],
              [1.5000, 1.6139, 1.6237, 7.6905, 7.6805, 7.6806, -0.088, 0.837, 1.6096, 7.6804],
              [1.5000, 1.6096, 1.6139, 7.6905, 7.6804, 7.6805, -0.092, 0.840, 1.6094, 7.6804]])
H = r.history
tabok = (H.shape[0] == 4
         and confere(H[:, [1, 2, 3, 4, 5, 6, 9, 10]], T[:, [0, 1, 2, 3, 4, 5, 8, 9]], 4)
         and confere(H[:, 7:9], T[:, 6:8], 3))
erros = np.abs(H[:, 9] - xs)
c = [tabok, r.nit == 4, r.nrefl == 0, r.nfev == 7,
     confere(r.x, 1.6094, 4), confere(r.fx, 7.6804, 4),
     confere(erros, [0.014, 0.0045, 0.0002, 5e-5], [3, 4, 4, 5]), bool(np.all(H[1:, 1] == 1.5))]
print("  tabela k = 0..3: %s | 4 it.: %s | sem reflexões: %s | n_f = 7: %s | x* = 1.6094: %s | f = 7.6804: %s | erros: %s | 1.5 fica no trio: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------ quando a parábola não serve
g = lambda x: x**4 / 4 - x**2 / 2
print("\nQuando a parábola não serve: f(x) = x^4/4 - x^2/2 (mínimos em -1 e 1)")
P = [[0.3, 0.45, 0.6], [0.55, 0.6, 0.65]]
xbP = []; a2P = []
for X in P:
    F = [g(xi) for xi in X]
    a1 = (F[1] - F[0]) / (X[1] - X[0])
    a2 = ((F[2] - F[0]) / (X[2] - X[0]) - a1) / (X[2] - X[1])
    xb = (X[0] + X[1]) / 2 - a1 / (2 * a2)
    xbP.append(xb); a2P.append(a2)
    print("  trio (%.2f, %.2f, %.2f): a1 = %.4f, a2 = %.4f, xb = %.2f, f(xb) = %.1f"
          % (*X, a1, a2, xb, g(xb)))
rr = powell_quadratic(g, 0.3, 0.15, 1e-3, 1e-3, 50)
print("  com a salvaguarda, de x1 = 0.3, Delta = 0.15: x* = %.4f, f = %.4f, %d reflexão(ões), %d it., n_f = %d"
      % (rr.x, rr.fx, rr.nrefl, rr.nit, rr.nfev))
c = [a2P[0] < 0, confere(xbP[0], -0.46, 2), a2P[1] > 0, confere(xbP[1], 5.31, 2),
     confere(g(xbP[1]), 184, 0), rr.flag == 0, confere(rr.x, 1, 2)]
print("  a2 < 0: %s | xb = -0.46: %s | a2 > 0: %s | xb = 5.31: %s | f(xb) = 184: %s | reflexão converge: %s | x* ~ 1: %s"
      % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
