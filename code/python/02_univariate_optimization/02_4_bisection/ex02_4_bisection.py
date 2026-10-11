"""Reproduz os exemplos do deck 2.4 (bisseção, paragem |f'(z)| <= eps).

Exemplo 2 (exemplo-guia, 1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x em [1,2],
          eps = 0.01: 6 iterações, z = 1.6094 (1.609375), f'(z) = -0.00001,
          n_g = 6 + 2 = 8, n_f = 0; com eps = 0.03 pararia na it. 3 (z = 1.625).
          Variante com 15 reduções: x2 - x1 = 1/2^15 = 3.1e-5 (secção áurea: 0.0007).
Exemplo 1: f(x) = -2 sin(x) + x^2/16 em [-1,3], eps = 0.03 (e eps = 0.05).
«Para resolver», eps = 0.01: x + 31/(x+3) em [0,5] e -3 sin(x) + x^2/8 em [0,4].
«Quantas iterações?» (variante): [0,2], tolx = 1e-3 -> k = 11; tolx = 1e-4 -> k = 15.

Imprime as tabelas das iterações (k = 1, 2, ...), os valores finais e as
contagens, e no fim compara com os slides (vírgula decimal nos slides, ponto
aqui). A coluna f(z) é só para leitura: o método não a usa.
Correr (de qualquer pasta):  python ex02_4_bisection.py

Otimização — deck 2.4.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from bisection import bisection


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 2.4: bisseção (paragem |f'(z)| <= eps)")

# ------------------------------------------------------------ Exemplo 2
f = lambda x: 8 * np.exp(1 - x) + 7 * np.log(x)
df = lambda x: -8 * np.exp(1 - x) + 7 / x
a, b, eps = 1, 2, 0.01
print("\nExemplo 2 (exemplo-guia, 1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x em [%g,%g], f'(%g) = %g < 0, f'(%g) = %.4f > 0, eps = %g"
      % (a, b, a, df(a), b, df(b), eps))
r = bisection(f, df, a, b, eps, verbose=True)
print("|f'(z)| = %.5f <= eps na iteração %d: x* ~ z = %.4f, f(z) = %.5f (exato: x* = 1.6093811)"
      % (abs(r.history[-1, 5]), r.nit, r.x, r.fx))
print("n_g = %d + 2 = %d, n_f = %d" % (r.nit, r.ngev, r.nfev))
r3 = bisection(f, df, a, b, 0.03)
print("Com eps = 0.03: termina na iteração %d, z = %.4f" % (r3.nit, r3.x))
rv = bisection(f, df, a, b, 0, tolx=1 / 2**15)
tau = (np.sqrt(5) - 1) / 2; Lgs = tau**15
print("Variante, 15 reduções (tolx = 1/2^15): %d iterações, x2 - x1 = %.1e, x = %.4f, n_g = %d"
      % (rv.nit, rv.L, rv.x, rv.ngev))
print("Secção áurea, mesmo [1,2], 15 reduções: L = %.4f, %.0f vezes maior" % (Lgs, Lgs / rv.L))
T2 = [[1, 1.0000, 2.0000, 1.5000, 7.69050, -0.18558],
      [2, 1.5000, 2.0000, 1.7500, 7.69624, 0.22107],
      [3, 1.5000, 1.7500, 1.6250, 7.68065, 0.02560],
      [4, 1.5000, 1.6250, 1.5625, 7.68227, -0.07826],
      [5, 1.5625, 1.6250, 1.5938, 7.68065, -0.02586],
      [6, 1.5938, 1.6250, 1.6094, 7.68045, -0.00001]]
D2 = [[0, 4, 4, 4, 5, 5]] * 6
c = [r.history.shape == (6, 6) and confere(r.history, T2, D2), r.flag == 0,
     confere(r.x, 1.6094, 4), r.ngev == 8, r.nfev == 0, r3.nit == 3 and r3.x == 1.625,
     rv.nit == 15 and rv.flag == 2, confere(rv.L, 3.1e-5, 6), rv.ngev == 17,
     confere(Lgs, 0.0007, 4), round(Lgs / rv.L) == 24]
print("  tabela k = 1..6: %s | parou por |f'| <= eps: %s | z = 1.6094: %s | n_g = 8: %s | n_f = 0: %s | eps = 0.03: 3 it.: %s"
      " | variante 15 it.: %s | 3.1e-5: %s | n_g = 17: %s | L áurea: %s | 24 vezes: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ Exemplo 1
f = lambda x: -2 * np.sin(x) + x**2 / 16
df = lambda x: -2 * np.cos(x) + x / 8
a, b, eps = -1, 3, 0.03
print("\nExemplo 1: f(x) = -2 sin(x) + x^2/16 em [%g,%g], f'(%g) = %.2f < 0, f'(%g) = %.2f > 0, eps = %g"
      % (a, b, a, df(a), b, df(b), eps))
r = bisection(f, df, a, b, eps, verbose=True)
xs = bisection(None, df, a, b, 1e-12).x                 # referência («exato»)
print("|f'(z)| = %.4f <= eps na iteração %d: x* ~ z = %.4f (exato: %.4f, erro %.4f)"
      % (abs(r.history[-1, 5]), r.nit, r.x, xs, abs(r.x - xs)))
r5 = bisection(f, df, a, b, 0.05)
print("Com eps = 0.05: termina na iteração %d, z = %.4f, |f'(z)| = %.4f, erro %.3f"
      % (r5.nit, r5.x, abs(r5.history[-1, 5]), abs(r5.x - xs)))
T1 = [[1, -1.0000, 3.0000, 1.0000, -1.6204, -0.9556],
      [2, 1.0000, 3.0000, 2.0000, -1.5686, 1.0823],
      [3, 1.0000, 2.0000, 1.5000, -1.8544, 0.0460],
      [4, 1.0000, 1.5000, 1.2500, -1.8003, -0.4744],
      [5, 1.2500, 1.5000, 1.3750, -1.8436, -0.2172],
      [6, 1.3750, 1.5000, 1.4375, -1.8531, -0.0861],
      [7, 1.4375, 1.5000, 1.4688, -1.8548, -0.0201]]
c = [r.history.shape == (7, 6) and confere(r.history, T1, 4), confere([df(-1), df(3)], [-1.21, 2.35], 2),
     confere(r.x, 1.4688, 4), confere(xs, 1.4783, 4),
     r5.nit == 3 and confere([r5.x, abs(r5.history[-1, 5]), abs(r5.x - xs)], [1.5, 0.0460, 0.022], [1, 4, 3])]
print("  tabela k = 1..7: %s | f'(-1), f'(3): %s | z = 1.4688: %s | exato 1.4783: %s | eps = 0.05: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------- Para resolver
print("\nPara resolver (eps = 0.01):")
f = lambda x: x + 31 / (x + 3); df = lambda x: 1 - 31 / (x + 3)**2
r3 = bisection(f, df, 0, 5, 0.01)
print("  x + 31/(x+3) em [0,5]:       %2d iterações, z = %.4f, f(z) = %.4f  (exato: sqrt(31) - 3 = %.4f; erro %.4f)"
      % (r3.nit, r3.x, r3.fx, np.sqrt(31) - 3, abs(r3.x - (np.sqrt(31) - 3))))
f = lambda x: -3 * np.sin(x) + x**2 / 8; df = lambda x: -3 * np.cos(x) + x / 4
r4 = bisection(f, df, 0, 4, 0.01)
x4 = bisection(None, df, 0, 4, 1e-12).x
r4v = bisection(None, df, 0, 4, 0, tolx=4 / 2**15)
print("  -3 sin(x) + x^2/8 em [0,4]:  %2d iterações, z = %.4f, f(z) = %.4f  (x* = %.4f; variante, 15 reduções: x2 - x1 = %.1e)"
      % (r4.nit, r4.x, r4.fx, x4, r4v.L))
c = [r3.nit == 6 and confere([r3.x, r3.fx], [2.5781, 8.1355], 4), confere(np.sqrt(31) - 3, 2.5678, 4),
     confere(abs(r3.x - (np.sqrt(31) - 3)), 0.0104, 4),
     r4.nit == 10 and confere([r4.x, r4.fx], [1.4492, -2.7153], 4), confere(x4, 1.4497, 4),
     r4v.nit == 15 and confere(r4v.L, 1.2e-4, 5)]
print("  1.º: %s | exato: %s | erro 0.0104: %s | 2.º: %s | x* = 1.4497: %s | 1.2e-4: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ Quantas iterações?
kk = [int(np.ceil(np.log2(2 / t))) for t in (1e-3, 1e-4)]
dq = lambda x: x - 0.3                                  # um f' qualquer com zero em [0,2]
nq = [bisection(None, dq, 0, 2, 0, tolx=t).nit for t in (1e-3, 1e-4)]
print("\nQuantas iterações? (variante) [0,2]: log2(2000) = %.2f -> k = %d (tolx = 1e-3); k = %d (tolx = 1e-4);"
      " a função faz %d e %d" % (np.log2(2000), kk[0], kk[1], nq[0], nq[1]))
c = [confere(np.log2(2000), 10.97, 2), kk[0] == 11, kk[1] == 15, nq == kk]
print("  log2 2000 = 10.97: %s | k = 11: %s | k = 15: %s | função: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
