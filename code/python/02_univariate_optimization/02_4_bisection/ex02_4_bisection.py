"""Reproduz os exemplos do deck 2.4 (bisseção, paragem |f'(z)| <= eps).

Exemplo 2 (exemplo-guia): f(x) = x^2 + 54/x em [1,4], eps = 0.05:
          7 iterações, z = 2.9922, f'(z) = -0.0470, n_g = 7 + 2 = 9, n_f = 0.
          Variante com 15 reduções: x2 - x1 = 3/2^15 = 9.2e-5 (secção áurea: 0.0022).
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
f = lambda x: x**2 + 54 / x
df = lambda x: 2 * x - 54 / x**2
a, b, eps = 1, 4, 0.05
print("\nExemplo 2: f(x) = x^2 + 54/x em [%g,%g], f'(%g) = %g < 0, f'(%g) = %g > 0, eps = %g"
      % (a, b, a, df(a), b, df(b), eps))
r = bisection(f, df, a, b, eps, verbose=True)
print("|f'(z)| = %.4f <= eps na iteração %d: x* ~ z = %.4f, f(z) = %.4f (exato: x* = 3)"
      % (abs(r.history[-1, 5]), r.nit, r.x, r.fx))
print("n_g = %d + 2 = %d, n_f = %d" % (r.nit, r.ngev, r.nfev))
rv = bisection(f, df, a, b, 0, tolx=3 / 2**15)
tau = (np.sqrt(5) - 1) / 2; Lgs = 3 * tau**15
print("Variante, 15 reduções (tolx = 3/2^15): %d iterações, x2 - x1 = %.1e, x = %.4f, n_g = %d"
      % (rv.nit, rv.L, rv.x, rv.ngev))
print("Secção áurea, mesmo [1,4], 15 reduções: L = %.4f, %.0f vezes maior" % (Lgs, Lgs / rv.L))
T2 = [[1, 1.0000, 4.0000, 2.5000, 27.8500, -3.6400],
      [2, 2.5000, 4.0000, 3.2500, 27.1779, 1.3876],
      [3, 2.5000, 3.2500, 2.8750, 27.0482, -0.7831],
      [4, 2.8750, 3.2500, 3.0625, 27.0116, 0.3674],
      [5, 2.8750, 3.0625, 2.9688, 27.0030, -0.1895],
      [6, 2.9688, 3.0625, 3.0156, 27.0007, 0.0933],
      [7, 2.9688, 3.0156, 2.9922, 27.0002, -0.0470]]
c = [r.history.shape == (7, 6) and confere(r.history, T2, 4), r.flag == 0,
     confere(r.x, 2.9922, 4), r.ngev == 9, r.nfev == 0,
     rv.nit == 15 and rv.flag == 2, confere(rv.L, 9.2e-5, 6), rv.ngev == 17,
     confere(Lgs, 0.0022, 4), round(Lgs / rv.L) == 24]
print("  tabela k = 1..7: %s | parou por |f'| <= eps: %s | z = 2.9922: %s | n_g = 9: %s | n_f = 0: %s"
      " | variante 15 it.: %s | 9.2e-5: %s | n_g = 17: %s | L áurea: %s | 24 vezes: %s" % simnao(c))
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
