"""Reproduz os exemplos do deck 2.4 (bisseção).

Exemplo 2 (exemplo-guia): f(x) = x^2 + 54/x em [1,4], 15 reduções:
          x ~ 3.0000, b - a = 3/2^15 = 9.2e-5, n_g = N + 2 = 17, n_f = 0.
Exemplo 1: f(x) = -2 sin(x) + x^2/16 em [-1,3], 15 reduções.
«Para resolver»: x + 31/(x+3) em [0,5] e -3 sin(x) + x^2/8 em [0,4].
«Quantas iterações?»: [0,2], tolx = 1e-3 -> k = 11; tolx = 1e-4 -> k = 15.

Imprime as tabelas das iterações (k = 0..14; os slides mostram k = 0..7),
os valores finais e as contagens, e no fim compara com os slides
(vírgula decimal nos slides, ponto aqui). A coluna f(x_k) é só para
leitura: o método não a usa.
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
print("Otimização — deck 2.4: bisseção")

# ------------------------------------------------------------ Exemplo 2
f = lambda x: x**2 + 54 / x
df = lambda x: 2 * x - 54 / x**2
a, b, N = 1, 4, 15
print("\nExemplo 2: f(x) = x^2 + 54/x em [%g,%g], f'(%g) = %g < 0, f'(%g) = %g > 0, N = %d"
      % (a, b, a, df(a), b, df(b), N))
r = bisection(f, df, a, b, N, verbose=True)
print("Após %d reduções, ponto médio: x* ~ %.4f, f(x*) ~ %.4f, b - a = %.1e (exato: x* = 3)"
      % (N, r.x, r.fx, r.L))
print("n_g = N + 2 = %d, n_f = %d" % (r.ngev, r.nfev))
tau = (np.sqrt(5) - 1) / 2; Lgs = 3 * tau**15
print("Secção áurea, mesmo [1,4], 15 reduções: L = %.4f, %.0f vezes maior" % (Lgs, Lgs / r.L))
T2 = [[1.0000, 4.0000, 2.5000, 27.8500, -3.6400],
      [2.5000, 4.0000, 3.2500, 27.1779, 1.3876],
      [2.5000, 3.2500, 2.8750, 27.0482, -0.7831],
      [2.8750, 3.2500, 3.0625, 27.0116, 0.3674],
      [2.8750, 3.0625, 2.9688, 27.0030, -0.1895],
      [2.9688, 3.0625, 3.0156, 27.0007, 0.0933],
      [2.9688, 3.0156, 2.9922, 27.0002, -0.0470],
      [2.9922, 3.0156, 3.0039, 27.0000, 0.0234]]
c = [confere(r.history[:8, 1:6], T2, 4), confere(r.x, 3.0000, 4),
     confere(r.L, 9.2e-5, 6), r.ngev == 17, r.nfev == 0,
     confere(Lgs, 0.0022, 4), round(Lgs / r.L) == 24]
print("  tabela k = 0..7: %s | x* ~ 3.0000: %s | b - a = 9.2e-5: %s | n_g = 17: %s | n_f = 0: %s | L áurea: %s | 24 vezes: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ Exemplo 1
f = lambda x: -2 * np.sin(x) + x**2 / 16
df = lambda x: -2 * np.cos(x) + x / 8
a, b, N = -1, 3, 15
print("\nExemplo 1: f(x) = -2 sin(x) + x^2/16 em [%g,%g], f'(%g) = %.2f < 0, f'(%g) = %.2f > 0, N = %d"
      % (a, b, a, df(a), b, df(b), N))
r = bisection(f, df, a, b, N, verbose=True)
print("Após %d reduções, ponto médio: x* ~ %.4f, f(x*) ~ %.4f, b - a = %.1e" % (N, r.x, r.fx, r.L))
T1 = [[-1.0000, 3.0000, 1.0000, -1.6204, -0.9556],
      [1.0000, 3.0000, 2.0000, -1.5686, 1.0823],
      [1.0000, 2.0000, 1.5000, -1.8544, 0.0460],
      [1.0000, 1.5000, 1.2500, -1.8003, -0.4744],
      [1.2500, 1.5000, 1.3750, -1.8436, -0.2172],
      [1.3750, 1.5000, 1.4375, -1.8531, -0.0861],
      [1.4375, 1.5000, 1.4688, -1.8548, -0.0201],
      [1.4688, 1.5000, 1.4844, -1.8548, 0.0129]]
c = [confere(r.history[:8, 1:6], T1, 4), confere([df(-1), df(3)], [-1.21, 2.35], 2),
     confere(r.x, 1.4783, 4), confere(r.fx, -1.8549, 4), confere(r.L, 1.2e-4, 5)]
print("  tabela k = 0..7: %s | f'(-1), f'(3): %s | x*: %s | f(x*): %s | b - a = 1.2e-4: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------- Para resolver
print("\nPara resolver (15 reduções, ponto médio final):")
f = lambda x: x + 31 / (x + 3); df = lambda x: 1 - 31 / (x + 3)**2
r3 = bisection(f, df, 0, 5, 15)
print("  x + 31/(x+3) em [0,5]:       x* = %.4f, f(x*) = %.4f  (exato: sqrt(31) - 3 = %.4f)"
      % (r3.x, r3.fx, np.sqrt(31) - 3))
f = lambda x: -3 * np.sin(x) + x**2 / 8; df = lambda x: -3 * np.cos(x) + x / 4
r4 = bisection(f, df, 0, 4, 15)
print("  -3 sin(x) + x^2/8 em [0,4]:  x* = %.4f, f(x*) = %.4f, b - a = %.1e" % (r4.x, r4.fx, r4.L))
c = [confere([r3.x, r3.fx], [2.5678, 8.1355], 4), confere(np.sqrt(31) - 3, 2.5678, 4),
     confere([r4.x, r4.fx], [1.4496, -2.7153], 4), confere(r4.L, 1.2e-4, 5)]
print("  1.º: %s | exato: %s | 2.º: %s | b - a = 1.2e-4: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ Quantas iterações?
kk = [int(np.ceil(np.log2(2 / t))) for t in (1e-3, 1e-4)]
print("\nQuantas iterações? [0,2]: log2(2000) = %.2f -> k = %d (tolx = 1e-3); k = %d (tolx = 1e-4)"
      % (np.log2(2000), kk[0], kk[1]))
c = [confere(np.log2(2000), 10.97, 2), kk[0] == 11, kk[1] == 15]
print("  log2 2000 = 10.97: %s | k = 11: %s | k = 15: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
