"""Reproduz os exemplos do deck 2.2 (secção áurea).

Exemplo 2 (exemplo-guia): f(x) = x^2 + 54/x em [1,4], 15 reduções.
Exemplo 1:                f(x) = -3 sin(x) + x^2/8 em [0,4], 15 reduções.
Slide «Quantas reduções?»: L0 = 1, tolx = 0.01  =>  N = 10, n_f = 12.

Imprime as tabelas das iterações (k = 0..15; os slides mostram k = 0..7),
os valores finais e as contagens, e no fim compara com os slides.
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex02_2_golden_section.py

Otimização — deck 2.2.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from golden_section import golden_section


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 2.2: método da secção áurea")

# ------------------------------------------------------------ Exemplo 2
f = lambda x: x**2 + 54 / x
a, b, N = 1, 4, 15
print("\nExemplo 2: f(x) = x^2 + 54/x em [%g,%g], N = %d reduções" % (a, b, N))
r = golden_section(f, a, b, N, verbose=True)
L15 = r.history[-1, 7]
print("Após %d reduções: L_%d = %.4f, x* ~ %.4f, f(x*) ~ %.4f  (exato: x* = 3, f = 27)"
      % (N, N, L15, r.x, r.fx))
print("n_f = N + 2 = %d para localizar x*; com f no ponto médio final, n_f = N + 3 = %d"
      % (r.nfev_loc, r.nfev))

# tabela do slide (k = 0..7): a  b  x1  x2  f(x1)  f(x2)  L_k
T2 = [[1.0000, 4.0000, 2.1459, 2.8541, 29.7692, 27.0660, 3.0000],
      [2.1459, 4.0000, 2.8541, 3.2918, 27.0660, 27.2403, 1.8541],
      [2.1459, 3.2918, 2.5836, 2.8541, 27.5761, 27.0660, 1.1459],
      [2.5836, 3.2918, 2.8541, 3.0213, 27.0660, 27.0014, 0.7082],
      [2.8541, 3.2918, 3.0213, 3.1246, 27.0014, 27.0453, 0.4377],
      [2.8541, 3.1246, 2.9574, 3.0213, 27.0055, 27.0014, 0.2705],
      [2.9574, 3.1246, 3.0213, 3.0608, 27.0014, 27.0109, 0.1672],
      [2.9574, 3.0608, 2.9969, 3.0213, 27.0000, 27.0014, 0.1033]]
c = [confere(r.history[:8, 1:8], T2, 4), confere(L15, 0.0022, 4),
     confere(r.x, 3.0002, 4), confere(r.fx, 27.0000, 4),
     r.nfev_loc == 17, r.nfev == 18]
print("  tabela k = 0..7: %s | L_15: %s | x*: %s | f(x*): %s | n_f = 17: %s | n_f = 18: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ Exemplo 1
f = lambda x: -3 * np.sin(x) + x**2 / 8
a, b, N = 0, 4, 15
print("\nExemplo 1: f(x) = -3 sin(x) + x^2/8 em [%g,%g], N = %d reduções" % (a, b, N))
r = golden_section(f, a, b, N, verbose=True)
L15 = r.history[-1, 7]
print("Após %d reduções: L_%d = %.4f, x* ~ %.4f, f(x*) ~ %.4f;  n_f = %d (+1 = %d)"
      % (N, N, L15, r.x, r.fx, r.nfev_loc, r.nfev))
T1 = [[0.0000, 4.0000, 1.5279, 2.4721, -2.7054, -1.0977, 4.0000],
      [0.0000, 2.4721, 0.9443, 1.5279, -2.3188, -2.7054, 2.4721],
      [0.9443, 2.4721, 1.5279, 1.8885, -2.7054, -2.4040, 1.5279],
      [0.9443, 1.8885, 1.3050, 1.5279, -2.6818, -2.7054, 0.9443],
      [1.3050, 1.8885, 1.5279, 1.6656, -2.7054, -2.6397, 0.5836],
      [1.3050, 1.6656, 1.4427, 1.5279, -2.7152, -2.7054, 0.3607],
      [1.3050, 1.5279, 1.3901, 1.4427, -2.7096, -2.7152, 0.2229],
      [1.3901, 1.5279, 1.4427, 1.4752, -2.7152, -2.7143, 0.1378]]
c = [confere(r.history[:8, 1:8], T1, 4), confere(L15, 0.0029, 4),
     confere(r.x, 1.4489, 4), confere(r.fx, -2.7153, 4)]
print("  tabela k = 0..7: %s | L_15: %s | x*: %s | f(x*): %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ «Quantas reduções?»
tau = (np.sqrt(5) - 1) / 2
L0, tolx = 1, 0.01
Nmin = np.log(tolx / L0) / np.log(tau)
N = int(np.ceil(Nmin))
r = golden_section(lambda x: (x - 0.3)**2, 0, L0, N)
print("\nQuantas reduções? L0 = %g, tolx = %g: N >= %.2f => N = %d, tau^N = %.4f, n_f = %d"
      % (L0, tolx, Nmin, N, tau**N, r.nfev_loc))
c = [confere(Nmin, 9.57, 2), N == 10, confere(tau**N, 0.0081, 4), r.nfev_loc == 12]
print("  N >= 9.57: %s | N = 10: %s | 0.618^10: %s | n_f = 12: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
