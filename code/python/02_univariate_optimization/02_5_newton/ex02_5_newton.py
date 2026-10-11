"""Reproduz os exemplos do deck 2.5 (Newton-Raphson).

Exemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x, x0 = 1:
    x1 = 2, x2 = 1.5331535, ...; 5 iterações até x* = 1.6093811 (erro 1.3e-13);
    n_g = n_H = 1 por iteração (+ f'(x0)). Paragem |f'| < tolg com tolg = 1e-6
    (o slide não fixa tolg; qualquer valor entre 2.2e-13 e 1.08e-6 dá as
    mesmas 5 iterações). Na pergunta de abertura, x0 = 1.5: x1 = 1.607.
«Convergência quadrática»: x^3 - 2x^2 + 4, x0 = 1, tolg = 1e-3 (3 it.).
«Newton procura f' = 0»: maximizar 8 + 5x - 3x^4 - 2x^6, x0 = 0.5,
    tolg = 1e-7 (5 it.).
«Os quatro métodos»: o exemplo-guia a partir de x0 = 1 (erro após 5 it.).

Imprime as tabelas como nos slides, os valores finais e as contagens, e
no fim compara com os slides (vírgula decimal nos slides, ponto aqui).
Correr (de qualquer pasta):  python ex02_5_newton.py

Otimização — deck 2.5.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from newton_1d import newton_1d


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 2.5: Newton-Raphson")

# ------------------------------------------------------- exemplo-guia
f = lambda x: 8 * np.exp(1 - x) + 7 * np.log(x)
df = lambda x: -8 * np.exp(1 - x) + 7 / x
ddf = lambda x: 8 * np.exp(1 - x) - 7 / x**2
xs = 1.6093810677230782                     # minimizante (Newton com muitas casas)
x0, tolg, kmax = 1, 1e-6, 20
print("\nExemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x, x0 = %g, tolg = %g" % (x0, tolg))
r = newton_1d(f, df, ddf, x0, tolg, kmax)
H = r.history
print("%3s %11s %11s %9s %11s %13s" % ("k", "x_k", "df(x_k)", "ddf(x_k)", "x_k+1", "|x_k+1 - x*|"))
for row in H:
    print("%3d %11.7f %+11.5f %9.4f %11.7f %13.1e" % (row[0], row[1], row[3], row[4], row[5], abs(row[5] - xs)))
print("x* = %.7f, f(x*) = %.7f; %d iterações, n_g = %d, n_H = %d, n_f = %d"
      % (r.x, r.fx, r.nit, r.ngev, r.nhev, r.nfev))
x1q = 1.5 - df(1.5) / ddf(1.5)
print("Pergunta de abertura: em x0 = 1.5, f'(1.5) = %.4f, f''(1.5) = %.4f, x1 = %.3f" % (df(1.5), ddf(1.5), x1q))
Tg = np.array([[1.0000000, -1.00000, 1.0000, 2.0000000],
               [2.0000000, 0.55696, 1.1930, 1.5331535],
               [1.5331535, -0.12826, 1.7160, 1.6078977],
               [1.6078977, -0.00244, 1.6484, 1.6093804],
               [1.6093804, -0.0000011, 1.6469, 1.6093811]])
eg = np.abs(H[:, 5] - xs)
c = [H.shape[0] == 5 and confere(H[:, [1, 5]], Tg[:, [0, 3]], 7) and confere(H[:4, 3], Tg[:4, 1], 5)
     and confere(H[4, 3], -1.1e-6, 7) and confere(H[:, 4], Tg[:, 2], 4),
     confere(eg, [3.9e-1, 7.6e-2, 1.5e-3, 6.6e-7, 1.3e-13], [2, 3, 4, 8, 14]),
     confere([df(1.5), ddf(1.5), x1q], [-0.1856, 1.7411, 1.607], [4, 4, 3]),
     r.nit == 5, r.ngev == 6, r.nhev == 5, r.nfev == 0]
print("  tabela k = 0..4: %s | erros: %s | abertura x1 = 1.607: %s | 5 it.: %s | n_g = 6: %s | n_H = 5: %s | n_f = 0: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------- convergência quadrática
f = lambda x: x**3 - 2 * x**2 + 4
df = lambda x: 3 * x**2 - 4 * x
ddf = lambda x: 6 * x - 4
x0, tolg = 1, 1e-3
print("\nConvergência quadrática: f(x) = x^3 - 2x^2 + 4, x0 = %g, tolg = %g" % (x0, tolg))
r = newton_1d(f, df, ddf, x0, tolg, 20)
H = r.history
print("%3s %9s %10s %9s %9s %11s" % ("k", "x_k", "df(x_k)", "ddf(x_k)", "x_k+1", "df(x_k+1)"))
for row in H:
    print("%3d %9.5f %+10.5f %9.4f %9.5f %+11.5f" % tuple(row[[0, 1, 3, 4, 5, 6]]))
e = np.abs(np.append(H[:, 1], r.x) - 4 / 3)
Hx = ddf(r.x)                        # classificar: 1 avaliação extra de f''
print("|f'(x_3)| = %.1e < tolg: para em x_3 = %.5f (x* = 4/3); f''(x_3) = %.4f > 0 (f''(4/3) = %g): mínimo"
      % (abs(H[-1, 6]), r.x, Hx, ddf(4 / 3)))
print("Erros e_k: %.4f -> %.4f -> %.4f -> %.4f;  e_3/e_2^2 = %.2f (C = 0.75)" % (*e, e[3] / e[2]**2))
Tc = np.array([[1.00000, -1.00000, 2.0000, 1.50000, 0.75000],
               [1.50000, 0.75000, 5.0000, 1.35000, 0.06750],
               [1.35000, 0.06750, 4.1000, 1.33354, 0.00081]])
c = [H.shape[0] == 3 and confere(H[:, [1, 3, 5, 6]], Tc[:, [0, 1, 3, 4]], 5) and confere(H[:, 4], Tc[:, 2], 4),
     confere(abs(H[-1, 6]), 8.1e-4, 5), confere(e, [0.33, 0.17, 0.017, 0.0002], [2, 2, 3, 4]),
     ddf(4 / 3) == 4, Hx > 0, r.ngev == 4, r.nhev == 3]
print("  tabela k = 0..2: %s | |f'(x_3)| = 8.1e-4: %s | erros: %s | f''(4/3) = 4: %s | mínimo: %s | n_g = 4: %s | n_H = 3: %s"
      % simnao(c))
ok = ok and all(c)
r03 = newton_1d(f, df, ddf, 0.3, 1e-8, 20)
print("De x0 = 0.3: x = %.6f, f''(x) = %g < 0: máximo (ponto estacionário errado)" % (r03.x, ddf(r03.x)))
c = [abs(r03.x) < 1e-6, ddf(r03.x) < 0]
print("  converge para 0: %s | máximo: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------ Newton procura f' = 0 (maximizar)
f = lambda x: 8 + 5 * x - 3 * x**4 - 2 * x**6
df = lambda x: 5 - 12 * (x**3 + x**5)
ddf = lambda x: -12 * (3 * x**2 + 5 * x**4)
x0, tolg = 0.5, 1e-7
print("\nMaximizar f(x) = 8 + 5x - 3x^4 - 2x^6, x0 = %g, tolg = %g" % (x0, tolg))
r = newton_1d(f, df, ddf, x0, tolg, 20)
H = r.history
print("%3s %11s %11s %11s %11s %11s" % ("k", "x_k", "f(x_k)", "df(x_k)", "ddf(x_k)", "x_k+1"))
for row in H:
    print("%3d %11.7f %11.6f %+11.6f %11.5f %11.7f" % tuple(row[:6]))
Hx = ddf(r.x)                        # classificar: 1 avaliação extra de f''
print("Após %d iterações |f'| < 1e-7: x* = %.7f, f(x*) = %.6f, f''(x*) = %.4f < 0: máximo"
      % (r.nit, r.x, r.fx, Hx))
Tm = np.array([[0.5000000, 10.281250, 3.125000, -12.75000, 0.7450980],
               [0.7450980, 10.458621, -2.719687, -38.47906, 0.6744184],
               [0.6744184, 10.563259, -0.355311, -28.78702, 0.6620756],
               [0.6620756, 10.565490, -0.009182, -27.30912, 0.6617394],
               [0.6617394, 10.565491, -0.000007, -27.26970, 0.6617391]])
c = [H.shape[0] == 5 and confere(H[:, [1, 5]], Tm[:, [0, 4]], 7) and confere(H[:, 2:4], Tm[:, 1:3], 6)
     and confere(H[:, 4], Tm[:, 3], 5),
     r.nit == 5, confere(r.x, 0.6617391, 7), confere(r.fx, 10.565491, 6), Hx < 0]
print("  tabela k = 0..4: %s | 5 it.: %s | x*: %s | f(x*): %s | máximo: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ os quatro métodos
f = lambda x: 8 * np.exp(1 - x) + 7 * np.log(x)
df = lambda x: -8 * np.exp(1 - x) + 7 / x
ddf = lambda x: 8 * np.exp(1 - x) - 7 / x**2
r5 = newton_1d(f, df, ddf, 1, 1e-6, 20)          # para após 5 it.
print("\nOs quatro métodos, Newton de x0 = 1: erro após %d it. = %.1e (n_g = %d, n_H = %d)"
      % (r5.nit, abs(r5.x - xs), r5.ngev, r5.nhev))
# o slide escreve 1,3e-13
c = [r5.nit == 5, round(abs(r5.x - xs) / 1e-14) == 13, r5.ngev == 6, r5.nhev == 5]
print("  5 it.: %s | erro 1,3e-13: %s | n_g = 6: %s | n_H = 5: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
