"""Reproduz os exemplos do deck 2.5 (Newton-Raphson).

Exemplo-guia: f(x) = x^2 + 54/x, x0 = 2.5: x1 = 2.908; 4 iterações até
    x* = 3 (erro 2.4e-12); n_g = n_H = 1 por iteração (+ f'(x0)).
    Paragem |f'| < tolg com tolg = 1e-8 (o slide não fixa tolg; qualquer
    valor entre 1.5e-11 e 1.6e-5 dá as mesmas 4 iterações).
«Convergência quadrática»: x^3 - 2x^2 + 4, x0 = 1, tolg = 1e-3 (3 it.).
«Newton procura f' = 0»: maximizar 8 + 5x - 3x^4 - 2x^6, x0 = 0.5,
    tolg = 1e-7 (5 it.).
«Os quatro métodos»: x^2 + 54/x a partir de x0 = 1 (erro após 5 e 7 it.).

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
f = lambda x: x**2 + 54 / x
df = lambda x: 2 * x - 54 / x**2
ddf = lambda x: 2 + 108 / x**3
x0, tolg, kmax = 2.5, 1e-8, 20
print("\nExemplo-guia: f(x) = x^2 + 54/x, x0 = %g, tolg = %g" % (x0, tolg))
r = newton_1d(f, df, ddf, x0, tolg, kmax)
H = r.history
print("%3s %11s %11s %9s %11s %12s" % ("k", "x_k", "df(x_k)", "ddf(x_k)", "x_k+1", "|x_k+1 - 3|"))
for row in H:
    print("%3d %11.7f %+11.5f %9.4f %11.7f %12.1e" % (row[0], row[1], row[3], row[4], row[5], abs(row[5] - 3)))
print("x* = %.7f, f(x*) = %.4f; %d iterações, n_g = %d, n_H = %d, n_f = %d"
      % (r.x, r.fx, r.nit, r.ngev, r.nhev, r.nfev))
r5 = newton_1d(None, df, ddf, x0, 1e-15, 5)
print("A 5.ª iteração chega à precisão de máquina: x_5 - 3 = %.1e (%d iterações)" % (r5.x - 3, r5.nit))
Tg = np.array([[2.5000000, -3.64000, 8.9120, 2.9084381],
               [2.9084381, -0.56685, 6.3898, 2.9971495],
               [2.9971495, -0.01712, 6.0114, 2.9999973],
               [2.9999973, -0.00002, 6.0000, 3.0000000]])
eg = np.abs(H[:, 5] - 3)
c = [H.shape[0] == 4 and confere(H[:, [1, 5]], Tg[:, [0, 3]], 7) and confere(H[:, 3], Tg[:, 1], 5)
     and confere(H[:, 4], Tg[:, 2], 4),
     confere(eg, [9.2e-2, 2.9e-3, 2.7e-6, 2.4e-12], [3, 4, 7, 13]), confere(H[0, 5], 2.908, 3),
     r.nit == 4, r.ngev == 5, r.nhev == 4, r.nfev == 0, r5.x == 3]
print("  tabela k = 0..3: %s | erros: %s | x1 = 2.908: %s | 4 it.: %s | n_g = 5: %s | n_H = 4: %s | n_f = 0: %s | x_5 = 3: %s"
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
f = lambda x: x**2 + 54 / x
df = lambda x: 2 * x - 54 / x**2
ddf = lambda x: 2 + 108 / x**3
r5 = newton_1d(f, df, ddf, 1, 1e-2, 20)          # para após 5 it.
r7 = newton_1d(f, df, ddf, 1, 1e-10, 20)         # para após 7 it.
print("\nOs quatro métodos, Newton de x0 = 1: erro após %d it. = %.1e (n_g = %d, n_H = %d); após %d it. = %.1e"
      % (r5.nit, abs(r5.x - 3), r5.ngev, r5.nhev, r7.nit, abs(r7.x - 3)))
# o slide escreve 4,5e-4
c = [r5.nit == 5, round(abs(r5.x - 3) / 1e-5) == 45, r5.ngev == 6, r5.nhev == 5,
     r7.nit == 7, np.floor(np.log10(abs(r7.x - 3))) == -15]
print("  5 it.: %s | erro 4,5e-4: %s | n_g = 6: %s | n_H = 5: %s | 7 it.: %s | erro da ordem de 1e-15: %s"
      % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
