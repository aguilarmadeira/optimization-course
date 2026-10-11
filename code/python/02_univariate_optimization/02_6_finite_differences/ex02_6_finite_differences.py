"""Reproduz os exemplos do deck 2.6 (derivadas numéricas; Newton numérico).

Exemplo-guia (1.º teste de 2019): f(x) = 8 e^(1-x) + 7 ln x.
Conta à mão: x = 1.5, h = 0.015 (= 0.01|x|): progressiva -0.17255,
    centrada -0.185605, f'' ~ 1.7411; passo de Newton x1 = 1.6066 (exato
    1.6066); com 6 casas f'' daria 1.7422, com 4 casas 1.7778.
«Então basta fazer h cada vez menor?»: tabela para h = 0.1 ... 1e-8.
Newton numérico de x0 = 1, h = 0.01|x|, tolg = 1e-4: estagna em 1.609407
    (x* = 1.609381), para em k = 4 e devolve x_5; n_f = 5*3 + 1 = 16.

Imprime as tabelas como nos slides, os valores finais e as contagens, e
no fim compara com os slides (vírgula decimal nos slides, ponto aqui).
Correr (de qualquer pasta):  python ex02_6_finite_differences.py

Otimização — deck 2.6.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from num_derivs import num_derivs
from newton_numeric import newton_numeric


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 2.6: derivadas numéricas; Newton numérico")

f = lambda x: 8 * np.exp(1 - x) + 7 * np.log(x)
dfex = lambda x: -8 * np.exp(1 - x) + 7 / x        # só para comparar
ddfex = lambda x: 8 * np.exp(1 - x) - 7 / x**2
xs = 1.6093810677230782                              # minimizante

# ---------------------------------------------------------- conta à mão
x = 1.5
df, ddf, info = num_derivs(f, x)
print("\nConta à mão: x = %g, h = %g (= 0.01|x|); exatos f' = %.5f, f'' = %.5f"
      % (x, info["h"], dfex(x), ddfex(x)))
print("  f(x-h) = %.8f, f(x) = %.8f, f(x+h) = %.8f   (n_f = %d)"
      % (info["fm"], info["f0"], info["fp"], info["nfev"]))
print("  progressiva %.5f (erro %.3f); regressiva %.5f; centrada %.6f (erro %.0e)"
      % (info["df_fwd"], abs(info["df_fwd"] - dfex(x)), info["df_bwd"], df, abs(df - dfex(x))))
print("  f''_num = %.4f (erro %.0e)" % (ddf, abs(ddf - ddfex(x))))
x1 = x - df / ddf
print("  passo de Newton com estes valores: x1 = %.4f (exato: %.4f)" % (x1, x - dfex(x) / ddfex(x)))
r6 = lambda v: np.round(v * 1e6) / 1e6    # valores com 6 casas
r4 = lambda v: np.round(v * 1e4) / 1e4    # valores com 4 casas
ddf6 = (r6(info["fp"]) - 2 * r6(info["f0"]) + r6(info["fm"])) / info["h"]**2
ddf4 = (r4(info["fp"]) - 2 * r4(info["f0"]) + r4(info["fm"])) / info["h"]**2
print("  com 6 casas (%.6f; %.6f; %.6f): f''_num = %.4f; com 4 casas: %.4f (cancelamento)"
      % (r6(info["fm"]), r6(info["f0"]), r6(info["fp"]), ddf6, ddf4))
c = [confere([info["fm"], info["f0"], info["fp"]], [7.69348098, 7.69050103, 7.68791283], 8),
     confere(info["df_fwd"], -0.17255, 5), confere(abs(info["df_fwd"] - dfex(x)), 0.013, 3),
     confere(df, -0.185605, 6), confere(abs(df - dfex(x)), 3e-5, 5),
     confere(ddf, 1.7411, 4), confere(x1, 1.6066, 4), confere(x - dfex(x) / ddfex(x), 1.6066, 4),
     confere(ddf6, 1.7422, 4), confere(ddf4, 1.7778, 4), info["nfev"] == 3]
print("  valores de f: %s | progressiva: %s | erro 0.013: %s | centrada: %s | erro 3e-5: %s | f'': %s | x1: %s | exato: %s | 6 casas: %s | 4 casas: %s | n_f = 3: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------- h cada vez menor?
print("\nh cada vez menor? (dupla precisão, x = 1.5)")
print("%8s %12s %9s %12s %9s" % ("h", "df central", "erro", "ddf central", "erro"))
hh = [0.1, 0.015, 1e-3, 1e-6, 1e-8]
R = np.zeros((len(hh), 4))
for i, h in enumerate(hh):
    d1, d2, _ = num_derivs(f, x, h)
    R[i] = [d1, abs(d1 - dfex(x)), d2, abs(d2 - ddfex(x))]
    print("%8.3g %12.6f %9.1e %12.4f %9.1e" % (h, *R[i]))
Sd = [-0.186738, -0.185605, -0.185579, -0.185579, -0.185579]
Se1 = [1.2e-3, 2.6e-5, 1.2e-7, 1e-10, 4e-8]; de1 = [4, 6, 8, 10, 8]
Sdd = [1.7382, 1.7411, 1.7411, 1.7408, 0.0000]
Se2 = [3e-3, 6e-5, 3e-7, 3e-4, 1.7]; de2 = [3, 5, 7, 4, 1]
c = [confere(R[:, 0], Sd, 6), confere(R[:, 1], Se1, de1), confere(R[:, 2], Sdd, 4),
     confere(R[:, 3], Se2, de2), R[4, 2] == 0]
print("  df: %s | erros df: %s | ddf: %s | erros ddf: %s | ddf = 0 com h = 1e-8: %s" % simnao(c))
ok = ok and all(c)
eps = np.finfo(float).eps
print("  ordens de grandeza do h ótimo: eps^(1/3) = %.0e (df), eps^(1/4) = %.1e (ddf)"
      % (eps ** (1 / 3), eps ** (1 / 4)))
c = [confere(eps ** (1 / 3), 6e-6, 6), confere(eps ** (1 / 4), 1.2e-4, 5)]
print("  eps^(1/3) ~ 6e-6: %s | eps^(1/4) ~ 1.2e-4: %s" % simnao(c))
ok = ok and all(c)

# ---------------------------------------------- Newton numérico, x0 = 1
x0, tolg, kmax = 1, 1e-4, 50
print("\nNewton numérico: f(x) = 8 e^(1-x) + 7 ln x, x0 = %g, h = max(0.01|x|, 1e-4), tolg = %g" % (x0, tolg))
r = newton_numeric(f, x0, tolg, kmax)
H = r.history
print("%3s %10s %8s %11s %10s %10s" % ("k", "x_k", "h", "df_num", "ddf_num", "x_k+1"))
for row in H:
    print("%3d %10.6f %8.4f %+11.4e %10.4f %10.6f" % tuple(row))
print("Para em k = %d e devolve x_%d = %.6f (estagna: x* = %.6f), f = %.7f" % (H[-1, 0], r.nit, r.x, xs, r.fx))
print("n_f = %d linhas x 3 = %d, + 1 (f(x_%d)) = %d" % (r.nit, 3 * r.nit, r.nit, r.nfev))
d3, _, _ = num_derivs(f, xs, 0.01 * xs)
print("Com h = 0.01 x*, f'_num(x*) = %.1e (não é 0)" % d3)
ra = newton_numeric(f, r.x, 1e-3, 1, 1e-4)
print("Reduzindo h para 1e-4|x| no fim, uma iteração: x = %.7f (n_f = %d)" % (ra.x, ra.nfev))
Tn = np.array([[1.000000, 0.0100, -0.9999, 0.9997, 2.000183],
               [2.000183, 0.0200, 0.5571, 1.1928, 1.533139],
               [1.533139, 0.0153, -0.1283, 1.7160, 1.607918],
               [1.607918, 0.0161, -0.0025, 1.6483, 1.609406],
               [1.609406, 0.0161, -1e-6, 1.6468, 1.609407]])
c = [H.shape[0] == 5 and confere(H[:, [1, 5]], Tn[:, [0, 4]], 6) and confere(H[:4, [2, 3, 4]], Tn[:4, [1, 2, 3]], 4)
     and confere(H[4, [2, 4]], Tn[4, [1, 3]], 4),
     confere(H[4, 3], -1e-6, 6), H[-1, 0] == 4, r.nit == 5, r.nfev == 16,
     confere(r.x, 1.609407, 6), confere(d3, -4e-5, 5), confere(ra.x, 1.6093811, 7)]
print("  tabela k = 0..4: %s | df_4 = -1e-6: %s | para em k = 4: %s | devolve x_5: %s | n_f = 16: %s | x = 1.609407: %s | f'_num(x*) = -4e-5: %s | 1.6093811: %s"
      % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
