"""Reproduz os exemplos do deck 3.2.3 (simplex de Nelder–Mead).

Notação e paragem das aulas: x_l, x_g, x_h; gamma = 2, beta = 0.5;
análise do erro Q = sqrt(sum (f(x_i) - f(x_c))^2 / (n+1)), parar se Q <= eps.
1. Uma iteração à mão: f = x1^2 + x2^2, simplex (-1,0), (1,2), (-1,3):
   reflexão aceite, Q = 2.380, 2 avaliações (x_c e x_r).
2. Himmelblau a partir de (0,0), simplex (0,0), (1.2,0), (0,1.2), eps = 1e-6:
   tabela k = 0..12 (com Q); 31 iterações, n_f = 92, x* = (3,2).
   Com o simplex (0,0), (-1.2,0), (0,-1.2) chega-se a (-2.81; 3.13).
3. Rosenbrock a partir de (-1.5, 2): n_f até f < 1e-4 = 244 (simplex de
   arestas 0.5) e 211 (arestas de 5 %); sem as avaliações de f(x_c), isto é,
   com a paragem do Nelder–Mead padrão (tolx, tolf, como o fminsearch):
   166 e 139.

O deck mostra também fminsearch (MATLAB) / scipy.optimize.minimize (Python);
aqui usa-se a implementação da UC, nelder_mead, com o algoritmo do deck.
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex03_2_3_nelder_mead.py

Otimização — deck 3.2.3.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from nelder_mead import nelder_mead


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 3.2.3: método do simplex de Nelder–Mead")

sq = lambda x: x[0]**2 + x[1]**2
himmel = lambda x: (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2
ros = lambda x: (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

# ------------------------------------------------------------ 1. à mão
print("\n1. Uma iteração à mão: f = x1^2 + x2^2, simplex (-1,0), (1,2), (-1,3)")
r = nelder_mead(sq, [[-1, 0], [1, 2], [-1, 3]], kmax=1, eps=0, verbose=True)
S = r.simplex; Fs = [sq(x) for x in S]
print("novo simplex: (%g,%g), (%g,%g), (%g,%g) com f = %g, %g, %g;  Q = %.3f;  avaliações nesta iteração: %d"
      % (*S[0], *S[1], *S[2], *Fs, r.Q, r.nfev - 3))
c = [r.ops[1] == "reflexão", confere(S, [[-1, 0], [1, -1], [1, 2]], 4),
     confere(Fs, [1, 2, 5], 4), confere(r.Q, 2.380, 3), r.nfev - 3 == 2]
print("  reflexão aceite: %s | novo simplex: %s | f = 1, 2, 5: %s | Q = 2.380: %s | 2 avaliações: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 2. Himmelblau
print("\n2. Himmelblau a partir de (0,0), simplex (0,0), (1.2,0), (0,1.2), eps = 1e-6")
r = nelder_mead(himmel, [[0, 0], [1.2, 0], [0, 1.2]], eps=1e-6)
print("primeiras 13 linhas de info.history [k, n_f, f(x_l), f(x_h), Q, x_l] e operação:")
print("%4s %5s %12s %12s %10s %10s %10s  %s" % (*r.cols, "operação"))
for row, o in zip(r.history[:13], r.ops[:13]):
    print("%4d %5d %12.4f %12.4f %10.4f %10.4f %10.4f  %s" % (*row, o))
print("... %d iterações, n_f = %d, x* = (%.6f, %.6f), f = %.1e  (%s)"
      % (r.nit, r.nfev, r.x[0], r.x[1], r.fx, r.message))
# tabela do slide (k = 0..12): operação, x_l, f(x_l), f(x_h), Q
Tops = ["inicial", "expansão", "reflexão", "reflexão", "contr. int.", "contr. int.", "contr. ext.",
        "contr. int.", "contr. int.", "contr. int.", "contr. int.", "contr. ext.", "contr. int."]
Txl = [[1.200, 0.000], [1.800, 1.800], [3.000, 0.600], [3.000, 0.600], [2.550, 1.650],
       [3.188, 1.762], [3.188, 1.762], [3.188, 1.762], [2.884, 1.921], [2.919, 2.050],
       [2.919, 2.050], [3.031, 1.983], [3.031, 1.983]]
Tfl = [125.0336, 39.3632, 15.2096, 15.2096, 11.0925, 1.3499, 1.3499, 1.3499, 0.7628,
       0.1972, 0.1972, 0.0303, 0.0303]
Tfh = [170.000, 126.954, 125.034, 39.363, 24.579, 15.210, 11.093, 2.965, 1.604, 1.350,
       0.763, 0.217, 0.197]
TQ = [57.286, 52.768, 9.943, 14.104, 7.211, 4.981, 2.095, 1.286, 0.643, 0.262, 0.130, 0.104]
H = r.history[:13]
c = [r.ops[:13] == Tops, confere(H[:, 5:7], Txl, 3), confere(H[:, 2], Tfl, 4),
     confere(H[:, 3], Tfh, 3), confere(H[1:, 4], TQ, 3), r.nit == 31, r.nfev == 92,
     confere(r.x, [3, 2], 3), r.fx < 1e-5]
print("  tabela k = 0..12: operação: %s | x_l: %s | f(x_l): %s | f(x_h): %s | Q: %s\n"
      "  31 it.: %s | n_f = 92: %s | x* = (3,2): %s | f ~ 0: %s" % simnao(c))
ok = ok and all(c)

r2 = nelder_mead(himmel, [[0, 0], [-1.2, 0], [0, -1.2]], eps=1e-6)
print("simplex (0,0), (-1.2,0), (0,-1.2): x* = (%.2f; %.2f), f = %.1e" % (r2.x[0], r2.x[1], r2.fx))
c = [confere(r2.x, [-2.81, 3.13], 2)]
print("  chega a (-2.81; 3.13): %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 3. Rosenbrock
x0 = np.array([-1.5, 2.0])
print("\n3. Rosenbrock a partir de (-1.5, 2): n_f até f < 1e-4 (info.nhit)")
ra = nelder_mead(ros, [x0, x0 + [0.5, 0], x0 + [0, 0.5]], ftarget=1e-4, eps=1e-6)
rb = nelder_mead(ros, x0, ftarget=1e-4, eps=1e-6)          # simplex por omissão: h_i = 0.05 max(1,|x0_i|)
print("simplex de arestas 0.5:     n_f até f < 1e-4 = %d  (final: %d it., n_f = %d, x = (%.4f, %.4f))"
      % (ra.nhit, ra.nit, ra.nfev, *ra.x))
print("simplex de arestas de 5 %%:  n_f até f < 1e-4 = %d  (final: %d it., n_f = %d, x = (%.4f, %.4f))"
      % (rb.nhit, rb.nit, rb.nfev, *rb.x))
pa = nelder_mead(ros, [x0, x0 + [0.5, 0], x0 + [0, 0.5]], ftarget=1e-4)     # paragem padrão, sem f(x_c)
pb = nelder_mead(ros, x0, ftarget=1e-4)
print("sem as avaliações de f(x_c) (paragem do Nelder–Mead padrão): %d e %d" % (pa.nhit, pb.nhit))
c = [ra.nhit == 244, rb.nhit == 211, pa.nhit == 166, pb.nhit == 139]
print("  n_f = 244: %s | n_f = 211: %s | sem f(x_c): 166: %s, 139: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
