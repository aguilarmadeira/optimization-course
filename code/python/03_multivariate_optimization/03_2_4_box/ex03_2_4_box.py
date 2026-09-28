"""Reproduz os exemplos do deck 3.2.4 (método de Box, EVOP simplificado).

«Uma iteração à mão»: Himmelblau em x = (0,0), Delta = (2,2): f nos 4 vértices.
Himmelblau a partir de (0,0), Delta = (2,2), tolx = 1e-4: tabela k = 0..9,
  20 iterações, n_f = 4*20 + 1 = 81, termina em (3,2) com f = 0.
Rosenbrock a partir de (-1.5, 2), Delta_0 = (1,1): n_f = 13 445 até f < 1e-4
  (tolx = 1e-6, como no caderno cap3_comparacao.ipynb; o deck não indica tolx);
  a caixa tem de encolher até Delta ~ 6e-5.
«Porque falha em vales curvos?»: em (0,0), passo 0.25, as 4 diagonais dão
  f entre 4.08 e 11.33 (> 1), mas f(0.25, 0) = 0.953.

Contagens: f(x0) conta e em cada iteração avaliam-se SEMPRE os 2^n vértices
(nota do deck: «nas contagens deste deck avaliam-se sempre os 2^n vértices»).
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex03_2_4_box.py

Otimização — deck 3.2.4.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from box_evo import box_evo


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
himmel = lambda x: (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2
rosen = lambda x: (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

print("Otimização — deck 3.2.4: método de Box (EVOP simplificado)")

# ------------------------------------------------------------ Himmelblau
x0, Delta, tolx = [0, 0], [2, 2], 1e-4
print("\nHimmelblau a partir de (%g,%g), Delta = (%g,%g), tolx = %g" % (*x0, *Delta, tolx))
r = box_evo(himmel, x0, Delta, tolx, verbose=True)
print(r.message)
print("x = (%.4f, %.4f), f(x) = %.4g;  %d iterações, n_f = %d (= 4 x %d + 1)"
      % (*r.x, r.fx, r.nit, r.nfev, r.nit))

# «Uma iteração à mão»: vértices (-1,-1), (1,-1), (-1,1), (1,1), pela ordem do slide
# (a ordem de avaliação de box_evo é (-1,-1), (-1,1), (1,-1), (1,1))
fv = r.ftrace[[1, 3, 2, 4]]
print("Iteração à mão: f(0,0) = %g; vértices (-1,-1), (1,-1), (-1,1), (1,1): %g, %g, %g, %g"
      % (r.ftrace[0], *fv))

# tabela do slide (k = 0..9): x_k(1) x_k(2)  Delta  f(x_k)  decisão (1 = move)  f novo
T = [[0.0, 0.0, 2.000, 170.000, 1, 106.000],
     [1.0, 1.0, 2.000, 106.000, 1, 26.000],
     [2.0, 2.0, 2.000, 26.000, 1, 10.000],
     [3.0, 1.0, 2.000, 10.000, 0, 10.000],
     [3.0, 1.0, 1.000, 10.000, 1, 9.125],
     [3.5, 1.5, 1.000, 9.125, 1, 0.000],
     [3.0, 2.0, 1.000, 0.000, 0, 0.000],
     [3.0, 2.0, 0.500, 0.000, 0, 0.000],
     [3.0, 2.0, 0.250, 0.000, 0, 0.000],
     [3.0, 2.0, 0.125, 0.000, 0, 0.000]]
H = r.history[:10, [1, 2, 3, 5, 6, 7]]
c = [confere(r.ftrace[0], 170, 0), confere(fv, [170, 146, 130, 106], 0),
     confere(H, T, 3), r.nit == 20, r.nfev == 81,
     confere(r.x, [3, 2], 3), confere(r.fx, 0, 3)]
print("  f(0,0) = 170: %s | vértices 170/146/130/106: %s | tabela k = 0..9: %s | "
      "20 it.: %s | n_f = 81: %s | x = (3,2): %s | f = 0: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ Rosenbrock
x0, Delta, tolx = [-1.5, 2], [1, 1], 1e-6
print("\nRosenbrock a partir de (%g,%g), Delta_0 = (%g,%g), tolx = %g" % (*x0, *Delta, tolx))
r = box_evo(rosen, x0, Delta, tolx)
print(r.message)
hit = int(np.argmax(r.ftrace < 1e-4)) + 1
kh = int(np.ceil((hit - 1) / 4))       # linha de history (iteração k = kh - 1), a contar de 1
Dh = r.history[kh - 1, 3]              # Delta_1 nessa iteração
print("primeira avaliação com f < 1e-4: n_f = %d (na iteração k = %d, com Delta = %.2e)"
      % (hit, kh - 1, Dh))
print("no fim: x = (%.6f, %.6f), f = %.2e;  %d iterações, n_f total = %d"
      % (*r.x, r.fx, r.nit, r.nfev))
c = [hit == 13445, confere(Dh, 6e-5, 5)]
print("  n_f = 13 445 até f < 1e-4: %s | Delta ~ 6e-5: %s" % simnao(c))
ok = ok and all(c)

# --------------------------------------- «Porque falha em vales curvos?»
r = box_evo(rosen, [0, 0], [0.5, 0.5], 0, kmax=1)       # uma só iteração
fd = r.ftrace[1:5]; fe1 = rosen([0.25, 0])
nomes = ("encolhe", "move")
print("\nRosenbrock em (0,0), f = %g, passo 0.25: diagonais f = %s -> %s"
      % (r.ftrace[0], "".join("%.2f " % v for v in fd), nomes[int(r.history[0, 6])]))
print("eixo +e1: f(0.25, 0) = %.3f < 1" % fe1)
c = [bool(np.all(fd > 1)), confere(fd.min(), 4.08, 2), confere(fd.max(), 11.33, 2),
     r.history[0, 6] == 0, confere(fe1, 0.953, 3)]
print("  nenhuma diagonal melhora: %s | 4.08: %s | 11.33: %s | encolhe: %s | 0.953: %s"
      % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
