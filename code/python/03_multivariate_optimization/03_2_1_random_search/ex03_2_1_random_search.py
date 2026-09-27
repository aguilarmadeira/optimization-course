"""Reproduz os exemplos do deck 3.2.1 (pesquisa aleatória pura e localizada).

1. f = x1^2 + x2^2 em [-5,5]^2, uma corrida (semente 1): melhor f com
   N = 100 e N = 1000 amostras (slide 1: 1.406 e 0.053).
2. Análise experimental: 30 corridas por N (gerador com semente 7, partilhado
   pelas corridas, como no script das figuras): mediana e quartis (slide 3).
3. Localizada no Rosenbrock de (-1.5, 2), m = 20, r0 = 1, gamma = 0.9:
   uma corrida (semente 3), 150 iterações, n_f = 3001 (slide 5);
   30 corridas (sementes 0-29): n_f até f < 1e-4 = 975 [867; 1052].
4. Pura vs. localizada, orçamento 3000, 30 corridas (slide 6).

Método estocástico: usa o gerador numpy.random.default_rng e as mesmas
sementes dos scripts das figuras e do caderno do capítulo 3, por isso
reproduz EXATAMENTE os números dos slides.
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex03_2_1_random_search.py

Otimização — deck 3.2.1.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from random_search import random_search
from local_random_search import local_random_search


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


def mq(v):
    """mediana e quartis [Q1, Q3] (percentis 50, 25, 75, como nos scripts)."""
    return np.percentile(v, [50, 25, 75])


ok = True
print("Otimização — deck 3.2.1: pesquisa aleatória pura e localizada")
print("(Python: mesmas sementes e mesmo gerador dos scripts das figuras -> reprodução exata)")

sq = lambda x: x[0]**2 + x[1]**2
ros = lambda x: (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

# ------------------------------------------------------------ 1. uma corrida
a, b = [-5, -5], [5, 5]
print("\n1. f = x1^2 + x2^2 em [-5,5]^2, uma corrida (semente 1), N = 1000")
r = random_search(sq, a, b, 1000, s=1, verbose=True)
H = r.history
f100 = H[H[:, 0] <= 100, -1][-1]          # melhor f com as primeiras 100 amostras
print("N = 100: melhor f = %.3f;  N = 1000: melhor f = %.3f em (%.4f, %.4f);  n_f = %d"
      % (f100, r.fx, r.x[0], r.x[1], r.nfev))
c = [confere(f100, 1.406, 3), confere(r.fx, 0.053, 3), r.nfev == 1000]
print("  N = 100, f = 1.406: %s | N = 1000, f = 0.053: %s | n_f = N: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 2. 30 corridas
print("\n2. Análise experimental: 30 corridas por N (gerador com semente 7)")
g = np.random.default_rng(7)
tab = {}
for N in (10, 30, 100, 300, 1000, 3000, 10000):   # a ordem do script (o gerador é partilhado)
    fb = [random_search(sq, a, b, N, s=g).fx for _ in range(30)]
    tab[N] = mq(fb)
# slide: N, melhor f (mediana), ||x_best|| mediana [Q1; Q3]
T = {10: (2.3, 1, 1.5, 1.2, 2.1, 1), 100: (0.29, 2, 0.54, 0.38, 0.76, 2),
     1000: (0.018, 3, 0.14, 0.09, 0.17, 2), 10000: (0.0028, 4, 0.053, 0.027, 0.070, 3)}
print("%6s %10s %10s %10s %10s" % ("N", "f (med)", "||x|| med", "Q1", "Q3"))
c = []
for N, (fs, df, xm, x1, x3, dx) in T.items():
    med, q1, q3 = tab[N]; nx = np.sqrt([med, q1, q3])   # ||x|| = sqrt(f): comuta com os percentis
    print("%6d %10.4g %10.3g %10.3g %10.3g" % (N, med, *nx))
    c.append(confere(med, fs, df) and confere(nx, [xm, x1, x3], dx))
print("  tabela N = 10, 100, 1000, 10000: %s | %s | %s | %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 3. localizada
X0 = [-1.5, 2.0]
print("\n3. Localizada no Rosenbrock de (-1.5, 2), m = 20, r0 = 1, gamma = 0.9")
r = local_random_search(ros, X0, 1.0, m=20, gamma=0.9, kmax=150, s=3)
print("uma corrida (semente 3): %d iterações, n_f = %d, x = (%.6f, %.6f), f = %.2e"
      % (r.nit, r.nfev, r.x[0], r.x[1], r.fx))
print("  primeiras e últimas linhas de info.history [k, n_f, r, f(x), x1, x2]:")
for row in np.r_[r.history[:4], r.history[-2:]]:
    print("  %5d %6d %10.3e %12.4e %10.6f %10.6f" % tuple(row))
c = [r.nit == 150, r.nfev == 3001, confere(r.x, [1, 1], 4), confere(r.fx / 1e-13, 5.9, 1)]
print("  150 it.: %s | n_f = 3001: %s | x -> (1,1): %s | f = 5.9e-13: %s" % simnao(c))
ok = ok and all(c)

nhit = [local_random_search(ros, X0, 1.0, m=20, gamma=0.9, tolx=1e-8, kmax=5000,
                            ftarget=1e-4, s=s).nhit for s in range(30)]
med, q1, q3 = mq(nhit)
print("30 corridas (sementes 0-29): n_f até f < 1e-4 = %.0f [%.0f; %.0f] (mediana [Q1; Q3]); falhas: %d"
      % (med, q1, q3, int(np.sum(np.isinf(nhit)))))
c = [confere([med, q1, q3], [975, 867, 1052], 0)]
print("  975 [867; 1052]: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 4. pura vs. localizada
print("\n4. Pura vs. localizada no Rosenbrock, orçamento 3000, 30 corridas (sementes 0-29)")
fp = [random_search(ros, [-2, -1], [2, 3], 3000, s=s).fx for s in range(30)]
fl = [local_random_search(ros, X0, 1.0, m=20, gamma=0.9, kmax=150, s=s).fx for s in range(30)]
mp, ml = mq(fp), mq(fl)
print("pura em [-2,2]x[-1,3], n_f = 3000:  melhor f = %.3f [%.3f; %.3f]" % tuple(mp))
print("localizada, 150 it., n_f = 3001:    melhor f = %.1e [%.1e; %.1e]" % tuple(ml))
c = [confere(mp, [0.014, 0.008, 0.022], 3), confere(ml / [1e-12, 1e-13, 1e-12], [2, 5, 8], 0)]
print("  pura 0.014 [0.008; 0.022]: %s | localizada 2e-12 [5e-13; 8e-12]: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
