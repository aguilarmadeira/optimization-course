"""Reproduz os exemplos do deck 3.2.2 (pesquisa em grelha).

Himmelblau em [-5,5]^2:
1. Grelha 10 x 10 (100 avaliações): melhor f = 4.51; os 4 melhores nós,
   f = 4.51, 4.70, 6.15, 11.88, ficam um em cada bacia.
2. Grosseira -> fina: mais 100 avaliações (grelha 10 x 10 em [x - h, x + h],
   h = 10/9) à volta do melhor nó: f = 0.32.
3. Grelha + Nelder–Mead (3.2.3) a partir de cada um dos 4 melhores nós
   (simplex do deck, tolerâncias 1e-6): n_f = 89, 94, 94, 87; os 4 mínimos
   globais, todos com f < 1e-11, em 464 avaliações no total; o 5.º melhor
   nó leva outra vez a (3,2), com mais 95 avaliações.

Usa nelder_mead da pasta 03_2_3_nelder_mead (posta no sys.path por
uc_setup; não há cópia do módulo nesta pasta).
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex03_2_2_grid_search.py

Otimização — deck 3.2.2.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from grid_search import grid_search
from nelder_mead import nelder_mead  # (3.2.3)


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 3.2.2: pesquisa em grelha")

himmel = lambda x: (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2
MIN = np.array([[3, 2], [-2.805118, 3.131312], [-3.779310, -3.283186], [3.584428, -1.848126]])
bacia = lambda p: int(np.argmin(np.linalg.norm(MIN - p, axis=1))) + 1   # mínimo mais próximo

# ------------------------------------------------------------ 1. grelha 10 x 10
print("\n1. Himmelblau em [-5,5]^2, grelha 10 x 10")
r = grid_search(himmel, [-5, -5], [5, 5], 10, verbose=True)
top4 = r.history[:4]
b4 = [bacia(p) for p in top4[:, :2]]
print("melhor nó (%.4f, %.4f), f = %.2f;  n_f = %d;  bacias dos 4 melhores: %d %d %d %d"
      % (*r.x, r.fx, r.nfev, *b4))
c = [confere(r.fx, 4.51, 2), r.nfev == 100, confere(top4[:, 2], [4.51, 4.70, 6.15, 11.88], 2),
     len(set(b4)) == 4]
print("  f = 4.51: %s | n_f = 100: %s | 4 melhores: %s | um em cada bacia: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 2. grosseira -> fina
h = r.h
r2 = grid_search(himmel, r.x - h, r.x + h, 10)
print("\n2. Grosseira -> fina: grelha 10 x 10 em [x - h, x + h], h = 10/9 = %.4f" % h[0])
print("melhor nó (%.4f, %.4f), f = %.2f;  mais %d avaliações" % (*r2.x, r2.fx, r2.nfev))
c = [confere(r2.fx, 0.32, 2), r2.nfev == 100]
print("  f = 0.32: %s | +100 avaliações: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 3. grelha + Nelder–Mead
print("\n3. Grelha + Nelder–Mead a partir dos 5 melhores nós (tolx = tolf = 1e-6)")
print("  nó      ponto inicial         f                 mínimo    f final   n_f")
tot = r.nfev; nf = []; ff = []; bm = []
for i in range(5):
    p = r.history[i, :2]
    rm = nelder_mead(himmel, p, tolx=1e-6, tolf=1e-6)   # simplex do deck a partir de p
    nf.append(rm.nfev); ff.append(rm.fx); bm.append(bacia(rm.x))
    print("%3d  (%7.4f, %7.4f) %9.2f  (%9.6f, %9.6f) %10.1e %5d"
          % (i + 1, *p, r.history[i, 2], *rm.x, rm.fx, rm.nfev))
    if i < 4:
        tot += rm.nfev
print("total grelha + 4 Nelder–Mead: %d avaliações;  o 5.º nó leva ao mínimo %d, com mais %d"
      % (tot, bm[4], nf[4]))
c = [nf[:4] == [89, 94, 94, 87], all(v < 1e-11 for v in ff[:4]), len(set(bm[:4])) == 4,
     tot == 464, bm[4] == 1, nf[4] == 95]
print("  n_f = 89, 94, 94, 87: %s | f < 1e-11: %s | 4 mínimos diferentes: %s | total 464: %s\n"
      "  5.º nó -> (3,2): %s | mais 95: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
