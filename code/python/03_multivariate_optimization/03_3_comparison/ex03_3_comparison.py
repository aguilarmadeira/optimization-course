"""Reproduz a tabela «Os métodos estudados no mesmo problema, com o custo em
avaliações» do fim do capítulo 3 (deck 3.3.3), com as funções da UC.

Regras (folha de convenções da UC):
* mesmo problema (Rosenbrock, fg.py), mesmo ponto inicial x0 = (-1.5; 2),
  mesmo critério: custo até à PRIMEIRA avaliação de f com f < 1e-4 (f nos
  iterandos dos métodos com derivadas também conta em n_f);
* custo em avaliações contadas à parte, n_f, n_g, n_H; equivalente em f com
  derivadas centrais (n = 2): n_f + 4 n_g + 4 n_H (a Hessiana, 2n^2+1 = 9
  pontos, reaproveita os 4 do gradiente e f(x_k));
* métodos de 3.3 com a pesquisa em linha precisa (Brent, tol 1e-10);
* pesquisa aleatória localizada (estocástica): 30 corridas, sementes 0..29
  (numpy.random.default_rng), mediana [quartis].

Valores dos slides (n_f / n_g / n_H / equiv.):
  aleatória localizada 975 [867; 1052]; Nelder-Mead 166; Box 13 445;
  Hooke-Jeeves 353; gradiente 60 900 / 3000 / - / 72 900 (não atinge
  f < 1e-4 em 3000 it.); FR (reinício n) 718 / 34 / - / 854;
  Newton puro 6 / 5 / 5 / 46; Newton amortecido 281 / 12 / 12 / 377.

Usa as funções das pastas 03_2_*, 03_3_* e common/ (postas no sys.path por
uc_setup).
Correr (de qualquer pasta):  python ex03_3_comparison.py   (alguns segundos)

Otimização — deck 3.3.3 (comparação do cap. 3).  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from box_evo import box_evo
from conj_grad import conj_grad
from fg import fg
from hooke_jeeves import hooke_jeeves
from local_random_search import local_random_search
from nelder_mead import nelder_mead
from newton_nd import newton_nd
from steepest_descent import steepest_descent


class Contador:
    """Conta n_f, n_g, n_H e guarda (n_f, n_g, n_H, f) em cada avaliação de f."""
    def __init__(self):
        self.nf = self.ng = self.nH = 0
        self.trace = []

    def F(self, x):
        v = fg(x); self.nf += 1
        self.trace.append((self.nf, self.ng, self.nH, v))
        return v

    def G(self, x):
        self.ng += 1
        return fg(x, 2)[1]

    def H(self, x):
        self.nH += 1
        return fg(x, 3)[2]

    def primeira(self, tol):
        """(n_f, n_g, n_H) na primeira avaliação com f < tol; None se nunca."""
        for nf, ng, nH, v in self.trace:
            if v < tol:
                return nf, ng, nH
        return None


def equiv(nf, ng, nH, n=2):
    """custo equivalente em f com derivadas centrais: n_f + 2n n_g + (2n^2 - 2n) n_H."""
    return nf + 2 * n * ng + (2 * n * n - 2 * n) * nH


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


X0 = np.array([-1.5, 2.0])
FTOL = 1e-4
ok = True
print("Otimização — capítulo 3: os métodos estudados no mesmo problema (Rosenbrock de (-1.5; 2))")
print("Critério: primeira avaliação com f < %g; custo equivalente em f = n_f + 4 n_g + 4 n_H" % FTOL)

metodos = [
    ("Nelder-Mead", 0, lambda C: nelder_mead(C.F, [X0, X0 + [0.5, 0], X0 + [0, 0.5]], 1e-8, 1e-10, 1000)),
    ("Box (Delta_0 = 1)", 0, lambda C: box_evo(C.F, X0, [1, 1], 1e-6, 10000)),
    ("Hooke-Jeeves (P_0 = 0.5)", 0, lambda C: hooke_jeeves(C.F, X0, 2, 0.5, 1e-6, 100000)),
    ("Gradiente", 1, lambda C: steepest_descent(C.F, C.G, X0, 1e-8, 3000)),
    ("Grad. conj. (FR, reinício n)", 1, lambda C: conj_grad(C.F, C.G, X0, 1e-8, 3000)),
    ("Newton puro", 2, lambda C: newton_nd(C.F, C.G, C.H, X0, 1e-8, 200)),
    ("Newton amortecido", 2, lambda C: newton_nd(C.F, C.G, C.H, X0, 1e-8, 200, damped=True, modify=True)),
]
slides = {"Nelder-Mead": (166, 0, 0, 166), "Box (Delta_0 = 1)": (13445, 0, 0, 13445),
          "Hooke-Jeeves (P_0 = 0.5)": (353, 0, 0, 353), "Gradiente": (60900, 3000, 0, 72900),
          "Grad. conj. (FR, reinício n)": (718, 34, 0, 854), "Newton puro": (6, 5, 5, 46),
          "Newton amortecido": (281, 12, 12, 377)}

# ------------------------------------------ pesquisa aleatória localizada
nhit = []
for s in range(30):
    C = Contador()
    local_random_search(C.F, X0, 1.0, 20, 0.9, 1e-8, 5000, s=s)
    hit = C.primeira(FTOL)
    nhit.append(hit[0] if hit else np.inf)
nhit = np.array(nhit)
q1, med, q3 = np.percentile(nhit[np.isfinite(nhit)], [25, 50, 75])
nfalhas = int(np.sum(~np.isfinite(nhit)))

print("\n%-30s %6s %18s %6s %5s %8s" % ("Método", "ordem", "n_f", "n_g", "n_H", "equiv."))
print("%-30s %6d %18s %6s %5s %8.0f" % ("Aleatória localizada (30x)", 0, "%.0f [%.0f; %.0f]" % (med, q1, q3),
                                         "--", "--", med))
linhas_ok = [round(med) == 975 and round(q1) == 867 and round(q3) == 1052 and nfalhas == 0]   # q3 = 1051.75

resultados = {}
for nome, ordem, corre in metodos:
    C = Contador()
    r = corre(C)
    hit = C.primeira(FTOL)
    # as contagens do contador coincidem com as da própria função?
    proprias = (r.nfev, getattr(r, "ngev", 0), getattr(r, "nhev", 0))
    if proprias != (C.nf, C.ng, C.nH):
        print("  AVISO: %s conta (%d, %d, %d), o contador (%d, %d, %d)" % (nome, *proprias, C.nf, C.ng, C.nH))
        ok = False
    if hit is None:                                 # não atinge: custo total da corrida
        nf, ng, nH = C.nf, C.ng, C.nH
        marca = "‡"
    else:
        nf, ng, nH = hit
        marca = ""
    eq = equiv(nf, ng, nH)
    resultados[nome] = (nf, ng, nH, eq, r)
    print("%-30s %6d %18d %6s %5s %8d" % (nome + marca, ordem, nf, ng if ng else "--", nH if nH else "--", eq))
    linhas_ok.append((nf, ng, nH, eq) == slides[nome])
print("‡ não atinge f < 1e-4 em 3000 iterações (f = %.1e); custo total da corrida"
      % resultados["Gradiente"][4].fx)

nd = resultados["Newton amortecido"]; fr = resultados["Grad. conj. (FR, reinício n)"]
print("\nDeck 3.3.3: até f < 1e-4, Newton amortecido n_g = n_H = %d, n_f = %d na pesquisa em linha (%d com f nos iterandos);"
      % (nd[1], nd[0] - nd[1], nd[0]))
print("            FR n_g = %d, n_f = %d na pesquisa em linha (%d)."
      % (fr[1], fr[0] - fr[1], fr[0]))
print("            (a primeira f < 1e-4 surge numa pesquisa em linha; até aí os iterandos avaliados são tantos quantos os gradientes)")
c2 = [nd[1] == 12 and nd[0] - nd[1] == 269, fr[1] == 34 and fr[0] - fr[1] == 684]
print("  Newton: 269 + 12 = 281: %s | FR: 684 + 34 = 718: %s" % simnao(c2))

print("\n  aleatória 975 [867; 1052], 30/30: %s | NM 166: %s | Box 13 445: %s | HJ 353: %s | gradiente 72 900: %s"
      " | FR 854: %s | Newton puro 46: %s | amortecido 377: %s" % simnao(linhas_ok))
ok = ok and all(linhas_ok) and all(c2)

# ------------------------------------------------ o mesmo com o SciPy (informativo)
try:
    from scipy.optimize import minimize
    print("\nO mesmo com o SciPy (informativo: critério e pesquisa em linha do SciPy, números diferentes):")
    f1 = lambda x: fg(x)
    g1 = lambda x: fg(x, 2)[1]
    h1 = lambda x: fg(x, 3)[2]
    for metodo, kw in (("Nelder-Mead", {}), ("CG", {"jac": g1}), ("BFGS", {"jac": g1}),
                       ("Newton-CG", {"jac": g1, "hess": h1})):
        r = minimize(f1, X0, method=metodo, tol=1e-8, **kw)
        print("  %-12s f = %.1e  nfev = %4d  njev = %4d  nhev = %4d"
              % (metodo, r.fun, r.nfev, getattr(r, "njev", 0), getattr(r, "nhev", 0)))
except ImportError:
    print("\n(SciPy não instalado: parte informativa omitida)")

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
