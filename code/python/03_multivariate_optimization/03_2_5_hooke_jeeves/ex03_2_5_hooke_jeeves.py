"""Reproduz os exemplos do deck 3.2.5 (Hooke–Jeeves).

f(x) = 3 x1^2 + x2^2 - 12 x1 - 8 x2  (mínimo em (2,4), f = -28):
  «Uma exploração à mão»: de xb = (1,1), f = -16, com P = (0.5, 0.5):
    x1: (1.5,1) -> -18.25, (0.5,1) -> -12.25: fica (1.5,1);
    x2: (1.5,1.5) -> -21, (1.5,0.5) -> -15: fica (1.5,1.5); 4 avaliações;
  «Exemplo completo»: x0 = (1,1), a = 2, P0 = (0.5,0.5), T = (0.1,0.1):
    padrão (2,2) -> (2,2.5), aceite; (2.5,3.5) -> (2,4), aceite;
    (2,5.5) -> (2,5), rejeitado; exploração em (2,4) falha com P = 0.5,
    0.25, 0.125; 0.0625 < T: para em (2,4), f = -28.
Rosenbrock a partir de (-1.5, 2), P0 = 0.5, a = 2, T = 1e-6 (como no
  caderno cap3_comparacao.ipynb): n_f = 430 até f < 1e-4.
Custo: 2n por exploração (avaliam-se +P_j e -P_j; Deb, 2012); 1 + 2n por
  movimento de padrão.

Contagens: f(x0) conta; f(xb) nunca é reavaliado (fica guardado).
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex03_2_5_hooke_jeeves.py

Otimização — deck 3.2.5.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from exploratory import exploratory
from hooke_jeeves import hooke_jeeves


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
quad = lambda x: 3 * x[0]**2 + x[1]**2 - 12 * x[0] - 8 * x[1]
rosen = lambda x: (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

print("Otimização — deck 3.2.5: método de Hooke–Jeeves")

# --------------------------------------------- «Uma exploração à mão»
xb = np.array([1.0, 1.0]); fb = quad(xb); P = [0.5, 0.5]
print("\nExploração à mão: xb = (%g,%g), f(xb) = %g, P = (%g,%g)" % (*xb, fb, *P))
r = exploratory(quad, xb, fb, P)
for h in r.history:
    print("  x%d %s P%d: (%g, %g) -> f = %g, %s"
          % (h[0], "+" if h[1] > 0 else "-", h[0], h[2], h[3], h[4],
             "a melhor: guarda-se" if h[5] else "não é a melhor"))
print("xe = (%g, %g), f = %g, com %d avaliações" % (*r.x, r.fx, r.nfev))
c = [confere(fb, -16, 0),
     confere(r.history[:, 2:5], [[1.5, 1, -18.25], [0.5, 1, -12.25], [1.5, 1.5, -21], [1.5, 0.5, -15]], 2),
     confere(r.x, [1.5, 1.5], 1), r.nfev == 4]
print("  f(1,1) = -16: %s | -18.25, -12.25, -21, -15: %s | xe = (1.5,1.5): %s | 4 avaliações: %s" % simnao(c))
ok = ok and all(c)

# --------------------------------------------- «Exemplo completo»
x0, a, P0, T = [1, 1], 2, [0.5, 0.5], [0.1, 0.1]
print("\nExemplo completo: x0 = (%g,%g), a = %g, P0 = (%g,%g), T = (%g,%g)" % (*x0, a, *P0, *T))
r = hooke_jeeves(quad, x0, a, P0, T, verbose=True)
print(r.message)
print("solução x = (%g, %g), f = %g;  n_f = %d" % (*r.x, r.fx, r.nfev))
H = r.history
# linhas 2 a 4 do slide (padrão): xt(1) xt(2) f(xt)  xe'(1) xe'(2) f(xe')  aceite
Tpad = [[2.0, 2.0, -24.00, 2.0, 2.5, -25.75, 1],
        [2.5, 3.5, -27.00, 2.0, 4.0, -28.00, 1],
        [2.0, 5.5, -25.75, 2.0, 5.0, -27.00, 0]]
c = [r.nit == 7 and list(H[:, 1]) == [0, 1, 1, 1, 0, 0, 0],   # explor., 3 padrões, 3 explor.
     confere(H[0, 7:10], [1.5, 1.5, -21], 2),
     confere(H[1:4][:, [4, 5, 6, 7, 8, 9, 13]], Tpad, 2),
     confere(H[4:7, 11], [0.5, 0.25, 0.125], 4) and bool(np.all(H[4:7, 13] == 0)),
     confere(r.P, [0.0625, 0.0625], 4),
     confere(r.x, [2, 4], 4) and confere(r.fx, -28, 4)]
print("  sequência de movimentos: %s | 1.ª exploração -> (1.5,1.5): %s | padrões (2,2), (2.5,3.5), (2,5.5): %s |\n"
      "  falha com P = 0.5, 0.25, 0.125: %s | P = 0.0625 < T: %s | (2,4), f = -28: %s" % simnao(c))
ok = ok and all(c)

# --------------------------------------------- Rosenbrock
x0, a, P0, T = [-1.5, 2], 2, [0.5, 0.5], [1e-6, 1e-6]
print("\nRosenbrock a partir de (%g,%g), a = %g, P0 = (%g,%g), T = (%g,%g)" % (*x0, a, *P0, *T))
r = hooke_jeeves(rosen, x0, a, P0, T)
print(r.message)
hit = int(np.argmax(r.ftrace < 1e-4)) + 1
print("primeira avaliação com f < 1e-4: n_f = %d" % hit)
print("no fim: x = (%.6f, %.6f), f = %.2e;  %d movimentos, n_f total = %d"
      % (*r.x, r.fx, r.nit, r.nfev))
# custo de cada movimento: 2n (exploração), 1 + 2n (padrão), n = 2
H = r.history; n = 2
dn = np.diff(np.concatenate(([1], H[:, -1])))
ce = dn[H[:, 1] == 0]; cp = dn[H[:, 1] == 1]
print("custo por exploração: %d a %d; por movimento de padrão: %d a %d"
      % (ce.min(), ce.max(), cp.min(), cp.max()))
r5 = hooke_jeeves(rosen, x0, a, P0, [1e-5, 1e-5])
print("(com T = 1e-5: n_f = %d até f < 1e-4 e n_f total = %d)"
      % (int(np.argmax(r5.ftrace < 1e-4)) + 1, r5.nfev))
c = [hit == 430, bool(np.all(ce == 2 * n)), bool(np.all(cp == 1 + 2 * n))]
print("  n_f = 430 até f < 1e-4: %s | n_f = 2n por exploração: %s | 1 + 2n por padrão: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
