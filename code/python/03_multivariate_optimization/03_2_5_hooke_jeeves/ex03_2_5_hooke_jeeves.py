"""Reproduz os exemplos do deck 3.2.5 (Hooke–Jeeves, algoritmo das aulas: Deb, 2012).

f(x) = 3 x1^2 + x2^2 - 12 x1 - 8 x2  (mínimo em (2,4), f = -28):
  «Uma exploração à mão»: de x0 = (1,1), f = -16, com Delta = (0.5, 0.5):
    x1: f+ = f(1.5,1) = -18.25, f- = f(0.5,1) = -12.25: fica (1.5,1);
    x2: f+ = f(1.5,1.5) = -21, f- = f(1.5,0.5) = -15: fica (1.5,1.5); 4 avaliações;
  «Exemplo completo»: x0 = (1,1), Delta = (0.5,0.5), R = 2, eps = 0.2:
    padrão (2,2) -> (2,2.5), sucesso; (2.5,3.5) -> (2,4), sucesso;
    (2,5.5) -> (2,5), insucesso: passo 3, Delta = (0.25,0.25); exploração em
    (2,4) falha; Delta = (0.125,0.125); falha; ||Delta|| = 0.177 < eps: pára em
    (2,4), f = -28, n_f = 28.
Himmelblau a partir de (0,0), Delta = (0.5,0.5), R = 2, eps = 0.2 (Deb, 2012,
  exercício 3.3.3): X1 = (0.5,0.5), X2 = (1.5,1.5), X3 = (3,2); o padrão
  (4.5,2.5) -> (4,2), f = 50, falha; pára em (3,2), f = 0.
Rosenbrock a partir de (-1.5, 2), Delta = 0.5, R = 2, eps = 1e-6: n_f = 627 até
  f < 1e-4 (com a versão de 1961 do código, variante="1961": 430).
Custo: 2n por exploração (avaliam-se +Delta_i e -Delta_i); 1 + 2n por
  movimento em padrão.

Contagens: f(x0) conta; f do ponto base nunca é reavaliado (fica guardado).
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
xb = np.array([1.0, 1.0]); fb = quad(xb); D = [0.5, 0.5]
print("\nExploração à mão: x0 = (%g,%g), f(x0) = %g, Delta = (%g,%g)" % (*xb, fb, *D))
r = exploratory(quad, xb, fb, D)
for h in r.history:
    print("  x%d %s Delta%d: (%g, %g) -> f%s = %g, %s"
          % (h[0], "+" if h[1] > 0 else "-", h[0], h[2], h[3], "+" if h[1] > 0 else "-", h[4],
             "a melhor: guarda-se" if h[5] else "não é a melhor"))
print("x1 = (%g, %g), f = %g, com %d avaliações" % (*r.x, r.fx, r.nfev))
c = [confere(fb, -16, 0),
     confere(r.history[:, 2:5], [[1.5, 1, -18.25], [0.5, 1, -12.25], [1.5, 1.5, -21], [1.5, 0.5, -15]], 2),
     confere(r.x, [1.5, 1.5], 1), r.nfev == 4]
print("  f(1,1) = -16: %s | -18.25, -12.25, -21, -15: %s | x1 = (1.5,1.5): %s | 4 avaliações: %s" % simnao(c))
ok = ok and all(c)

# --------------------------------------------- «Exemplo completo»
x0, D0, R, eps = [1, 1], [0.5, 0.5], 2, 0.2
print("\nExemplo completo: x0 = (%g,%g), Delta = (%g,%g), R = %g, eps = %g" % (*x0, *D0, R, eps))
r = hooke_jeeves(quad, x0, 2, D0, eps, R=R, verbose=True)
print(r.message)
print("solução x = (%g, %g), f = %g;  n_f = %d" % (*r.x, r.fx, r.nfev))
H = r.history
# movimentos 2 a 4 do slide (padrão): x^P(1) x^P(2) f(x^P)  x(1) x(2) f(x)  sucesso
Tpad = [[2.0, 2.0, -24.00, 2.0, 2.5, -25.75, 1],
        [2.5, 3.5, -27.00, 2.0, 4.0, -28.00, 1],
        [2.0, 5.5, -25.75, 2.0, 5.0, -27.00, 0]]
c = [r.nit == 6 and list(H[:, 1]) == [0, 1, 1, 1, 0, 0],   # explor., 3 padrões, 2 explor.
     confere(H[0, 7:10], [1.5, 1.5, -21], 2),
     confere(H[1:4][:, [4, 5, 6, 7, 8, 9, 13]], Tpad, 2),
     confere(H[4:6, 11], [0.25, 0.125], 4) and bool(np.all(H[4:6, 13] == 0)),
     confere(np.linalg.norm(r.P), 0.177, 3),
     confere(r.x, [2, 4], 4) and confere(r.fx, -28, 4), r.nfev == 28]
print("  sequência de movimentos: %s | 1.ª exploração -> (1.5,1.5): %s | padrões (2,2), (2.5,3.5), (2,5.5): %s |\n"
      "  falha com Delta = 0.25, 0.125: %s | ||Delta|| = 0.177 < eps: %s | (2,4), f = -28: %s | n_f = 28: %s" % simnao(c))
ok = ok and all(c)

# --------------------------------------------- Himmelblau (Deb, 2012, exercício 3.3.3)
himmel = lambda x: (x[0]**2 + x[1] - 11)**2 + (x[0] + x[1]**2 - 7)**2
print("\nHimmelblau a partir de (0,0), Delta = (0.5,0.5), R = 2, eps = 0.2")
r = hooke_jeeves(himmel, [0, 0], 2, [0.5, 0.5], 0.2, verbose=True)
print(r.message)
H = r.history
c = [confere(H[0, 7:10], [0.5, 0.5, 144.125], 3), confere(H[1, 7:10], [1.5, 1.5, 63.125], 3),
     confere(H[2, 7:10], [3, 2, 0], 3), confere(H[3, 4:10], [4.5, 2.5, 152.125, 4, 2, 50], 3) and H[3, 13] == 0,
     confere(r.x, [3, 2], 4) and r.fx == 0]
print("  X1 = (0.5,0.5): %s | X2 = (1.5,1.5): %s | X3 = (3,2): %s | padrão (4.5,2.5) -> (4,2), f = 50, falha: %s | (3,2), f = 0: %s" % simnao(c))
ok = ok and all(c)

# --------------------------------------------- Rosenbrock
x0, D0, eps = [-1.5, 2], [0.5, 0.5], 1e-6
print("\nRosenbrock a partir de (%g,%g), Delta = (%g,%g), R = 2, eps = %g" % (*x0, *D0, eps))
r = hooke_jeeves(rosen, x0, 2, D0, eps)
print(r.message)
hit = int(np.argmax(r.ftrace < 1e-4)) + 1
print("primeira avaliação com f < 1e-4: n_f = %d" % hit)
print("no fim: x = (%.6f, %.6f), f = %.2e;  %d movimentos, n_f total = %d"
      % (*r.x, r.fx, r.nit, r.nfev))
# custo de cada movimento: 2n (exploração), 1 + 2n (padrão), n = 2
H = r.history; n = 2
dn = np.diff(np.concatenate(([1], H[:, -1])))
ce = dn[H[:, 1] == 0]; cp = dn[H[:, 1] == 1]
print("custo por exploração: %d a %d; por movimento em padrão: %d a %d"
      % (ce.min(), ce.max(), cp.min(), cp.max()))
r61 = hooke_jeeves(rosen, x0, 2, D0, [1e-6, 1e-6], variante="1961")
hit61 = int(np.argmax(r61.ftrace < 1e-4)) + 1
print("(versão de 1961, variante=\"1961\": n_f = %d até f < 1e-4)" % hit61)
c = [hit == 627, hit61 == 430, bool(np.all(ce == 2 * n)), bool(np.all(cp == 1 + 2 * n))]
print("  n_f = 627 até f < 1e-4: %s | 1961: 430: %s | n_f = 2n por exploração: %s | 1 + 2n por padrão: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
