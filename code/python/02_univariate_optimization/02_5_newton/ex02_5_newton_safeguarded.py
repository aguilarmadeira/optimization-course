"""Reproduz o exemplo «Newton salvaguardado: bisseção + Newton» do deck 2.5.

f'(x) = arctan x (mínimo de f(x) = x arctan x - ln(1 + x^2)/2 em 0),
intervalo [-1, 2], x0 = 1.5 -- o caso em que Newton puro diverge
(página «Quando falha»).  Um passo de bisseção, depois três de Newton:

  k   x_k        passo de Newton               x_k+1        novo intervalo
  0   1.5        -1.694: fora de (a,b) -> bis.  0.5          [-1; 0.5]
  1   0.5        -0.0796: aceite                -0.0796      [-0.0796; 0.5]
  2   -0.0796    3.4e-4: aceite                 3.4e-4       [-0.0796; 3.4e-4]
  3   3.4e-4     -2.5e-11: aceite               -2.5e-11

Paragem |f'| < tolg = 1e-10 (a do script das figuras).
Imprime a tabela, as contagens e a comparação com Newton puro, e no fim
compara com os slides (vírgula decimal nos slides, ponto aqui).
Correr (de qualquer pasta):  python ex02_5_newton_safeguarded.py

Otimização — deck 2.5.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from newton_safeguarded import newton_safeguarded


def confere_sig(v, s, m):
    """v confere com s mostrado com m algarismos significativos?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    m = np.ravel(np.asarray(m, float))
    e = np.floor(np.log10(np.abs(s)))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (e - m + 1) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 2.5: Newton salvaguardado (bisseção + Newton)")

f = lambda x: x * np.arctan(x) - 0.5 * np.log(1 + x**2)
df = lambda x: np.arctan(x)
ddf = lambda x: 1 / (1 + x**2)

# ---------------------------------------------- Newton puro diverge
x0 = 1.5
print("\nNewton puro para f'(x) = arctan x a partir de x0 = %g (3 iterações):" % x0)
xs = [x0]
for _ in range(3):                    # x_{k+1} = x_k - f'(x_k)/f''(x_k)
    xs.append(xs[-1] - df(xs[-1]) / ddf(xs[-1]))
print("  x_k: " + " -> ".join("%.4f" % v for v in xs) + "  (|x_k| cresce: diverge)")
c = [abs(xs[3]) > abs(xs[2]) > abs(xs[1]) > abs(xs[0])]
print("  diverge: %s" % simnao(c))
ok = ok and all(c)

# ---------------------------------------------- Newton salvaguardado
a, b, tolg = -1.0, 2.0, 1e-10
print("\nNewton salvaguardado em [%g, %g], x0 = %g, tolg = %g:" % (a, b, x0, tolg))
r = newton_safeguarded(f, df, ddf, a, b, x0, tolg, 50, verbose=True)
H = r.history
print("%s; x* ~ %.1e, f(x*) = %.1e" % (r.message, r.x, r.fx))
print("n_g = %d (f'(a), f'(b), f'(x0) e uma por iteração), n_H = %d, n_f = %d" % (r.ngev, r.nhev, r.nfev))
c = [H.shape[0] == 4,
     list(H[:, 3]) == [1, 0, 0, 0],
     confere_sig(H[:, 1], [1.5, 0.5, -0.0796, 3.4e-4], [2, 1, 3, 2]),
     confere_sig(H[:, 2], [-1.694, -0.0796, 3.4e-4, -2.5e-11], [4, 3, 2, 2]),
     confere_sig(H[:, 4], [0.5, -0.0796, 3.4e-4, -2.5e-11], [1, 3, 2, 2]),
     confere_sig(H[:3, 5:7], [[-1, 0.5], [-0.0796, 0.5], [-0.0796, 3.4e-4]], [[1, 1], [3, 1], [3, 2]]),
     r.nbis == 1, r.ngev == 7, r.nhev == 4]
print("  4 iterações: %s | 1 bisseção + 3 Newton: %s | x_k: %s | passos de Newton: %s | x_k+1: %s"
      " | intervalos: %s | 1 bisseção: %s | n_g = 7: %s | n_H = 4: %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
