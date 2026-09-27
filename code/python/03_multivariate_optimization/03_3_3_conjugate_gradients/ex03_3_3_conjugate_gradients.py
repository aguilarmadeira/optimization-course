"""Reproduz os exemplos do deck 3.3.3 (gradientes conjugados, Fletcher-Reeves).

Quadrática f(x) = (x1^2 + 9 x2^2)/2 de x0 = (9, 1), pesquisa em linha exata:
    o 1.º passo é o do gradiente, x1 = (7.2; -0.8); beta_1 = 103.68/162 = 0.64,
    d1 = (-12.96; 1.44), d0^T A d1 = 0 (A-conjugadas), alpha_1 = 5/9,
    x2 = (0, 0): 2 iterações, n_g = 3.
Rosenbrock de (-1.5; 2), ||g|| < 1e-4: com reinício a cada n = 2, 36 it.
    (n_g = 37); sem reinício, 102 it.
«No computador»: com uma pesquisa do tipo fminbnd de tolerância por
    omissão, FR no Rosenbrock pode não convergir (informativo).

Imprime as tabelas, os valores finais e as contagens, e no fim compara com
os slides (vírgula decimal nos slides, ponto aqui).
Correr (de qualquer pasta):  python ex03_3_3_conjugate_gradients.py

Otimização — deck 3.3.3.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from conj_grad import conj_grad


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 3.3.3: gradientes conjugados (Fletcher-Reeves)")

# ------------------------------------------------ quadrática
A = np.diag([1.0, 9.0])
quad = lambda x: 0.5 * (x[0]**2 + 9 * x[1]**2)     # = x^T A x / 2
gq = lambda x: np.array([x[0], 9 * x[1]])
x0 = np.array([9.0, 1.0])
print("\nQuadrática (x1^2 + 9 x2^2)/2 de (9, 1), tolg = 1e-6, pesquisa em linha de Brent (tol 1e-10):")
r = conj_grad(quad, gq, x0, 1e-6, 100, verbose=True)
H = r.history                          # [k, x1, x2, f, ||g||, alpha, beta, d1, d2]
g0, g1 = gq(H[0, 1:3]), gq(H[1, 1:3])
d0, d1 = H[0, 7:9], H[1, 7:9]
print("||g1||^2/||g0||^2 = %.2f/%.0f = %.4f; d0^T A d1 = %.1e; alpha_1 = %.6f (5/9 = %.6f)"
      % (g1 @ g1, g0 @ g0, H[1, 6], d0 @ A @ d1, H[2, 5], 5 / 9))
print("%s; n_g = %d, n_f = %d" % (r.message, r.ngev, r.nfev))
c = [confere(H[1, 1:3], [7.2, -0.8], 4), confere(g1 @ g1, 103.68, 2), confere(H[1, 6], 0.64, 4),
     confere(d1, [-12.96, 1.44], 2), abs(d0 @ A @ d1) < 1e-4, confere(H[2, 5], 5 / 9, 6),
     confere(r.x, [0, 0], 6), r.nit == 2, r.ngev == 3]
print("  x1: %s | ||g1||^2 = 103.68: %s | beta_1 = 0.64: %s | d1: %s | d0'Ad1 = 0: %s | alpha_1 = 5/9: %s"
      " | x2 = (0, 0): %s | 2 it.: %s | n_g = 3: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ Rosenbrock
rosen = lambda x: (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2
grosen = lambda x: np.array([-2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0]**2), 200 * (x[1] - x[0]**2)])
xr = np.array([-1.5, 2.0])
print("\nRosenbrock de (-1.5; 2), tolg = 1e-4:")
rc = conj_grad(rosen, grosen, xr, 1e-4, 500)                # reinício a cada n = 2
rs = conj_grad(rosen, grosen, xr, 1e-4, 500, restart=0)     # sem reinício periódico
print("FR com reinício a cada n = 2: %d it., n_g = %d, n_f = %d, f = %.1e (%d reinícios)"
      % (rc.nit, rc.ngev, rc.nfev, rc.fx, rc.nrestart))
print("FR sem reinício periódico:    %d it., n_g = %d, n_f = %d, f = %.1e (%d reinícios por g^T d >= 0)"
      % (rs.nit, rs.ngev, rs.nfev, rs.fx, rs.nrestart))
c = [rc.nit == 36, rc.ngev == 37, rs.nit == 102, rc.flag == 0 and rs.flag == 0]
print("  36 it. com reinício: %s | n_g = 37: %s | 102 it. sem: %s | convergem: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------- pesquisa tipo fminbnd (informativo)
rb = conj_grad(rosen, grosen, xr, 1e-4, 500, ls="fminbnd")
print("\nCom a pesquisa do tipo fminbnd (TolX = 1e-4) e reinício a cada n: %s; f = %.2g"
      % (rb.message, rb.fx))
print("  (o slide diz que «pode não convergir»: %s)" % ("não convergiu" if rb.flag else "convergiu"))

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
