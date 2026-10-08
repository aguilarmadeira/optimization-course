"""Reproduz os exemplos do deck 3.3.1 (método do gradiente) e o do gradiente
numérico do deck 3.3.

«Uma iteração à mão»: f(x) = (x1^2 + 9 x2^2)/2, x0 = (9, 1): alpha_0 = 0.2,
    x1 = (7.2; -0.8), f = 28.8, g1 = (7.2; -7.2), g1^T d0 = 0.
«O percurso em ziguezague»: tabela k = 0..8; 74 iterações para ||g|| < 1e-6,
    n_g = 75; f multiplica por 0.64 em cada iteração.
«No computador»: pesquisa em linha de alta precisão (a das figuras, cerca de
    24 avaliações de f por pesquisa) e do tipo fminbnd (6 por pesquisa): as
    mesmas 74 iterações.
Deck 3.3: gradiente de Rosenbrock em (-1.5; 2) = (-155, -50); com h = 1e-4
    o erro das diferenças centrais é 6e-6.
«Rosenbrock»: de (-1.5; 2), 3000 iterações, ainda f = 2.4e-4, x ~ (0.98; 0.97).

Imprime as tabelas, os valores finais e as contagens, e no fim compara com
os slides (vírgula decimal nos slides, ponto aqui).
Correr (de qualquer pasta):  python ex03_3_1_steepest_descent.py

Otimização — deck 3.3.1.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from grad_fd import grad_fd
from steepest_descent import steepest_descent


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


ok = True
print("Otimização — deck 3.3.1: método do gradiente")

f = lambda x: 0.5 * (x[0]**2 + 9 * x[1]**2)       # = x^T A x / 2, A = diag(1, 9)
grad = lambda x: np.array([x[0], 9 * x[1]])
x0 = np.array([9.0, 1.0])

# ------------------------------------------------ uma iteração à mão
g0 = grad(x0); d0 = -g0
r1 = steepest_descent(f, grad, x0, 1e-6, 1)          # uma só iteração
x1 = r1.x; g1 = grad(x1); alpha0 = r1.history[1, 5]
print("\nUma iteração à mão: g0 = (%g, %g), d0 = (%g, %g)" % (*g0, *d0))
print("alpha_0 = %.4f (exato: 162/810 = 0.2), x1 = (%.4f; %.4f), f(x1) = %.4f (era %g)"
      % (alpha0, *x1, f(x1), f(x0)))
print("g1 = (%.4f; %.4f), g1^T d0 = %.1e (ortogonais)" % (*g1, g1 @ d0))
c = [confere(alpha0, 0.2, 4), confere(x1, [7.2, -0.8], 4), confere(f(x1), 28.8, 4),
     confere(g1, [7.2, -7.2], 4), abs(g1 @ d0) < 1e-6]
print("  alpha_0 = 0.2: %s | x1: %s | f = 28.8: %s | g1: %s | g1^T d0 = 0: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ percurso em ziguezague
tolg = 1e-6
print("\nZiguezague: f(x) = (x1^2 + 9 x2^2)/2, x0 = (9, 1), tolg = %g, pesquisa em linha de Brent (tol 1e-10)" % tolg)
r = steepest_descent(f, grad, x0, tolg, 1000)
H = r.history
print("%4s %22s %12s %10s %9s" % ("k", "x_k", "f(x_k)", "||g_k||", "alpha_k-1"))
for row in H[:9]:
    print("%4d  (%9.4f; %9.4f) %12.5f %10.4f %9.4f" % tuple(row))
print("...")
print("%4d  (%9.2e; %9.2e) %12.2e %10.2e %9.4f" % tuple(H[-1]))
razao = H[1:, 3] / H[:-1, 3]
print("%s: %d iterações, n_g = %d, n_f = %d (%d nas pesquisas em linha + %d iteradas)"
      % (r.message, r.nit, r.ngev, r.nfev, r.nfev_ls, r.nit + 1))
print("média de %.1f avaliações de f por pesquisa em linha (o slide diz «cerca de 24»)" % (r.nfev_ls / r.nit))
print("f_{k+1}/f_k: mín. %.4f, máx. %.4f (cota de Kantorovich ((9-1)/(9+1))^2 = 0.64)" % (razao.min(), razao.max()))
T = np.array([[9.0000, 1.0000, 45.00000, 12.7279, np.nan],
              [7.2000, -0.8000, 28.80000, 10.1823, 0.2000],
              [5.7600, 0.6400, 18.43200, 8.1459, 0.2000],
              [4.6080, -0.5120, 11.79648, 6.5167, 0.2000],
              [3.6864, 0.4096, 7.54975, 5.2134, 0.2000],
              [2.9491, -0.3277, 4.83184, 4.1707, 0.2000],
              [2.3593, 0.2621, 3.09238, 3.3365, 0.2000],
              [1.8874, -0.2097, 1.97912, 2.6692, 0.2000],
              [1.5099, 0.1678, 1.26664, 2.1354, 0.2000]])
c = [confere(H[:9, 1:3], T[:, 0:2], 4) and confere(H[:9, 3], T[:, 2], 5) and confere(H[:9, 4], T[:, 3], 4)
     and confere(H[1:9, 5], T[1:, 4], 4),
     r.nit == 74, r.ngev == 75, r.flag == 0, confere([razao.min(), razao.max()], [0.64, 0.64], 4),
     round(r.nfev_ls / r.nit) == 24]
print("  tabela k = 0..8: %s | 74 it.: %s | n_g = 75: %s | ||g|| < 1e-6: %s | razão 0.64: %s | cerca de 24 por pesquisa: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------ pesquisa em linha tipo fminbnd
rb = steepest_descent(f, grad, x0, tolg, 1000, ls="fminbnd")
print("\nCom a pesquisa do tipo fminbnd em [0, 1] (TolX = 1e-4): %d iterações, n_g = %d, n_f = %d; %.1f avaliações por pesquisa"
      % (rb.nit, rb.ngev, rb.nfev, rb.nfev_ls / rb.nit))
c = [rb.nit == 74, rb.ngev == 75, rb.nfev_ls / rb.nit == 6]
print("  74 it.: %s | n_g = 75: %s | 6 avaliações por pesquisa: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ gradiente numérico
rn = steepest_descent(f, None, x0, tolg, 1000)
print("\nCom gradiente numérico (grad = None: grad_fd, 2n = 4 avaliações de f por gradiente):")
print("%d iterações, n_g = %d, n_f = %d = %d (pesquisas) + %d (iteradas) + 4 x %d (gradientes)"
      % (rn.nit, rn.ngev, rn.nfev, rn.nfev_ls, rn.nit + 1, rn.nit + 1))
c = [rn.nit == 74, rn.ngev == 0, rn.nfev == rn.nfev_ls + (rn.nit + 1) + 4 * (rn.nit + 1)]
print("  74 it.: %s | n_g = 0: %s | 2n por gradiente em n_f: %s" % simnao(c))
ok = ok and all(c)

rosen = lambda x: (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2
grosen = lambda x: np.array([-2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0]**2), 200 * (x[1] - x[0]**2)])
xr = np.array([-1.5, 2.0])
ge = grosen(xr)
print("\nDeck 3.3 — gradiente de Rosenbrock em (-1.5; 2): exato (%g, %g)" % tuple(ge))
erros = []
for h in (1e-2, 1e-4, 1e-6):
    gn, nf = grad_fd(rosen, xr, h)
    erros.append(np.linalg.norm(gn - ge))
    print("  h = %g: (%.6f, %.6f), erro %.1e, %d avaliações de f" % (h, *gn, erros[-1], nf))
gd, nf = grad_fd(rosen, xr)
print("  h por omissão (max(0.01|x_i|, 1e-4)): erro %.1e" % np.linalg.norm(gd - ge))
c = [confere(ge, [-155, -50], 0), confere(erros[1], 6e-6, 6), nf == 4]
print("  grad = (-155, -50): %s | erro 6e-6 com h = 1e-4: %s | 2n = 4 avaliações: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ Rosenbrock
print("\nRosenbrock de (-1.5; 2), tolg = 1e-4, kmax = 3000:")
rr = steepest_descent(rosen, grosen, xr, 1e-4, 3000)
print("%s\nx = (%.4f; %.4f), f = %.2e, ||g|| no último gradiente = %.2e; n_g = %d, n_f = %d"
      % (rr.message, *rr.x, rr.fx, rr.history[-2, 4], rr.ngev, rr.nfev))
c = [rr.nit == 3000, rr.flag == 1, confere(rr.fx, 2.4e-4, 5), confere(rr.x, [0.98, 0.97], 2)]
print("  3000 it.: %s | não converge: %s | f = 2.4e-4: %s | x ~ (0.98; 0.97): %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
