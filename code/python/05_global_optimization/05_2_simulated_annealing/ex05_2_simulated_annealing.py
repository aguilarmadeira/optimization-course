"""Reproduz os exemplos do deck 5.2 (simulated annealing).

Função de referência do cap. 5: f(x) = 0.1 x^2 + 3 sin 2x + 2 cos(3.3x + 1)
em X = [-8, 8], mínimo global x* = 2.4681 (f = -4.2381).

1. «Uma execução passo a passo»: x0 = 6, T0 = 3, c = 0.98, sigma = 0.5,
   Nmax = 400, semente 3: T e x_k em k = 30, 100, 200; o salto em k = 315
   (de 4.73 para 3.07, a 3.3 sigma); x = 2.464, n_f = 401.
   E «SA -> método local»: de x_SA, uma pesquisa direta chega a |x - x*| < 1e-6
   com 37 avaliações.
2. «Uma execução não chega: 100 corridas»: a tabela de 6 linhas (sucesso
   |x - x*| < 0.3, sempre de x0 = 6, semente 99): 35, 19, 74, 20, 88, 73 %.
3. «Em duas dimensões: Rastrigin»: de (4, -4), T0 = 20, c = 0.995,
   sigma = 0.5, Nmax = 3000, semente 3: f = 0.068, |x| = 0.02, n_f = 3001;
   30 corridas de pontos aleatórios (semente 11), sucesso f < 1: 97 % com
   c = 0.9, 0.95, 0.99; 100 % com 0.995; 90 % com 0.999.
   E, a título informativo, a chamada do slide «No computador» (sementes 1..30).

Os números aleatórios são os dos scripts das figuras (make_figs_5.py): o
mesmo gerador numpy (default_rng), as mesmas sementes e a mesma ordem de
consumo, por isso os valores dos slides reproduzem-se EXATAMENTE.
Demora alguns segundos (cerca de 10^6 iterações ao todo).
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex05_2_simulated_annealing.py

Otimização — deck 5.2.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from simulated_annealing import simulated_annealing


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


def rastrigin(x):
    x = np.asarray(x, float)
    return 10 * x.shape[-1] + np.sum(x**2 - 10 * np.cos(2 * np.pi * x), axis=-1)


f1 = lambda x: 0.1 * x**2 + 3 * np.sin(2 * x) + 2 * np.cos(3.3 * x + 1)
f = lambda x: f1(x[0])          # f de um vetor com 1 componente (como no script)
LB, UB = -8.0, 8.0
XS = 2.4680679685               # mínimo global da função de referência

ok = True
print("Otimização — deck 5.2: simulated annealing")

# ------------------------------------------------ 1. Execução passo a passo
x0, T0, c, sigma, Nmax = 6.0, 3.0, 0.98, 0.5, 400
print("\n1. Execução passo a passo: x0 = %g, T0 = %g, c = %g, sigma = %g, Nmax = %d, semente 3"
      % (x0, T0, c, sigma, Nmax))
r = simulated_annealing(f, x0, T0, c, sigma, Nmax, LB, UB, rng=np.random.default_rng(3))
H = r.history                   # [k T aceite f(x') f(x_k) f_best x' x_k]
print("%5s %9s %9s %9s %9s" % ("k", "T", "x_k", "f(x_k)", "f_best"))
for k in (0, 30, 100, 200, 314, 315, 400):
    print("%5d %9.4f %9.4f %9.4f %9.4f" % (k, H[k, 1], H[k, 7], H[k, 4], H[k, 5]))
# saltos: propostas aceites a mais de 1 de distância do ponto corrente
xk = H[:, 7]
saltos = [k for k in range(1, Nmax + 1) if H[k, 2] == 1 and abs(xk[k] - xk[k - 1]) > 1]
kj = saltos[-1]
dz = (xk[kj] - xk[kj - 1]) / sigma
print("último salto aceite: k = %d, de %.2f para %.2f (%.1f sigma)" % (kj, xk[kj - 1], xk[kj], dz))
print("x = %.4f, f(x) = %.4f, x* = %.4f;  n_f = %d, aceites = %d, T final = %.2e"
      % (r.x[0], r.fx, XS, r.nfev, r.nacc, r.T))
c1 = [confere(H[30, 1], 1.7, 1), confere(H[100, 1], 0.4, 1), confere(H[200, 1], 0.05, 2),
      abs(xk[100] - 4.72) < 0.5 and abs(xk[200] - 4.72) < 0.5, confere(H[100, 5], 0.89, 2),
      kj == 315 and confere([xk[314], xk[315]], [4.73, 3.07], 2) and confere(abs(dz), 3.3, 1),
      confere(r.x[0], 2.464, 3), confere(XS, 2.468, 3), r.nfev == 401]
print("  T(30) = 1.7: %s | T(100) = 0.4: %s | T(200) = 0.05: %s | bacia de 4.7 em k = 100, 200: %s |"
      " f = 0.89: %s\n  salto em k = 315 (4.73 -> 3.07, 3.3 sigma): %s | x = 2.464: %s |"
      " x* = 2.468: %s | n_f = 401: %s" % simnao(c1))
ok = ok and all(c1)

# SA -> método local: pesquisa direta 1D (poll {+1, -1}), passo inicial 0.01,
# mantido no sucesso e reduzido a metade no insucesso, até ser <= 1e-6
# (a de make_figs_54.py; não conta f(x_SA), já avaliado pelo SA)
x = float(r.x[0]); fx = f1(x); nf = 0; a = 0.01
while a > 1e-6:
    sucesso = False
    for d in (1, -1):
        y = x + a * d; fy = f1(y); nf += 1
        if fy < fx:
            x, fx, sucesso = y, fy, True
            break
    if not sucesso:
        a /= 2
print("SA -> método local: de %.4f, pesquisa direta: x = %.7f, |x - x*| = %.1e, %d avaliações"
      " (total %d)" % (r.x[0], x, abs(x - XS), nf, r.nfev + nf))
c1 = [abs(x - XS) < 1e-6, nf == 37]
print("  |x - x*| < 1e-6: %s | 37 avaliações: %s" % simnao(c1))
ok = ok and all(c1)

# ------------------------------------------------ 2. 100 corridas
print("\n2. 100 corridas de x0 = 6 (semente 99), sucesso |x - x*| < 0.3")
print("%6s %6s %7s %6s %10s %8s %8s" % ("sigma", "T0", "c", "Nmax", "T_Nmax", "sucesso", "slide"))
linhas = [(0.5, 3.0, 0.98, 400, 35), (0.5, 3.0, 0.90, 400, 19), (0.5, 3.0, 0.995, 400, 74),
          (0.5, 0.3, 0.98, 400, 20), (1.0, 3.0, 0.98, 400, 88), (0.5, 3.0, 0.98, 2000, 73)]
Tslide = [9.3e-4, 1.5e-18, 4.0e-1, 9.3e-5, 9.3e-4, 8.5e-18]
c2 = []
for (sg, T0_, c_, N_, s_), Ts in zip(linhas, Tslide):
    rng = np.random.default_rng(99)
    suc = 0
    for _ in range(100):
        rr = simulated_annealing(f, 6.0, T0_, c_, sg, N_, LB, UB, rng=rng)
        suc += abs(rr.x[0] - XS) < 0.3
    print("%6.1f %6.1f %7.3f %6d %10.1e %7d%% %7d%%" % (sg, T0_, c_, N_, rr.T, suc, s_))
    c2 += [suc == s_, "%.1e" % rr.T == "%.1e" % Ts]     # T_Nmax com 2 algarismos
print("  taxas 35/19/74/20/88/73 %%: %s | T_Nmax: %s"
      % ("sim" if all(c2[0::2]) else "não", "sim" if all(c2[1::2]) else "não"))
ok = ok and all(c2)

# ------------------------------------------------ 3. Rastrigin 2D
print("\n3. Rastrigin 2D em [-5.12, 5.12]^2: T0 = 20, c = 0.995, sigma = 0.5, Nmax = 3000")
r = simulated_annealing(rastrigin, [4.0, -4.0], 20.0, 0.995, 0.5, 3000, -5.12, 5.12,
                        rng=np.random.default_rng(3))
print("de (4, -4), semente 3: x = (%.4f; %.4f), f = %.4f, |x| = %.4f, n_f = %d, aceites = %d"
      % (r.x[0], r.x[1], r.fx, np.linalg.norm(r.x), r.nfev, r.nacc))
c3 = [confere(r.fx, 0.068, 3), confere(np.linalg.norm(r.x), 0.02, 2), r.nfev == 3001]
print("  f = 0.068: %s | |x| = 0.02: %s | n_f = 3001: %s" % simnao(c3))
ok = ok and all(c3)

print("30 corridas de pontos aleatórios (semente 11), sucesso f < 1:")
c3 = []
for c_, s_ in ((0.90, 97), (0.95, 97), (0.99, 97), (0.995, 100), (0.999, 90)):
    rng = np.random.default_rng(11)
    suc = 0
    for _ in range(30):
        x0 = rng.uniform(-5.12, 5.12, 2)
        rr = simulated_annealing(rastrigin, x0, 20.0, c_, 0.5, 3000, -5.12, 5.12, rng=rng)
        suc += rr.fx < 1.0
    print("  c = %.3f: sucesso %3.0f%%  (slide %d%%), T_3000 = %.1e"
          % (c_, 100 * suc / 30, s_, rr.T))
    c3.append(round(100 * suc / 30) == s_)
print("  taxas 97/97/97/100/90 %%: %s" % ("sim" if all(c3) else "não"))
ok = ok and all(c3)

# a chamada do slide «No computador», com default_rng(s) no lugar de rng(s)
# (informativo: o slide não mostra o resultado desta chamada)
fs = np.zeros(30)
for s in range(1, 31):
    rng = np.random.default_rng(s)
    x0 = -5.12 + 10.24 * rng.random(2)
    rr = simulated_annealing(rastrigin, x0, 20, 0.995, 0.5, 3000, -5.12, 5.12, rng=rng)
    fs[s - 1] = rr.fx
q1, q3 = np.percentile(fs, [25, 75])
print("chamada do slide (sementes 1..30): mediana f = %.3f, [Q1; Q3] = [%.3f; %.3f], sucesso %.0f%%,"
      " nfev = %d  (informativo)" % (np.median(fs), q1, q3, 100 * np.mean(fs < 1), rr.nfev))

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
