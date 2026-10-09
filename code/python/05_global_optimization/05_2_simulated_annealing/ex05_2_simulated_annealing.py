"""Reproduz os exemplos do deck 5.2 (simulated annealing, versão das aulas).

Função de referência do cap. 5: f(x) = 0.1 x^2 + 3 sin 2x + 2 cos(3.3x + 1)
em X = [-8, 8], mínimo global x* = 2.4681 (f = -4.2381).

1. «Uma execução passo a passo»: x0 = 6, T0 = 3, alpha = 0.8, n = 10,
   Ts = 0.05, sigma = 0.5, semente 3: T e x_t em n_f = 60, 150, 300; o salto
   em n_f = 325 (de 4.59 para 2.92, a 3.3 sigma, T = 0.084); termina com
   T = 0.043 < Ts em x_t = 2.457, n_f = 507, 191 aceites.
   E «SA -> método local»: de x_t, uma pesquisa direta chega a |x - x*| < 1e-6
   com 38 avaliações.
2. «Uma execução não chega: 100 corridas»: a tabela de 6 linhas (sempre de
   x0 = 6, semente 99): sucesso |x - x*| < 0.3 do ponto final e do melhor
   visitado, e a mediana de n_f.
3. «Em duas dimensões: Rastrigin»: de (4, -4), T0 = 20, alpha = 0.5, n = 10,
   Ts = 0.1, sigma = 0.5, semente 3: f = 0.18, n_f = 5421, 81 aceites;
   30 corridas de pontos aleatórios (semente 11), sucesso f < 1: 90 %
   (melhor visitado 100 %), n_f ~ 4800; alpha = 0.8: 97 %, n_f ~ 18 600;
   alpha = 0.3: 97 %. E, a título informativo, a chamada do slide
   «No computador» (sementes 1..30).

Os números aleatórios são os dos scripts das figuras (make_figs_5.py): o
mesmo gerador numpy (default_rng), as mesmas sementes e a mesma ordem de
consumo, por isso os valores dos slides reproduzem-se EXATAMENTE.
Demora cerca de um minuto (cerca de 2 milhões de avaliações ao todo).
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
print("Otimização — deck 5.2: simulated annealing (versão das aulas)")

# ------------------------------------------------ 1. Execução passo a passo
x0, T0, alpha, n, Ts, sigma = 6.0, 3.0, 0.8, 10, 0.05, 0.5
print("\n1. Execução passo a passo: x0 = %g, T0 = %g, alpha = %g, n = %d, Ts = %g, sigma = %g, semente 3"
      % (x0, T0, alpha, n, Ts, sigma))
r = simulated_annealing(f, x0, T0, alpha, n, Ts, sigma, LB, UB, rng=np.random.default_rng(3))
H = r.history                   # [n_f t T aceite f(x') f(x_t) f_best x' x_t]
print("%6s %5s %8s %9s %9s %9s" % ("n_f", "t", "T", "x_t", "f(x_t)", "f_best"))
for k in (1, 60, 150, 300, 324, 325, r.nfev):
    i = k - 1
    print("%6d %5d %8.4f %9.4f %9.4f %9.4f" % (H[i, 0], H[i, 1], H[i, 2], H[i, 8], H[i, 5], H[i, 6]))
xt = H[:, 8]; ac = H[:, 3] == 1
saltos = [i for i in range(1, len(H)) if ac[i] and xt[i - 1] > 3.68 and xt[i] < 3.68]
ij = saltos[-1]
dz = (xt[ij] - xt[ij - 1]) / sigma
print("salto: n_f = %d, de %.2f para %.2f (%.1f sigma), T = %.3f" % (H[ij, 0], xt[ij - 1], xt[ij], dz, H[ij, 2]))
print("ponto final x_t = %.4f (f = %.4f), melhor visitado %.4f;  n_f = %d, aceites = %d, T final = %.3f"
      % (r.x[0], r.fx, r.xbest[0], r.nfev, r.nacc, r.T))
c1 = [confere([H[59, 2], H[149, 2], H[299, 2]], [1.2, 0.32, 0.11], [1, 2, 2]),
      all(abs(xt[[59, 149, 299]] - 4.72) < 0.5), confere(np.min(f1(np.linspace(4.5, 5.0, 5001))), 0.89, 2),
      len(saltos) == 1 and H[ij, 0] == 325 and confere([xt[ij - 1], xt[ij]], [4.59, 2.92], 2)
      and confere(abs(dz), 3.3, 1) and confere(H[ij, 2], 0.084, 3),
      r.flag == 0 and confere(r.T, 0.043, 3), confere(r.x[0], 2.457, 3), r.nfev == 507, r.nacc == 191]
print("  T(60, 150, 300) = 1.2, 0.32, 0.11: %s | bacia de 4.7 em n_f = 60, 150, 300: %s | mínimo local f = 0.89: %s |\n"
      "  salto em n_f = 325 (4.59 -> 2.92, 3.3 sigma, T = 0.084): %s | T final 0.043 < Ts: %s |"
      " x_t = 2.457: %s | n_f = 507: %s | 191 aceites: %s" % simnao(c1))
ok = ok and all(c1)

# SA -> método local: pesquisa direta 1D (poll {+1, -1}), passo inicial 0.01,
# mantido no sucesso e reduzido a metade no insucesso, até ser <= 1e-6
# (não conta f(x_SA), já avaliado pelo SA)
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
c1 = [abs(x - XS) < 1e-6, nf == 38]
print("  |x - x*| < 1e-6: %s | 38 avaliações: %s" % simnao(c1))
ok = ok and all(c1)

# ------------------------------------------------ 2. 100 corridas
print("\n2. 100 corridas de x0 = 6 (semente 99), Ts = 0.05, sucesso |x - x*| < 0.3")
print("%6s %5s %6s %4s %22s %8s %8s %16s" % ("sigma", "T0", "alpha", "n", "n_f med [Q1; Q3]", "final", "melhor", "slide"))
linhas = [(0.5, 3.0, 0.80, 10, 693, 60, 63), (0.5, 3.0, 0.50, 10, 148, 11, 11), (0.5, 3.0, 0.95, 10, 3318, 69, 93),
          (0.5, 0.3, 0.80, 10, 404, 15, 15), (1.0, 3.0, 0.80, 10, 1518, 60, 92), (0.5, 3.0, 0.80, 50, 3942, 64, 93)]
c2 = []
for (sg, T0_, a_, n_, nf_s, s_, sb_) in linhas:
    rng = np.random.default_rng(99)
    suc = sucb = 0; nfs = []
    for _ in range(100):
        rr = simulated_annealing(f, 6.0, T0_, a_, n_, 0.05, sg, LB, UB, rng=rng)
        suc += abs(rr.x[0] - XS) < 0.3; sucb += abs(rr.xbest[0] - XS) < 0.3; nfs.append(rr.nfev)
    q1, q2, q3 = np.percentile(nfs, [25, 50, 75])
    print("%6.1f %5.1f %6.2f %4d %8.1f [%5.0f; %5.0f] %7d%% %7d%%   %d%%, %d%%, %d"
          % (sg, T0_, a_, n_, q2, q1, q3, suc, sucb, s_, sb_, nf_s))
    c2 += [suc == s_ and sucb == sb_, "%.0f" % q2 == "%d" % nf_s]
print("  taxas (ponto final / melhor): %s | medianas de n_f: %s"
      % ("sim" if all(c2[0::2]) else "não", "sim" if all(c2[1::2]) else "não"))
ok = ok and all(c2)

# ------------------------------------------------ 3. Rastrigin 2D
print("\n3. Rastrigin 2D em [-5.12, 5.12]^2: T0 = 20, alpha = 0.5, n = 10, Ts = 0.1, sigma = 0.5")
r = simulated_annealing(rastrigin, [4.0, -4.0], 20.0, 0.5, 10, 0.1, 0.5, -5.12, 5.12,
                        rng=np.random.default_rng(3))
print("de (4, -4), semente 3: x_t = (%.4f; %.4f), f = %.4f, n_f = %d, aceites = %d; melhor visitado f = %.4f"
      % (r.x[0], r.x[1], r.fx, r.nfev, r.nacc, r.fbest))
c3 = [confere(r.fx, 0.18, 2), r.nfev == 5421, r.nacc == 81]
print("  f = 0.18: %s | n_f = 5421: %s | 81 aceites: %s" % simnao(c3))
ok = ok and all(c3)

print("30 corridas de pontos aleatórios (semente 11), sucesso f < 1 (ponto final / melhor visitado):")
c3 = []
for a_, s_, sb_, nf_s in ((0.3, 97, 97, 3630), (0.5, 90, 100, 4814), (0.8, 97, 100, 18572)):
    rng = np.random.default_rng(11)
    suc = sucb = 0; nfs = []
    for _ in range(30):
        x0 = rng.uniform(-5.12, 5.12, 2)
        rr = simulated_annealing(rastrigin, x0, 20.0, a_, 10, 0.1, 0.5, -5.12, 5.12, rng=rng)
        suc += rr.fx < 1.0; sucb += rr.fbest < 1.0; nfs.append(rr.nfev)
    q2 = np.median(nfs)
    print("  alpha = %.1f: sucesso %3.0f%% / %3.0f%%, n_f mediana %.1f  (slide %d%% / %d%%, %d)"
          % (a_, 100 * suc / 30, 100 * sucb / 30, q2, s_, sb_, nf_s))
    c3.append(round(100 * suc / 30) == s_ and round(100 * sucb / 30) == sb_ and abs(q2 - nf_s) <= 0.5)
print("  alpha = 0.3, 0.5, 0.8: %s" % ("sim" if all(c3) else "não"))
ok = ok and all(c3)

# a chamada do slide «No computador», com default_rng(s) no lugar de rng(s)
# (informativo: o slide não mostra o resultado desta chamada)
fs = np.zeros(30); nfs = np.zeros(30)
for s in range(1, 31):
    rng = np.random.default_rng(s)
    x0 = -5.12 + 10.24 * rng.random(2)
    rr = simulated_annealing(rastrigin, x0, 20, 0.5, 10, 0.1, 0.5, -5.12, 5.12, rng=rng)
    fs[s - 1] = rr.fx; nfs[s - 1] = rr.nfev
q1, q3 = np.percentile(fs, [25, 75])
print("chamada do slide (sementes 1..30): mediana f = %.3f, [Q1; Q3] = [%.3f; %.3f], sucesso %.0f%%,"
      " n_f mediana %.0f  (informativo)" % (np.median(fs), q1, q3, 100 * np.mean(fs < 1), np.median(nfs)))

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
