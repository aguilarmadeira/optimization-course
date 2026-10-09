"""Reproduz os exemplos do deck 5.3 (algoritmos genéticos).

1. «Uma geração à mão»: min x^2 em [-5, 5], 8 bits, N = 6, roleta com
   F = 1/(1+x^2), cruzamento a 1 ponto (pc = 0.8), mutação por bit
   (pm = 0.05); semente 5 do script. População inicial, roleta (xi e pais),
   cortes, descendentes, a mutação e as médias de f: 11.9 -> 10.4 (o melhor
   piora: 0.111 -> 0.323). É o AG binário do slide, não genetic_algorithm
   (que usa codificação real); por isso está escrito aqui, passo a passo.
2. «Rastrigin: explorar ou já convergir?»: N = 40, G = 60, pc = 0.9,
   pm = 0.1, sigma = 0.5, 2 elites (2440 avaliações), semente 2: melhor e
   média nas gerações 0, 5, 20, 60, f < 1e-12 no fim. E «A mesma função,
   várias estratégias»: 30 corridas (semente 21, a amostra do script, com o
   SA do 5.2 intercalado): AG 100 % com n_f = 2440, SA 97 % (ponto final) com
   n_f mediana 5323 [4196; 6070].
   E, a título informativo, a chamada do slide «No computador» (sementes 1..30).
3. «O que fazem os parâmetros?»: função de referência do cap. 5, 300 corridas
   (semente 11), sucesso |x - x*| < 0.3: referência 75 %, s = 6 61 %, pc = 0
   60 %, N = 6 34 %, N = 60 96 %; pm = 0 76 %, pm = 0.8 78 %; mesmo orçamento
   n_f = 420: N = 60 (G = 6) 96 %, N = 6 (G = 69) 46 %.

Os números aleatórios são os dos scripts das figuras (make_figs_5.py e
make_figs_54.py): o mesmo gerador numpy (default_rng), as mesmas sementes e a
mesma ordem de consumo, por isso os valores dos slides reproduzem-se
EXATAMENTE. Uma única diferença de ordem, na parte 3, é tratada com um
adaptador explicado abaixo (OrdemMakeFigs54).
Demora cerca de 20 s (2700 corridas do AG na parte 3). Com COMPLETO = False
a parte 3 faz só 30 corridas por linha e mostra as taxas sem as conferir
(os valores exatos do slide só saem com COMPLETO = True).
Usa simulated_annealing da pasta 05_2_simulated_annealing (parte 2).
Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex05_3_genetic_algorithms.py

Otimização — deck 5.3.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from genetic_algorithm import genetic_algorithm
from simulated_annealing import simulated_annealing  # (5.2)

COMPLETO = True      # False: parte 3 com 30 corridas por linha, sem conferir as taxas


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


ok = True
print("Otimização — deck 5.3: algoritmos genéticos")

# ------------------------------------------------ 1. Uma geração à mão
print("\n1. Uma geração à mão: min x^2 em [-5, 5], 8 bits, N = 6, roleta, pc = 0.8, pm = 0.05 (semente 5)")
L, U, nb, N, pc, pm = -5.0, 5.0, 8, 6, 0.8, 0.05
decode = lambda b: L + (U - L) * int("".join(str(v) for v in b), 2) / (2**nb - 1)
cadeia = lambda b: "".join(str(v) for v in b)
rng = np.random.default_rng(5)
pop = rng.integers(0, 2, (N, nb))                      # população inicial
x = np.array([decode(b) for b in pop]); fx = x**2
Fa = 1 / (1 + fx); p = Fa / Fa.sum(); acum = np.cumsum(p)
print("%2s %9s %7s %7s %6s %6s %6s" % ("i", "cromos.", "x", "f", "F", "p_i", "acum."))
for i in range(N):
    print("%2d %9s %7.3f %7.3f %6.3f %6.3f %6.3f" % (i + 1, cadeia(pop[i]), x[i], fx[i], Fa[i], p[i], acum[i]))
xi = rng.random(N); pais = np.searchsorted(acum, xi)  # roleta: 1.º i com acumulada >= xi
print("roleta: xi =", " ".join("%.3f" % v for v in xi), "-> pais", " ".join(str(v + 1) for v in pais))
rc = rng.random(N // 2); cortes = rng.integers(1, nb, N // 2)
filhos = []
for j in range(N // 2):
    a, b = pop[pais[2 * j]], pop[pais[2 * j + 1]]
    if rc[j] < pc:                                     # cruzamento a 1 ponto
        cj = cortes[j]
        c1 = np.concatenate([a[:cj], b[cj:]]); c2 = np.concatenate([b[:cj], a[cj:]])
        nota = "corte após bit %d" % cj
    else:
        c1, c2 = a.copy(), b.copy(); nota = "sem cruzamento"
    filhos += [c1, c2]
    print("pais %d, %d: %s, %s  %s  -> %s, %s" % (pais[2 * j] + 1, pais[2 * j + 1] + 1, cadeia(a),
                                                cadeia(b), nota, cadeia(c1), cadeia(c2)))
filhos = np.array(filhos)
mut = rng.random(filhos.shape) < pm                    # mutação: inverter o bit
filhos_m = np.where(mut, 1 - filhos, filhos)
x1 = np.array([decode(b) for b in filhos_m]); fx1 = x1**2
print("mutações: %d (descendente %s, bit %s)" % (mut.sum(), ", ".join(str(i + 1) for i in np.argwhere(mut)[:, 0]),
                                                ", ".join(str(j + 1) for j in np.argwhere(mut)[:, 1])))
for i in range(N):
    print("%2d %9s -> %9s %7.3f %7.3f" % (i + 1, cadeia(filhos[i]), cadeia(filhos_m[i]), x1[i], fx1[i]))
print("média de f: %.1f -> %.1f;  melhor: %.3f -> %.3f" % (fx.mean(), fx1.mean(), fx.min(), fx1.min()))
c1 = [[cadeia(b) for b in pop] == ["11010110", "10001000", "00010110", "00010111", "00101110", "11100100"],
      confere(x, [3.392, 0.333, -4.137, -4.098, -3.196, 3.941], 3),
      confere(fx, [11.507, 0.111, 17.117, 16.794, 10.215, 15.533], 3),
      confere(Fa, [0.080, 0.900, 0.055, 0.056, 0.089, 0.060], 3),
      confere(p, [0.064, 0.725, 0.044, 0.045, 0.072, 0.049], 3),
      confere(acum, [0.064, 0.790, 0.834, 0.879, 0.951, 1.000], 3),
      confere(xi, [0.679, 0.870, 0.227, 0.895, 0.872, 0.019], 3) and list(pais + 1) == [2, 4, 2, 5, 4, 1],
      list(cortes) == [3, 4, 7] and all(rc < pc),
      [cadeia(b) for b in filhos_m] == ["10010111", "00001100", "10001110", "00101000", "00010110", "11010111"]
      and mut.sum() == 1,
      confere(x1, [0.922, -4.529, 0.569, -3.431, -4.137, 3.431], 3),
      confere(fx1, [0.849, 20.516, 0.323, 11.774, 17.117, 11.774], 3),
      confere([fx.mean(), fx1.mean()], [11.9, 10.4], 1) and confere([fx.min(), fx1.min()], [0.111, 0.323], 3)]
print("  cromossomas: %s | x: %s | f: %s | F: %s | p_i: %s | acumulada: %s | roleta: %s |\n"
      "  cortes 3, 4, 7: %s | descendentes e mutação: %s | x: %s | f: %s | médias 11.9 -> 10.4: %s" % simnao(c1))
ok = ok and all(c1)

# ------------------------------------------------ 2. Rastrigin 2D
print("\n2. Rastrigin 2D em [-5.12, 5.12]^2: N = 40, G = 60, pc = 0.9, pm = 0.1, sigma = 0.5, 2 elites")
r = genetic_algorithm(rastrigin, 2, -5.12, 5.12, 40, 60, 0.9, 0.1, 0.5, 2, rng=np.random.default_rng(2))
H = r.history                                          # [t melhor média x]
for t in (0, 5, 20, 60):
    print("  geração %2d: melhor f = %.3g, média = %.1f" % (t, H[t, 1], H[t, 2]))
print("semente 2: x = (%.1e; %.1e), f = %.2e, n_f = %d" % (r.x[0], r.x[1], r.fx, r.nfev))
c2 = [H[20, 1] < 1, r.fx < 1e-12, r.nfev == 2440,
      confere(H[[0, 5, 20, 60], 2], [31.3, 15.2, 2.0, 1.4], 1) and confere(H[[0, 5], 1], [6.479, 0.581], 3)]
print("  na bacia do global na geração 20: %s | f < 1e-12: %s | n_f = 2440: %s |"
      " melhor/média da figura: %s" % simnao(c2))
ok = ok and all(c2)

# a amostra de 30 pontos do script (semente 21): SA e AG intercalados no mesmo gerador
rng = np.random.default_rng(21)
fsa = np.zeros(30); fag = np.zeros(30); nsa = np.zeros(30)
for k in range(30):
    x0 = rng.uniform(-5.12, 5.12, 2)
    rs = simulated_annealing(rastrigin, x0, 20.0, 0.5, 10, 0.1, 0.5, -5.12, 5.12, rng=rng)   # SA das aulas
    ra = genetic_algorithm(rastrigin, 2, -5.12, 5.12, 40, 60, 0.9, 0.1, 0.5, 2, rng=rng)
    fsa[k] = rs.fx; fag[k] = ra.fx; nsa[k] = rs.nfev
q1, q2, q3 = np.percentile(nsa, [25, 50, 75])
print("30 corridas (semente 21), sucesso f < 1: AG %.0f%% (n_f = %d), SA %.0f%% (ponto final; n_f mediana"
      " %.0f [%.0f; %.0f]); medianas f: AG %.1e, SA %.3f"
      % (100 * np.mean(fag < 1), ra.nfev, 100 * np.mean(fsa < 1), q2, q1, q3, np.median(fag), np.median(fsa)))
c2 = [np.all(fag < 1), ra.nfev == 2440, round(100 * np.mean(fsa < 1)) == 97,
      ["%.0f" % v for v in (q2, q1, q3)] == ["5323", "4196", "6070"]]
print("  AG 100 %%: %s | 2440: %s | SA 97 %%: %s | n_f do SA 5323 [4196; 6070]: %s" % simnao(c2))
ok = ok and all(c2)

# a chamada do slide «No computador», com default_rng(s) no lugar de rng(s) (informativo)
fs = np.zeros(30)
for s in range(1, 31):
    rr = genetic_algorithm(rastrigin, 2, -5.12, 5.12, 40, 60, 0.9, 0.1, 0.5, 2,
                           rng=np.random.default_rng(s))
    fs[s - 1] = rr.fx
q1, q3 = np.percentile(fs, [25, 75])
print("chamada do slide (sementes 1..30): mediana f = %.1e, [Q1; Q3] = [%.1e; %.1e], sucesso %.0f%%,"
      " nfev = %d  (informativo)" % (np.median(fs), q1, q3, 100 * np.mean(fs < 1), rr.nfev))

# ------------------------------------------------ 3. O que fazem os parâmetros?
f1 = lambda x: 0.1 * x**2 + 3 * np.sin(2 * x) + 2 * np.cos(3.3 * x + 1)
# f de um vetor com 1 componente, avaliada como vetor (o script avalia a
# população de uma vez; assim os valores coincidem bit a bit)
f = lambda x: f1(x)[0]
LB, UB = -8.0, 8.0
XS = 2.4680679685


class OrdemMakeFigs54:
    """Adaptador do gerador, só para reproduzir make_figs_54.py bit a bit.

    genetic_algorithm sorteia os candidatos do torneio como uma matriz s x N
    (linha r = r-ésimo candidato de cada pai: i1, i2 do slide); make_figs_54.py
    sorteia-os como N x s (linha = pai). A distribuição é a mesma; só muda a
    ordem em que os números saem do gerador. Este adaptador faz o sorteio na
    ordem do script e devolve-o transposto; tudo o resto passa ao gerador.
    """
    def __init__(self, rng):
        self.rng = rng

    def integers(self, low, high, size):
        s, N = size
        return self.rng.integers(low, high, (N, s)).T

    def __getattr__(self, nome):
        return getattr(self.rng, nome)


def taxa(R, N=20, G=20, pc=0.9, pm=0.1, sigma=0.5, s=2):
    rng = OrdemMakeFigs54(np.random.default_rng(11))
    suc = 0
    for _ in range(R):
        rr = genetic_algorithm(f, 1, LB, UB, N, G, pc, pm, sigma, 2, s=s, rng=rng)
        suc += abs(rr.x[0] - XS) < 0.3
    return 100 * suc / R, rr.nfev


R = 300 if COMPLETO else 30
print("\n3. Função de referência, %d corridas (semente 11), sucesso |x - x*| < 0.3%s"
      % (R, "" if COMPLETO else "  [COMPLETO = False: taxas só indicativas]"))
linhas = [("referência: N = 20, s = 2, pm = 0.1", {}, 75, 420),
          ("pressão forte: s = 6", dict(s=6), 61, 420),
          ("sem cruzamento: pc = 0", dict(pc=0.0), 60, 420),
          ("população pequena: N = 6", dict(N=6), 34, 126),
          ("população maior: N = 60", dict(N=60), 96, 1260),
          ("sem mutação: pm = 0", dict(pm=0.0), 76, 420),
          ("mutação alta: pm = 0.8", dict(pm=0.8), 78, 420),
          ("mesmo orçamento: N = 60, G = 6", dict(N=60, G=6), 96, 420),
          ("mesmo orçamento: N = 6, G = 69", dict(N=6, G=69), 46, 420)]
c3t = []; c3n = []
for nome, kw, s_, nf_ in linhas:
    tx, nf = taxa(R, **kw)
    print("  %-36s sucesso %5.1f%%  (slide %d%%)  n_f = %4d (slide %d)" % (nome, tx, s_, nf, nf_))
    c3t.append(round(tx) == s_); c3n.append(nf == nf_)
if COMPLETO:
    print("  taxas 75/61/60/34/96, 76/78, 96/46 %%: %s | n_f: %s" % simnao([all(c3t), all(c3n)]))
    ok = ok and all(c3t) and all(c3n)
else:
    print("  n_f: %s (taxas não conferidas: só com COMPLETO = True, 300 corridas)" % simnao([all(c3n)]))
    ok = ok and all(c3n)

print("\nconfere com os slides: %s%s" % ("sim" if ok else "não",
                                        "" if COMPLETO else " (sem as taxas da parte 3)"))
