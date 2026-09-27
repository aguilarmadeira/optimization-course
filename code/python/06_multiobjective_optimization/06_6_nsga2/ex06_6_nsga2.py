"""Reproduz os exemplos do deck 6.6 (NSGA-II) e o hipervolume do deck 6.2.

1. Oito soluções (frames «rank» e «crowding»): F_1 = {1,2,4,8},
   F_2 = {3,5,6,7}; d_2 = 1.4, d_4 = 1.233; só cabem 3: sai o 4.
2. Hipervolume (deck 6.2): três pontos com r = (1.1; 1.3): 0.549;
   convergência/diversidade com r = (1.1; 1.65): 0.654, 0.656, 1.095;
   frentes exatas com r = (1.1; 1.65): convexa 1.648, não convexa 1.148.
3. NSGA-II no exemplo convexo da UC, N = 40, G = 50 (semente 4):
   geração 0 com 6 não dominados, 2040 avaliações, HV 1.637.
4. Tabela «DMS vs. NSGA-II» (exemplo não convexo): coluna do NSGA-II,
   20 corridas (sementes 100-119), mediana e quartis do HV por orçamento.

Método estocástico: usa o gerador numpy.random.default_rng e as mesmas
sementes do script das figuras do capítulo 6, por isso reproduz EXATAMENTE
os números dos slides.  Os slides usam vírgula decimal; aqui usa-se o ponto.
Correr (de qualquer pasta):  python ex06_6_nsga2.py   (cerca de 7 s)

Otimização — deck 6.6.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from non_dominated_sort import non_dominated_sort
from crowding_distance import crowding_distance
from hypervolume_2d import hypervolume_2d
from nsga2 import nsga2


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


def um(fr):                              # índices base 0 -> numeração do slide
    return "{" + ",".join(str(i + 1) for i in sorted(fr)) + "}"


ok = True
print("Otimização — deck 6.6: NSGA-II (e hipervolume, deck 6.2)")

# ------------------------------------------------ 1. oito soluções
Pex = np.array([[1, 6], [2, 3], [3, 5], [4, 1.5], [5, 4], [6, 2.5], [2.5, 6.5], [7, 1]])
print("\n1. Oito soluções (f1, f2): rank e crowding")
fronts, rank = non_dominated_sort(Pex)
d = crowding_distance(Pex, fronts[0])
dp = dict(zip(fronts[0], d))
print("%3s %12s %5s %9s" % ("p", "(f1; f2)", "rank", "d_p"))
for i in range(len(Pex)):
    print("%3d %12s %5d %9s" % (i + 1, "(%g; %g)" % tuple(Pex[i]), rank[i],
                                "%.3f" % dp[i] if i in dp else "--"))
print("F_1 = %s, F_2 = %s" % (um(fronts[0]), um(fronts[1])))
fica = sorted(range(len(d)), key=lambda k: -d[k])[:3]          # só cabem 3 dos 4
fica = [fronts[0][k] for k in fica]
sai = [i for i in fronts[0] if i not in fica]
print("Só cabem 3 dos 4 de F_1: ficam %s, sai %s" % (um(fica), um(sai)))
c = [fronts[0] == [0, 1, 3, 7] and sorted(fronts[1]) == [2, 4, 5, 6],
     confere(dp[1], 1.4, 1), confere(dp[3], 1.233, 3),
     bool(np.isinf(dp[0]) and np.isinf(dp[7])), sorted(fica) == [0, 1, 7] and sai == [3]]
print("  frentes: %s | d_2 = 1.4: %s | d_4 = 1.233: %s | extremos d = inf: %s | sai o 4: %s"
      % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 2. hipervolume (6.2)
print("\n2. Hipervolume (deck 6.2)")
P3 = np.array([[0.2, 0.96], [0.5, 0.75], [0.8, 0.36]])
hv3 = hypervolume_2d(P3, (1.1, 1.3))
REF = (1.1, 1.65)
ff = lambda s: np.stack([s, 1 - s**2], -1)          # frente não convexa f2 = 1 - f1^2
A = ff(np.linspace(0.40, 0.55, 10))                  # perto, concentrados
B = ff(np.linspace(0, 1, 10)) + np.array([0.12, 0.25])  # espalhados, longe
C = ff(np.linspace(0, 1, 10))                        # perto e espalhados
hvA, hvB, hvC = (hypervolume_2d(Z, REF) for Z in (A, B, C))
t = np.linspace(0, 1, 2000)
hv_cx = hypervolume_2d(np.stack([t**2, (1 - t)**2], -1), REF)   # frente convexa, 2000 pontos
hv_nc = hypervolume_2d(ff(t), REF)
hv_cx_ex = 1.1 * 1.65 - 1 / 6        # exato: r1*r2 - integral_0^1 (1 - sqrt(f1))^2 df1
hv_nc_ex = 1.1 * 1.65 - 2 / 3        # exato: r1*r2 - integral_0^1 (1 - f1^2) df1
print("três pontos, r = (1.1; 1.3): HV = %.3f" % hv3)
print("r = (1.1; 1.65): perto/concentrados %.3f, espalhados/longe %.3f, perto/espalhados %.3f"
      % (hvA, hvB, hvC))
print("frente exata convexa: %.4f (2000 pontos: %.4f); não convexa: %.4f (2000 pontos: %.4f)"
      % (hv_cx_ex, hv_cx, hv_nc_ex, hv_nc))
c = [confere(hv3, 0.549, 3), confere(hvA, 0.654, 3), confere(hvB, 0.656, 3), confere(hvC, 1.095, 3),
     confere([hv_cx_ex, hv_cx], 1.648, 3), confere([hv_nc_ex, hv_nc], 1.148, 3)]
print("  0.549: %s | 0.654: %s | 0.656: %s | 1.095: %s | 1.648: %s | 1.148: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 3. NSGA-II, exemplo convexo
f = lambda x: np.array([x[0]**2 + x[1]**2, (x[0] - 1)**2 + x[1]**2])
cv = lambda x: max(0.0, x[0]**2 + x[1]**2 - 1.0)    # x1^2 + x2^2 <= 1 (regra de Deb)
print("\n3. NSGA-II no exemplo convexo: N = 40, G = 50, semente 4, r = (1.1; 1.65)")
r = nsga2(f, [0, 0], [1, 1], 40, 50, 4, cv=cv, ref=REF, verbose=True)
H = r.history
print("geração 0: %d não dominados (%d inadmissíveis); geração 50: %d pontos, "
      "n_f = %d, HV = %.4f (exata %.4f)" % (H[0, 2], H[0, 3], len(r.F), r.nfev, H[-1, 4], hv_cx_ex))
print(r.message)
c = [H[0, 2] == 6, r.nfev == 2040, confere(H[-1, 4], 1.637, 3)]
print("  6 não dominados na geração 0: %s | n_f = 2040: %s | HV = 1.637: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------ 4. DMS vs. NSGA-II
fn = lambda x: np.array([x[0], 1 - x[0]**2 + x[1]])  # exemplo não convexo, x em [0,1]x[0,0.6]
print("\n4. DMS vs. NSGA-II no exemplo não convexo (r = (1.1; 1.65)), coluna do NSGA-II")
r = nsga2(fn, [0, 0], [1, 0.6], 40, 50, 5, ref=REF)
print("uma corrida (semente 5, G = 50): %d pontos, n_f = %d, HV = %.4f"
      % (len(r.F), r.nfev, r.history[-1, 4]))
c = [len(r.F) == 40, r.nfev == 2040]
BUDGET = 2000
HV = []
for s in range(100, 120):                            # 20 corridas, G = 2000/40 - 1 = 49
    rs = nsga2(fn, [0, 0], [1, 0.6], 40, BUDGET // 40 - 1, s, ref=REF)
    HV.append(rs.history[:, 4])
ev = rs.history[:, 1]                                # n_f de cada geração: 40, 80, ..., 2000
q1, med, q3 = np.percentile(np.array(HV), [25, 50, 75], axis=0)
orc = [100, 200, 400, 800, 2000]
dms = [1.064, 1.094, 1.123, 1.139, 1.144]            # coluna do DMS (slide; ver 06_5_dms)
Tn = [[0.986, 0.968, 1.006], [1.083, 1.066, 1.096], [1.125, 1.117, 1.129],
      [1.132, 1.131, 1.133], [1.132, 1.132, 1.133]]  # slide: mediana [Q1; Q3]
print("%8s %8s %10s %18s %24s" % ("n_f", "DMS", "NSGA med.", "[Q1; Q3]", "slide: med. [Q1; Q3]"))
obt = []
for e, dm, T in zip(orc, dms, Tn):
    k = np.searchsorted(ev, e, side="right") - 1     # última geração completa que não excede
    obt.append([med[k], q1[k], q3[k]])
    print("%8d %8.3f %10.3f   [%.3f; %.3f] %12.3f [%.3f; %.3f]" % (e, dm, med[k], q1[k], q3[k], *T))
c.append(confere(obt, Tn, 3))
print("  uma corrida: 40 pontos: %s | 2040 aval.: %s | tabela NSGA-II (20 corridas): %s" % simnao(c))
ok = ok and all(c)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
