"""Reproduz os números do deck 4.3.3 (restrições na prática).

Um só problema de referência, resolvido pelo MESMO Hooke–Jeeves da UC (3.2.5),
com a restrição comunicada ao algoritmo de cinco formas:

    min f(x) = (x1 - 2)^2 + (x2 - 1)^2   s.a.  g(x) = 2 - x1 - x2 >= 0,
    x* = (1.5, 0.5), f* = 0.5, u* = 1;  (2, 1) é inadmissível.
Convenção do capítulo 4: g >= 0, <g> = min(g, 0); violação v = -<g> = max(0, -g).

  A  ignorar a restrição:         f
  B  penalização quadrática:      f + R <g>^2 = f + R v^2,  R = 1, 10, 100, 1000
  C  penalização exata (L1):      f + R |<g>| = f + R v,  R = 0.5, 1, 2
  D  barreira extrema:            f se g >= 0, senão +inf (f não é avaliada)
  E  regra de admissibilidade:    compara (v, f), v = max(0, -g); nenhuma soma
     (hooke_jeeves com o argumento opcional menor)

Hooke–Jeeves: x0 = (0, 0), a = 2, P0 = 0.5, T = 1e-6 (todos os casos).
Também: (i) B analítico, v(x_R) = 1/(1 + 2R); (ii) restrição oculta: um
«simulador» que falha (NaN) sem dar g, tratado com +inf = D; (iii) E e D a
partir de x0 = (3, 3), inadmissível; (iv) o bloqueio em (1, 1) e a correção
com as direções rodadas (HJ nas variáveis y, x = Q y, Q = [1 1; 1 -1]/sqrt 2);
(v) informativo: o Nelder–Mead do 3.2.3 em C (R = 2) e D.
Funções: extreme_barrier.py (D e a restrição oculta); hooke_jeeves (3.2.5) e
nelder_mead (3.2.3), postas no caminho por uc_setup.
Os números dos slides saem daqui; os slides usam vírgula decimal.
Correr (de qualquer pasta):  python ex04_3_3_constraint_handling.py

Otimização — deck 4.3.3.  J. F. A. Madeira — Licença MIT.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)

import numpy as np

from hooke_jeeves import hooke_jeeves      # (3.2.5)
from nelder_mead import nelder_mead        # (3.2.3)
from extreme_barrier import extreme_barrier


def confere(v, s, d):
    """v confere com o valor s do slide, mostrado com d casas decimais?"""
    v = np.ravel(np.asarray(v, float)); s = np.ravel(np.asarray(s, float))
    d = np.ravel(np.asarray(d, float))
    return bool(np.all(np.abs(v - s) <= 0.5 * 10.0 ** (-d) * (1 + 1e-6)))


def simnao(c):
    return tuple("sim" if ci else "não" for ci in c)


# --- o problema (o mesmo em todo o deck) -----------------------------------
f = lambda x: (x[0] - 2)**2 + (x[1] - 1)**2
g = lambda x: 2 - x[0] - x[1]                      # g(x) >= 0 (convenção do cap. 4)
v = lambda x: max(0.0, -g(x))                      # violação v = -<g>
xs = np.array([1.5, 0.5])

# --- os cinco tratamentos: só muda a função (ou a comparação) --------------
FA = f                                                   # A: ignorar
FB = lambda R: (lambda x: f(x) + R * v(x)**2)            # B: quadrática
FC = lambda R: (lambda x: f(x) + R * v(x))               # C: exata (L1)
FD = lambda x: extreme_barrier(f, g, x)                  # D: barreira extrema
FE = lambda x: (v(x), f(x))                              # E: o par (v, f)
menorE = lambda a, b: a[0] < b[0] or (a[0] == b[0] and a[1] < b[1])   # E: regra

x0 = [0.0, 0.0]; a = 2; P0 = 0.5; T = 1e-6
HJ = lambda F, x0=x0, **k: hooke_jeeves(F, x0, a, P0, T, **k)

ok = True
print("Otimização — deck 4.3.3: restrições na prática (um problema, cinco tratamentos)")
print("min (x1-2)^2 + (x2-1)^2 s.a. g = 2 - x1 - x2 >= 0;  x* = (1.5, 0.5), f* = 0.5, u* = 1")
print("Hooke–Jeeves (3.2.5): x0 = (0, 0), a = 2, P0 = 0.5, T = 1e-6\n")


def linha(nome, r):
    x = r.x
    adm = "sim" if g(x) >= -1e-12 else "não"
    print("  %-26s x = (%.6f; %.6f)  f = %.6f  g = %+.2e  admissível: %-3s  n_f = %4d"
          % (nome, x[0], x[1], f(x), g(x), adm, r.nfev))
    return x


# ------------------------------------------------------------ 1. a tabela
print("1. A tabela final (slide «O mesmo problema, cinco tratamentos»)")
rA = HJ(FA); xA = linha("A ignorar", rA)
rB = {R: HJ(FB(R)) for R in (1, 10, 100, 1000)}
xB = {R: linha("B quadrática, R = %g" % R, rB[R]) for R in rB}
rC = {R: HJ(FC(R)) for R in (0.5, 1, 2)}
xC = {R: linha("C exata, R = %g" % R, rC[R]) for R in rC}
rD = HJ(FD); xD = linha("D barreira extrema", rD)
rE = HJ(FE, menor=menorE); xE = linha("E regra (v, f)", rE)

c = [confere(xA, [2, 1], 4),
     confere([xB[R] for R in (1, 10, 100, 1000)],
             [[1.667, 0.667], [1.524, 0.524], [1.502, 0.502], [1.500, 0.500]], 3),
     confere([v(xB[R]) for R in (1, 10, 100, 1000)], [0.33, 0.048, 0.0050, 0.00050], [2, 3, 4, 5]),
     confere([xC[R] for R in (0.5, 1, 2)], [[1.75, 0.75], [1.5, 0.5], [1, 1]], 4),
     confere([xD, xE], [[1, 1], [1, 1]], 4)]
print("  A = (2; 1) | B -> x* por fora, v = 0.33 ... 0.0005 | C: (1.75; 0.75), x*, (1; 1)"
      " | D = E = (1; 1): %s" % (simnao(c),))
ok = ok and all(c)
nfora = int(np.sum(np.isinf(rD.ftrace)))
c = [nfora == 45]
print("  D: n_F = %d chamadas, %d fora de X (f não avaliada): %s" % (rD.nfev, nfora, simnao(c)[0]))
ok = ok and all(c)
c = [np.array_equal(rD.history[:, 2:4], rE.history[:, 2:4]) and rD.nfev == rE.nfev]
print("  de x0 admissível, E passa pelos mesmos pontos base que D, com o mesmo n_F: %s" % simnao(c))
print("  (as tentativas em torno de pontos inadmissíveis diferem: D compara +inf com +inf, E compara v)")
ok = ok and all(c)

# ------------------------------------------------------------ 2. B analítico
print("\n2. B: x_R = (2, 1) - R v_R (1, 1), v_R = 1/(1 + 2R) (slide «Quanto deve valer R?»)")
print("  %6s %22s %10s %10s %12s %6s" % ("R", "x_R (HJ)", "v (HJ)", "1/(1+2R)", "||x_R-x*||", "n_f"))
for R in rB:
    x = xB[R]
    print("  %6g   (%.4f; %.4f) %10.5f %10.5f %12.5f %6d"
          % (R, x[0], x[1], v(x), 1 / (1 + 2 * R), np.linalg.norm(x - xs), rB[R].nfev))
c = [confere([v(xB[R]) for R in rB], [1 / (1 + 2 * R) for R in rB], 5)]
print("  HJ confere com 1/(1+2R) a 5 casas: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 3. restrição oculta
print("\n3. Restrição oculta: o «simulador» falha (NaN) e não diz porquê (não há g)")


def simulador(x):
    """Caixa negra: devolve um número, ou NaN quando «a simulação não converge»
    (aqui, quando x1 + x2 > 2), sem dizer porquê: quem a usa não conhece g."""
    return np.nan if x[0] + x[1] > 2 else (x[0] - 2)**2 + (x[1] - 1)**2


FH = lambda x: extreme_barrier(simulador, None, x)       # falhou -> +inf
rH = HJ(FH); xH = linha("oculta, +inf se falha", rH)
c = [np.array_equal(rH.history, rD.history, equal_nan=True)]
print("  a mesma corrida que D, sem nunca conhecer g: %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 4. x0 inadmissível
print("\n4. De x0 = (3, 3), inadmissível (g = -4)")
rD3 = HJ(FD, [3.0, 3.0]); linha("D barreira extrema", rD3)
rE3 = HJ(FE, [3.0, 3.0], menor=menorE); xE3 = linha("E regra (v, f)", rE3)
c = [np.isinf(rD3.fx) and np.array_equal(rD3.x, [3, 3]), confere(xE3, [1, 1], 4)]
print("  D não sai de x0 (tudo +inf) | E chega à fronteira e para em (1; 1): %s" % (simnao(c),))
ok = ok and all(c)

# ------------------------------------------------------------ 5. o bloqueio em (1, 1)
print("\n5. Porque para em (1, 1)? As 4 tentativas +-P_j (P = 0.5) a partir de (1, 1):")
xb = np.array([1.0, 1.0]); P = 0.5
for d, nome in (([1, 0], "+e1"), ([-1, 0], "-e1"), ([0, 1], "+e2"), ([0, -1], "-e2")):
    xp = xb + P * np.array(d)
    print("  %s: x = (%.1f; %.1f)  f = %.2f  g = %+.1f  C(R=2) = %.2f  D = %s"
          % (nome, xp[0], xp[1], f(xp), g(xp), FC(2)(xp), "%.2f" % FD(xp) if g(xp) >= 0 else "+inf"))
print("  em (1, 1): f = %.2f; nenhuma direção coordenada melhora (qualquer P > 0);"
      " a direção (1, -1)/sqrt 2, ao longo da fronteira, melhora." % f(xb))
xq = xb + 0.25 * np.array([1, -1])                # (1.25, 0.75): na fronteira
c = [all(FC(2)(xb + p * np.array(d)) > f(xb) for p in (0.5, 1e-3, 1e-6)
         for d in ([1, 0], [-1, 0], [0, 1], [0, -1])),
     FD(xq) < f(xb)]
print("  nenhuma tentativa +-P e_j melhora em C (R = 2), P = 0.5, 1e-3, 1e-6 | (1,-1) melhora: %s"
      % (simnao(c),))
ok = ok and all(c)

# ------------------------------------------------------------ 6. a correção: rodar as direções
print("\n6. O mesmo Hooke–Jeeves nas variáveis y, x = Q y, Q = [1 1; 1 -1]/sqrt 2"
      " (direções ao longo e na normal à fronteira)")
Q = np.array([[1.0, 1.0], [1.0, -1.0]]) / np.sqrt(2)
rot = lambda F: (lambda y: F(Q @ y))
y0 = Q.T @ np.array(x0)
rotres = {}
for nome, F, k in (("C exata, R = 2", FC(2), {}), ("D barreira extrema", FD, {}),
                   ("E regra (v, f)", FE, {"menor": menorE})):
    r = hooke_jeeves(rot(F), y0, a, P0, T, **k)
    r.x = Q @ r.x
    rotres[nome] = linha(nome + " (rodado)", r)
c = [confere(list(rotres.values()), [xs] * 3, 4)]
print("  com as direções rodadas, C (R = 2), D e E chegam a x* = (1.5; 0.5): %s" % simnao(c))
ok = ok and all(c)

# ------------------------------------------------------------ 7. informativo: Nelder–Mead
print("\n7. Informativo: o Nelder–Mead do 3.2.3 (simplex inicial de x0 = (0, 0), tolx = 1e-8,"
      " tolf = 1e-10) chega a x* aqui, mas sem garantia geral")
for nome, F in (("C exata, R = 2", FC(2)), ("D barreira extrema", FD)):
    r = nelder_mead(F, np.array(x0), tolx=1e-8, tolf=1e-10, kmax=2000)
    linha(nome + " (Nelder–Mead)", r)

print("\nconfere com os slides: %s" % ("sim" if ok else "não"))
