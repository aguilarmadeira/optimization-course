"""NSGA-II: algoritmo genético multiobjetivo (rank + crowding, elitismo P ∪ Q).

Otimização — deck 6.6 (NSGA-II), frame «O algoritmo».
Reproduz o exemplo dos slides: ver ex06_6_nsga2.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não trata valores de f que sejam NaN; o SBX é a
versão simples, com todas as variáveis cruzadas).  É a implementação do
script das figuras do capítulo 6, com as mesmas chamadas ao gerador
(numpy.random.default_rng), por isso reproduz exatamente os números do deck.
Alternativa com biblioteca: nsga2_exemplo_uc.py (pymoo).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np

from non_dominated_sort import non_dominated_sort
from crowding_distance import crowding_distance
from hypervolume_2d import hypervolume_2d

COLS = ("t", "n_f", "|F_1|", "inadm", "HV")


@dataclass
class NSGA2Result:
    """Resultado de `nsga2` (os mesmos campos que [X, F, info] no MATLAB)."""
    X: np.ndarray             # soluções da frente F_1 de P_G (uma por linha)
    F: np.ndarray             # objetivos dessas soluções
    nit: int                  # gerações (= G)
    nfev: int                 # avaliações do vetor f: N + G*N
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    P: np.ndarray = None      # população final P_G
    FP: np.ndarray = None     # objetivos de P_G
    cvP: np.ndarray = None    # violação das restrições em P_G
    rank: np.ndarray = None   # rank em P_G (1 = F_1)
    crowd: np.ndarray = None  # crowding em P_G (calculado em cada frente)
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 5)))
    cols: tuple = COLS
    flag: int = 0
    message: str = ""


def nsga2(f, lb, ub, N, G, s=None, cv=None, pc=0.9, eta_c=15, pm=None,
          eta_m=20, ref=None, verbose=False):
    """NSGA-II com torneio binário por ≺_n, SBX, mutação polinomial e P ∪ Q.

    Parâmetros
    ----------
    f : função de x (vetor numpy de n componentes) -> vetor de m objetivos
    lb, ub : limites da caixa [x^(L), x^(U)] (vetores de n componentes)
    N : tamanho da população (par)
    G : número de gerações
    s : semente (inteiro) ou um numpy.random.Generator já criado
    cv : (opcional) função de x -> violação das restrições (>= 0; 0 = admissível),
         tratada pela regra de Deb; sem cv, só a caixa (garantida por projeção)
    pc : probabilidade de cruzamento (0.9)
    eta_c, eta_m : índices de distribuição do SBX (15) e da mutação polinomial (20)
    pm : probabilidade de mutação por variável (por omissão 1/n)
    ref : (opcional) ponto de referência r para o HV no histórico (m = 2)
    verbose : se True, imprime o histórico de 5 em 5 gerações (G = 50)

    Devolve
    -------
    NSGA2Result com X, F (frente F_1 de P_G), nfev = N + G*N (a população
    inicial conta), nit = G, a população final e history com uma linha por
    população P_t, t = 0, ..., G: [t, n_f, |F_1|, inadm, HV]
    (|F_1| = não dominados admissíveis; HV da frente F_1 com r = ref).
    """
    lb = np.asarray(lb, float); ub = np.asarray(ub, float); n = lb.size
    pm = pm or 1.0 / n
    rng = s if isinstance(s, np.random.Generator) else np.random.default_rng(s)
    nfev = 0

    def avaliar(X):
        nonlocal nfev
        FX = np.array([np.ravel(f(x)) for x in X])
        CX = np.zeros(len(X)) if cv is None else np.array([float(cv(x)) for x in X])
        nfev += len(X)
        return FX, CX

    def ordenar(FX, CX):
        return non_dominated_sort(FX, None if cv is None else CX)

    def sbx(a, b):                       # cruzamento binário simulado (Deb)
        u = rng.random(n)
        beta = np.where(u <= 0.5, (2 * u)**(1 / (eta_c + 1)),
                        (1 / (2 - 2 * u))**(1 / (eta_c + 1)))
        c1 = 0.5 * ((1 + beta) * a + (1 - beta) * b)
        c2 = 0.5 * ((1 - beta) * a + (1 + beta) * b)
        return np.clip(c1, lb, ub), np.clip(c2, lb, ub)

    def mutar(x):                        # mutação polinomial (Deb)
        mk = rng.random(n) < pm
        u = rng.random(n)
        delta = np.where(u < 0.5, (2 * u)**(1 / (eta_m + 1)) - 1,
                         1 - (2 - 2 * u)**(1 / (eta_m + 1)))
        return np.clip(x + mk * delta * (ub - lb), lb, ub)

    def registo(t, FX, CX, fronts):
        F1 = [i for i in fronts[0] if CX[i] <= 0]      # não dominados admissíveis
        hv = np.nan if ref is None else hypervolume_2d(FX[F1], ref) if F1 else 0.0
        return [t, nfev, len(F1), int((CX > 0).sum()), hv]

    hist = []
    # --- núcleo (o do slide «O algoritmo», com as contagens) --------------
    P = rng.uniform(lb, ub, (N, n))                    # P_0 na caixa
    FP, CP = avaliar(P)                                # avaliar f em P_0
    for t in range(G):
        fronts, rank = ordenar(FP, CP)                 # frentes de P_t
        d = np.zeros(N)
        for fr in fronts:
            d[fr] = crowding_distance(FP, fr)          # d_p em cada frente
        hist.append(registo(t, FP, CP, fronts))

        def vence(i, j):                               # torneio binário por ≺_n
            if rank[i] != rank[j]:
                return i if rank[i] < rank[j] else j
            return i if d[i] >= d[j] else j

        Q = []
        while len(Q) < N:                              # N filhos
            i1, i2 = rng.integers(0, N, 2); j1, j2 = rng.integers(0, N, 2)
            a = P[vence(i1, i2)]; b = P[vence(j1, j2)]
            if rng.random() < pc:
                c1, c2 = sbx(a, b)
            else:
                c1, c2 = a.copy(), b.copy()
            Q += [mutar(c1), mutar(c2)]
        Q = np.array(Q[:N])
        FQ, CQ = avaliar(Q)                            # avaliar Q_t

        R = np.vstack([P, Q]); FR = np.vstack([FP, FQ]); CR = np.concatenate([CP, CQ])
        fronts, _ = ordenar(FR, CR)                    # R_t = P_t ∪ Q_t: novo rank
        novos = []
        for fr in fronts:                              # encher P_{t+1} frente a frente
            if len(novos) + len(fr) <= N:
                novos += fr
            else:                                      # truncar pelos maiores d_p
                dfr = crowding_distance(FR, fr)
                ordem = sorted(range(len(fr)), key=lambda k: -dfr[k])
                novos += [fr[k] for k in ordem[:N - len(novos)]]
                break
        P = R[novos]; FP = FR[novos]; CP = CR[novos]
    # ----------------------------------------------------------------------

    fronts, rank = ordenar(FP, CP)
    d = np.zeros(N)
    for fr in fronts:
        d[fr] = crowding_distance(FP, fr)
    hist.append(registo(G, FP, CP, fronts))
    hist = np.array(hist, float)
    F1 = fronts[0]

    if verbose:
        passo = max(1, G // 10)
        print("%5s %7s %7s %7s %9s" % COLS)
        for row in hist:
            if row[0] % passo == 0 or row[0] == G:
                print("%5d %7d %7d %7d %9.4f" % tuple(row))

    return NSGA2Result(X=P[F1], F=FP[F1], nit=G, nfev=nfev, P=P, FP=FP, cvP=CP,
                       rank=rank, crowd=d, history=hist, flag=0,
                       message="fez as %d gerações pedidas: %d não dominados em P_G"
                               % (G, len(F1)))
