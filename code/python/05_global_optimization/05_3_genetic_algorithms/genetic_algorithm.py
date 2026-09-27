"""Algoritmo genético com codificação real: min f(x) em X = [lb, ub].

Seleção por torneio de tamanho s, cruzamento BLX-0.5 (com probabilidade pc,
aos pares), mutação gaussiana por gene (probabilidade pm, desvio sigma),
projeção em X e elitismo (as e melhores soluções passam à geração seguinte).
Cada geração avalia os N descendentes: n_f = N (G + 1) no total.

Otimização — deck 5.3 (Algoritmos genéticos).
Reproduz os exemplos dos slides: ver ex05_3_genetic_algorithms.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. reavalia as elites em vez de guardar os valores
já calculados, não tem critério de paragem por estagnação nem por
diversidade, e não trata restrições além dos limites).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class GAResult:
    """Resultado de `genetic_algorithm` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # melhor solução da população final
    fx: float                 # f(x)
    nit: int                  # gerações feitas (= G)
    nfev: int                 # avaliações de f: N (G + 1)
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    P: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))   # população final
    F: np.ndarray = field(default_factory=lambda: np.zeros(0))        # f na população final
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0
    message: str = ""


def cruzar_blx(pais, pc, rng, beta=0.5):
    """Cruzamento BLX-beta aos pares (1-2, 3-4, ...), cada par com probabilidade pc.

    x = a x_a + (1 - a) x_b, com a ~ U(-beta, 1 + beta) por coordenada (pode
    sair do segmento entre os pais); o segundo descendente troca os papéis.
    Sem cruzamento, os descendentes são cópias dos pais. Com N ímpar, o último
    pai passa sem par.
    """
    N, n = pais.shape
    Q = pais.copy()
    for j in range(0, N - 1, 2):
        if rng.random() < pc:
            a = rng.uniform(-beta, 1 + beta, n)
            Q[j] = a * pais[j] + (1 - a) * pais[j + 1]
            Q[j + 1] = (1 - a) * pais[j] + a * pais[j + 1]
    return Q


def genetic_algorithm(f, n, lb, ub, N, G, pc, pm, sigma, e, s=2, rng=None):
    """AG real (assinatura do slide, com o tamanho do torneio s opcional).

    Parâmetros
    ----------
    f : função objetivo de x (vetor numpy com n componentes)
    n : número de variáveis
    lb, ub : limites de X (escalares ou vetores com n componentes)
    N : tamanho da população;  G : número de gerações
    pc : probabilidade de cruzamento (por par);  pm : probabilidade de mutação (por gene)
    sigma : desvio-padrão da mutação gaussiana
    e : número de elites (0 <= e <= N)
    s : tamanho do torneio (por omissão 2, como no slide)
    rng : gerador numpy (np.random.Generator) ou semente inteira; None: novo gerador

    Devolve
    -------
    GAResult com x, fx (o melhor da população final), nfev = N (G + 1),
    nit = G, P e F (população final e os seus valores de f) e history com uma
    linha por geração t = 0, 1, ..., G: [t, melhor f, média de f, x_melhor_1..n].
    """
    if rng is None or np.isscalar(rng):            # semente -> gerador
        rng = np.random.default_rng(rng)
    aval = lambda P: np.array([f(P[i]) for i in range(P.shape[0])])
    hist = np.zeros((G + 1, 3 + n))

    P = rng.uniform(lb, ub, (N, n))                # população inicial em X
    F = aval(P)                                    # F: valores de f (menor = melhor)
    nfev = N
    ib = np.argmin(F); hist[0] = np.concatenate(([0, F[ib], F.mean()], P[ib]))

    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    for t in range(1, G + 1):
        o = np.argsort(F); E = P[o[:e]].copy()     # elites
        I = rng.integers(0, N, (s, N))             # torneio: s candidatos por pai
        i1 = I[0].copy()                           # (s = 2: i1, i2 do slide)
        for i2 in I[1:]:
            m = F[i2] < F[i1]
            i1[m] = i2[m]                          # fica o melhor
        pais = P[i1]
        Q = cruzar_blx(pais, pc, rng)              # BLX-0.5
        M = rng.random((N, n)) < pm                # mutação
        Q = Q + M * sigma * rng.normal(size=(N, n))  # gaussiana
        Q = np.minimum(np.maximum(Q, lb), ub)      # projetar em X
        Q[:e] = E                                  # elites
        P = Q; F = aval(P)
        nfev += N                                  # N avaliações por geração
        ib = np.argmin(F); hist[t] = np.concatenate(([t, F[ib], F.mean()], P[ib]))
    # ----------------------------------------------------------------------

    cols = ("t", "melhor f", "média f") + tuple("x_%d" % (i + 1) for i in range(n))
    return GAResult(x=P[ib].copy(), fx=F[ib], nit=G, nfev=nfev, P=P, F=F,
                    history=hist, cols=cols, flag=0,
                    message="fez as %d gerações pedidas: melhor f = %.4g" % (G, F[ib]))
