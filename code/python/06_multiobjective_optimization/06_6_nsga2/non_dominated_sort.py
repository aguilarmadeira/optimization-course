"""Ordenação não dominada (fast non-dominated sort) e regra de Deb.

Otimização — deck 6.6 (NSGA-II), frame «Ordenação não dominada: o rank»;
dominância: deck 6.2.  Usada por `nsga2` (ver ex06_6_nsga2.py).
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não trata valores de f que sejam NaN).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import numpy as np


def dominates(a, b):
    """a domina b (minimização): a <= b em todos os objetivos e < em pelo menos um."""
    return bool(np.all(a <= b) and np.any(a < b))


def _fast_sort(F):
    """Ordenação de Deb et al. (2002): S_p = quem p domina, n_p = quantos dominam p."""
    N = len(F)
    S = []                                  # S[p]: índices dominados por p
    n = np.zeros(N, int)                    # n[p]: número de pontos que dominam p
    fronts = [[]]
    for p in range(N):
        domina = np.all(F[p] <= F, axis=1) & np.any(F[p] < F, axis=1)       # p domina q
        dominado = np.all(F <= F[p], axis=1) & np.any(F < F[p], axis=1)     # q domina p
        S.append(list(np.flatnonzero(domina)))
        n[p] = int(dominado.sum())
        if n[p] == 0:
            fronts[0].append(p)             # ninguém domina p: frente F_1
    i = 0
    while fronts[i]:                        # retira F_i e repete nas restantes
        nxt = []
        for p in fronts[i]:
            for q in S[p]:
                n[q] -= 1
                if n[q] == 0:
                    nxt.append(int(q))
        i += 1
        fronts.append(nxt)
    return fronts[:-1]


def non_dominated_sort(F, cv=None):
    """Ordena os pontos em frentes F_1, F_2, ... (rank 1 é o melhor).

    Parâmetros
    ----------
    F : matriz N x m dos valores dos objetivos (uma linha por ponto)
    cv : (opcional) vetor com a violação das restrições de cada ponto
         (cv = 0: admissível).  Com cv, aplica a regra de Deb:
         admissível vence inadmissível; entre inadmissíveis, menor violação;
         entre admissíveis, dominância normal.  Os admissíveis são ordenados
         em frentes; cada inadmissível forma uma frente sozinho, depois, por
         violação crescente (empates pela ordem dos índices).

    Devolve
    -------
    fronts : lista de frentes, cada uma uma lista de índices (base 0)
    rank : vetor numpy com o rank de cada ponto (1 = frente F_1)
    """
    F = np.atleast_2d(np.asarray(F, float))
    N = len(F)
    if cv is None:
        fronts = _fast_sort(F)
    else:
        cv = np.ravel(np.asarray(cv, float))
        adm = np.flatnonzero(cv <= 0)
        inadm = np.flatnonzero(cv > 0)
        fronts = []
        if adm.size:
            fronts = [[int(adm[i]) for i in fr] for fr in _fast_sort(F[adm])]
        for i in inadm[np.argsort(cv[inadm], kind="stable")]:
            fronts.append([int(i)])
    rank = np.zeros(N, int)
    for r, fr in enumerate(fronts, start=1):
        rank[fr] = r
    return fronts, rank
