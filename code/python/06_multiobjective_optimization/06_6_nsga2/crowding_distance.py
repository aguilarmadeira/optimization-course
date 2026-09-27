"""Distância de aglomeração (crowding distance) numa frente.

Otimização — deck 6.6 (NSGA-II), frame «Distância de aglomeração: o crowding».
Usada por `nsga2` (ver ex06_6_nsga2.py).
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import numpy as np


def crowding_distance(F, idx):
    """Distância de aglomeração dos pontos idx (uma frente) de F.

    Os extremos recebem d = inf; para os outros, ordenando por cada objetivo i,
        d_p = soma_i (f_i(seguinte) - f_i(anterior)) / (f_i^max - f_i^min).
    Frentes com 1 ou 2 pontos: d = inf para todos.

    Parâmetros
    ----------
    F : matriz N x m dos objetivos
    idx : índices (base 0) dos pontos da frente

    Devolve
    -------
    d : vetor numpy com d_p, pela ordem de idx (maior d = região menos povoada)
    """
    F = np.atleast_2d(np.asarray(F, float))
    idx = list(idx)
    m = F.shape[1]
    d = np.zeros(len(idx))
    if len(idx) <= 2:
        d[:] = np.inf
        return d
    for i in range(m):
        o = np.argsort(F[idx, i])            # ordenar a frente pelo objetivo i
        v = F[idx, i][o]
        amp = v[-1] - v[0]                   # f_i^max - f_i^min na frente
        d[o[0]] = d[o[-1]] = np.inf          # extremos
        if amp > 0:
            for k in range(1, len(idx) - 1):
                d[o[k]] += (v[k + 1] - v[k - 1]) / amp
    return d
