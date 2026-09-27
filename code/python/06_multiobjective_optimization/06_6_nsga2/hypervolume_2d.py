"""Hipervolume em 2D (minimização): área dominada pelos pontos e limitada por r.

Otimização — deck 6.2 (frame «Medir a qualidade: o hipervolume») e deck 6.6.
Usada por `nsga2` e por ex06_6_nsga2.py.
Implementação didática — só para m = 2 objetivos.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import numpy as np


def hypervolume_2d(F, r):
    """Hipervolume 2D de F (N x 2) em relação ao ponto de referência r.

    Só contam os pontos que dominam estritamente r (f_1 < r_1 e f_2 < r_2);
    retiram-se os dominados, ordena-se por f_1 e somam-se os retângulos
        HV = soma_k (r_1 - f_1^(k)) * (f_2^(k-1) - f_2^(k)),  f_2^(0) = r_2.
    Devolve 0 se nenhum ponto dominar r.
    """
    F = np.atleast_2d(np.asarray(F, float))
    r = np.asarray(r, float)
    F = F[np.all(F < r, axis=1)]               # só os que dominam r
    if F.size == 0:
        return 0.0
    nd = np.array([not np.any(np.all(F <= z, axis=1) & np.any(F < z, axis=1))
                   for z in F])                # não dominados
    F = F[nd]
    F = F[np.argsort(F[:, 0])]
    hv = 0.0
    prev = r[1]
    for f1, f2 in F:
        hv += (r[0] - f1) * (prev - f2)
        prev = f2
    return float(hv)
