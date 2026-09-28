"""Movimento exploratório de Hooke–Jeeves em torno de xb.

Percorre as coordenadas j = 1, ..., n, sequencialmente: tenta x_j + P_j;
se não melhora, tenta x_j - P_j; se nenhuma melhora, deixa x_j como estava.
Cada teste parte do melhor ponto encontrado até aí.

Otimização — deck 3.2.5 (Método de Hooke–Jeeves).
Usada por hooke_jeeves.py; exemplo: ex03_2_5_hooke_jeeves.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class ExploratoryResult:
    """Resultado de `exploratory` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # ponto explorado xe (= xb se a exploração falhou)
    fx: float                 # f(xe)
    nit: int                  # coordenadas percorridas (= n)
    nfev: int                 # entre n (todas as tentativas +P_j funcionam) e 2n
    ngev: int = 0
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    ftrace: np.ndarray = field(default_factory=lambda: np.zeros(0))
    flag: int = 0             # 0 se melhorou (fe < fb); 1 se falhou
    message: str = ""


def exploratory(f, xb, fb, P):
    """Exploração eixo a eixo em torno de xb, com perturbações P.

    Parâmetros
    ----------
    f : função de um vetor x (array de dimensão n)
    xb : ponto base
    fb : f(xb), já conhecido (NÃO é avaliado aqui nem contado)
    P : perturbações, uma por coordenada (escalar = igual em todas)

    Devolve
    -------
    ExploratoryResult com x = xe, fx = f(xe), nfev (n a 2n) e history com uma
    linha por avaliação: [j, sinal, xp(1:n), f(xp), aceite], sinal = +1 ou -1
    (x_j + P_j ou x_j - P_j), aceite = 1 ou 0; ftrace com os valores de f pela
    ordem em que foram avaliados; flag = 0 se melhorou, 1 se falhou.
    """
    xb = np.array(xb, dtype=float).ravel()
    n = xb.size
    P = np.array(P, dtype=float).ravel() * np.ones(n)

    xe = xb.copy(); fe = fb               # xe = cópia de xb
    nfev = 0
    hist = []
    ftrace = []

    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    for j in range(n):                    # exploração (xe = cópia de xb)
        xp = xe.copy(); xp[j] = xp[j] + P[j]; fp = f(xp)   # xp: teste
        nfev += 1; hist.append([j + 1, +1, *xp, fp, fp < fe]); ftrace.append(fp)
        if fp >= fe:
            xp[j] = xe[j] - P[j]; fp = f(xp)
            nfev += 1; hist.append([j + 1, -1, *xp, fp, fp < fe]); ftrace.append(fp)
        if fp < fe:
            xe = xp; fe = fp
    # ----------------------------------------------------------------------

    cols = ("j", "sinal") + tuple("xp%d" % (j + 1) for j in range(n)) + ("f(xp)", "aceite")
    if fe < fb:
        flag = 0
        message = "a exploração melhorou: f = %.6g -> %.6g (%d avaliações)" % (fb, fe, nfev)
    else:
        flag = 1
        message = "a exploração falhou: nenhuma das tentativas +-P_j melhorou (%d avaliações)" % nfev

    return ExploratoryResult(x=xe, fx=fe, nit=n, nfev=nfev,
                             history=np.array(hist, dtype=float).reshape(-1, n + 4),
                             cols=cols, ftrace=np.array(ftrace), flag=flag, message=message)
