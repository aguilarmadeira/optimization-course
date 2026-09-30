"""Movimento exploratório de Hooke–Jeeves em torno de xb.

Percorre as coordenadas j = 1, ..., n, sequencialmente: avalia x_j + P_j e
x_j - P_j e fica com o melhor dos três pontos (o atual e as duas tentativas),
como em Deb (2012, sec. 3.3.3); se nenhuma tentativa melhora, deixa x_j como
estava. Cada coordenada parte do melhor ponto encontrado até aí.

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
    nfev: int                 # 2n (duas tentativas por coordenada)
    ngev: int = 0
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    ftrace: np.ndarray = field(default_factory=lambda: np.zeros(0))
    flag: int = 0             # 0 se melhorou (fe < fb); 1 se falhou
    message: str = ""


def exploratory(f, xb, fb, P, menor=None):
    """Exploração eixo a eixo em torno de xb, com perturbações P.

    Parâmetros
    ----------
    f : função de um vetor x (array de dimensão n)
    xb : ponto base
    fb : f(xb), já conhecido (NÃO é avaliado aqui nem contado)
    P : perturbações, uma por coordenada (escalar = igual em todas)
    menor : (opcional) função menor(u, v) que diz se o valor u é melhor do que
        v; por omissão, u < v. O método só compara valores, por isso pode usar
        outra ordem: p. ex. a regra de admissibilidade (deck 4.3.3), com f a
        devolver o par (v(x), f(x)); history e ftrace registam o último
        componente.

    Devolve
    -------
    ExploratoryResult com x = xe, fx = f(xe), nfev = 2n e history com uma
    linha por avaliação: [j, sinal, xp(1:n), f(xp), aceite], sinal = +1 ou -1
    (x_j + P_j ou x_j - P_j), aceite = 1 se foi a tentativa escolhida; ftrace
    com os valores de f pela ordem em que foram avaliados; flag = 0 se
    melhorou, 1 se falhou. Empate entre +P_j e -P_j: fica +P_j.
    """
    if menor is None:
        menor = lambda u, v: u < v        # a ordem habitual
    reg = lambda y: y[-1] if np.ndim(y) else y   # valor registado
    xb = np.array(xb, dtype=float).ravel()
    n = xb.size
    P = np.array(P, dtype=float).ravel() * np.ones(n)

    xe = xb.copy(); fe = fb               # xe = cópia de xb
    nfev = 0
    hist = []
    ftrace = []

    # --- núcleo (o dos slides, com as contagens; menor(u, v) = u < v) -----
    for j in range(n):                    # exploração (xe = cópia de xb)
        xp = xe.copy(); xp[j] = xp[j] + P[j]; fp = f(xp)   # +P_j
        xm = xe.copy(); xm[j] = xm[j] - P[j]; fm = f(xm)   # -P_j
        nfev += 2; ftrace += [reg(fp), reg(fm)]
        mais = menor(fp, fe) and not menor(fm, fp)        # fp < fe e fp <= fm
        menos = not mais and menor(fm, fe)                # fm < fe e fm < fp
        hist.append([j + 1, +1, *xp, reg(fp), mais])
        hist.append([j + 1, -1, *xm, reg(fm), menos])
        if mais:                          # o melhor dos três: +P_j
            xe = xp; fe = fp
        elif menos:                       # o melhor dos três: -P_j
            xe = xm; fe = fm
    # ----------------------------------------------------------------------

    cols = ("j", "sinal") + tuple("xp%d" % (j + 1) for j in range(n)) + ("f(xp)", "aceite")
    if menor(fe, fb):
        flag = 0
        message = "a exploração melhorou: f = %.6g -> %.6g (%d avaliações)" % (reg(fb), reg(fe), nfev)
    else:
        flag = 1
        message = "a exploração falhou: nenhuma das tentativas +-P_j melhorou (%d avaliações)" % nfev

    return ExploratoryResult(x=xe, fx=fe, nit=n, nfev=nfev,
                             history=np.array(hist, dtype=float).reshape(-1, n + 4),
                             cols=cols, ftrace=np.array(ftrace), flag=flag, message=message)
