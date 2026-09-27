"""Método de Box (versão determinística simplificada de EVOP).

Em cada iteração avalia f nos 2^n vértices x + (1/2)(±Delta_1, ..., ±Delta_n)
da caixa centrada no ponto atual; se o melhor vértice melhora f(x), move-se
para ele (Delta fica igual); senão encolhe a caixa, Delta <- Delta/2.
Para quando ||Delta|| < tolx.

Otimização — deck 3.2.4 (Método de Box — versão simplificada de EVOP).
Reproduz os exemplos dos slides: ver ex03_2_4_box.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não reaproveita vértices já avaliados e o
custo 2^n por iteração torna-o impraticável para n grande).
Uma caixa pequena não certifica estacionariedade a Delta finito.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class BoxResult:
    """Resultado de `box_evo` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # último centro da caixa
    fx: float                 # f(x) (já avaliado; não há avaliação extra no fim)
    nit: int                  # iterações (cada uma: 2^n avaliações e mover ou encolher)
    nfev: int                 # 1 + 2^n * nit (f(x0) conta)
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    ftrace: np.ndarray = field(default_factory=lambda: np.zeros(0))
    Delta: np.ndarray = field(default_factory=lambda: np.zeros(0))
    flag: int = 0
    message: str = ""


def box_evo(f, x0, Delta, tolx, kmax=10000, verbose=False):
    """Método de Box a partir de x0 com lados Delta, até ||Delta|| < tolx.

    Parâmetros
    ----------
    f : função de um vetor x (array de dimensão n)
    x0 : ponto inicial
    Delta : lados da caixa, um por coordenada (escalar = igual em todas)
    tolx : tolerância: para quando ||Delta|| < tolx
    kmax : número máximo de iterações
    verbose : se True, imprime a tabela das iterações

    Devolve
    -------
    BoxResult com x, fx, nit, nfev = 1 + 2^n * nit (f(x0) conta; avaliam-se
    SEMPRE os 2^n vértices, sem reaproveitar pontos — é a convenção das
    contagens do deck 3.2.4), history com uma linha por iteração
    k = 0, ..., nit-1 (o centro x_k e a decisão tomada a partir dele):
    [k, x_k(1:n), Delta(1:n), f(x_k), decisão, f_novo, nfev], com
    decisão = 1 (mover) ou 0 (encolher) e nfev acumulado no fim da iteração;
    ftrace com os valores de f pela ordem em que foram avaliados
    (p. ex. np.argmax(ftrace < 1e-4) + 1 = avaliações até f < 1e-4);
    flag = 0 se ||Delta|| < tolx, 1 se atingiu kmax.
    """
    x = np.array(x0, dtype=float).ravel()
    n = x.size
    Delta = np.array(Delta, dtype=float).ravel() * np.ones(n)
    if np.any(Delta <= 0):
        raise ValueError("box_evo: os lados Delta têm de ser positivos.")

    fx = f(x); nfev = 1                   # f(x0) conta
    ftrace = [fx]
    hist = []

    # linhas de +-1, pela ordem de dec2bin(0:2^n-1) do MATLAB
    # (em 2D: (-1,-1), (-1,1), (1,-1), (1,1))
    S = np.array([[(i >> (n - 1 - j)) & 1 for j in range(n)] for i in range(2**n)], dtype=float)
    S = 2 * S - 1
    k = 0
    while np.linalg.norm(Delta) >= tolx and k < kmax:
        xk, Dk, fk = x.copy(), Delta.copy(), fx
        # --- núcleo (o mesmo dos slides, com as contagens) ----------------
        V = x + 0.5 * S * Delta           # vértices
        Fv = np.array([f(V[i]) for i in range(2**n)])
        i = int(np.argmin(Fv)); fmin = Fv[i]
        if fmin < fx:
            x = V[i].copy(); fx = fmin    # mover
            dec = 1
        else:
            Delta = Delta / 2             # encolher
            dec = 0
        # ------------------------------------------------------------------
        ftrace.extend(Fv)
        nfev += 2**n                      # avaliam-se sempre os 2^n vértices
        hist.append([k, *xk, *Dk, fk, dec, fx, nfev])
        k += 1

    cols = (("k",) + tuple("x%d" % (j + 1) for j in range(n))
            + tuple("Delta%d" % (j + 1) for j in range(n))
            + ("f(x_k)", "decisão", "f_novo", "nfev"))
    hist = np.array(hist, dtype=float).reshape(-1, 2 * n + 5)

    if np.linalg.norm(Delta) < tolx:
        flag = 0
        message = "norm(Delta) = %.3g < tolx = %.3g ao fim de %d iterações" % (
            np.linalg.norm(Delta), tolx, k)
    else:
        flag = 1
        message = "atingiu kmax = %d iterações (norm(Delta) = %.3g)" % (kmax, np.linalg.norm(Delta))

    if verbose:
        nomes = ("encolhe", "move")
        print("%5s" % "k" + "".join("%10s" % c for c in cols[1:2 * n + 1])
              + " %11s   decisão %12s %7s" % ("f(x_k)", "f_novo", "n_f"))
        for r in hist:
            print("%5d" % r[0] + "".join("%10.4f" % v for v in r[1:2 * n + 1])
                  + "%12.4f %9s %12.4f %7d" % (r[2 * n + 1], nomes[int(r[2 * n + 2])],
                                               r[2 * n + 3], r[2 * n + 4]))

    return BoxResult(x=x, fx=fx, nit=k, nfev=nfev, history=hist, cols=cols,
                     ftrace=np.array(ftrace), Delta=Delta, flag=flag, message=message)
