"""Método da bisseção sobre f': min f(x) em [a, b], com N reduções.

Otimização — deck 2.4 (Bisseção).
Reproduz o exemplo dos slides: ver ex02_4_bisection.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não verifica a unicidade do zero de f').

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np

COLS = ("k", "a_k", "b_k", "x_k", "f(x_k)", "df(x_k)")


@dataclass
class BisectionResult:
    """Resultado de `bisection` (os mesmos campos que `info` no MATLAB)."""
    x: float                  # ponto médio do intervalo final
    fx: float                 # f(x) (nan se f is None) -- só para leitura
    nit: int                  # reduções feitas (N, ou menos se f'(x_k) = 0)
    ngev: int                 # avaliações de f': N + 2
    L: float                  # comprimento do intervalo final
    nfev: int = 0             # o método só usa o sinal de f'
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 6)))
    cols: tuple = COLS
    flag: int = 0             # 0: fez as N reduções; 2: f'(x_k) = 0
    message: str = ""


def bisection(f, df, a, b, N, verbose=False):
    """Bisseção sobre f' em [a, b] com N reduções; exige df(a) < 0 < df(b).

    f só é usada para leitura (coluna f(x_k) da tabela e fx); o método não a
    usa e essas avaliações não entram em nfev (n_f = 0, como nos slides).
    Com f=None nada disso é calculado.

    Para garantir b_N - a_N <= tolx: N = ceil(log2((b - a)/tolx)).

    Devolve um BisectionResult com x = ponto médio final, ngev = N + 2 e
    history com uma linha por k = 0, ..., N-1: [k, a_k, b_k, x_k, f(x_k), f'(x_k)].
    """
    if not a < b:
        raise ValueError("bisection: é preciso a < b.")
    temf = f is not None

    # --- verificar o enquadramento: 2 avaliações de f' --------------------
    da = df(a); db = df(b)
    ngev = 2
    if not (da < 0 and db > 0):
        raise ValueError("enquadramento")

    hist = []
    flag = 0
    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    for k in range(1, N + 1):
        x = (a + b) / 2; d = df(x)
        ngev += 1
        fk = f(x) if temf else np.nan        # só para leitura
        hist.append([k - 1, a, b, x, fk, d])
        if d == 0:
            a = x; b = x; flag = 2           # terminou
            break
        elif d < 0:
            a = x                            # mínimo à direita
        else:
            b = x                            # mínimo à esquerda
    x = (a + b) / 2
    # ----------------------------------------------------------------------

    fx = f(x) if temf else np.nan            # só para leitura
    hist = np.array(hist, float).reshape(-1, 6)
    nit = hist.shape[0]
    L = b - a
    if flag == 2:
        message = "f'(x_k) = 0 na iteração k = %d: terminou" % (nit - 1)
    else:
        message = "fez as %d reduções pedidas: b - a = %.2g" % (N, L)

    if verbose:
        print("%3s %9s %9s %9s %10s %10s" % COLS)
        for row in hist:
            print("%3d %9.4f %9.4f %9.4f %10.4f %+10.4f" % tuple(row))

    return BisectionResult(x=x, fx=fx, nit=nit, ngev=ngev, L=L,
                           history=hist, flag=flag, message=message)
