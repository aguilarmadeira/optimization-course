"""Método de Newton-Raphson em 1D: procura um zero de f'.

Otimização — deck 2.5 (Newton-Raphson).
Reproduz o exemplo dos slides: ver ex02_5_newton.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. passo amortecido ou salvaguarda por bisseção).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

COLS = ("k", "x_k", "f(x_k)", "df(x_k)", "ddf(x_k)", "x_k+1", "df(x_k+1)")


@dataclass
class NewtonResult:
    """Resultado de `newton_1d` (os mesmos campos que `info` no MATLAB)."""
    x: float                  # última iterada
    fx: float                 # f(x) (nan se f is None) -- só para leitura
    nit: int                  # passos de Newton
    ngev: int                 # avaliações de f': 1 + nit
    nhev: int                 # avaliações de f'': nit
    nfev: int = 0             # o método não usa f
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 7)))
    cols: tuple = COLS
    flag: int = 0             # 0: |f'| < tolg; 1: atingiu kmax
    message: str = ""


def newton_1d(f, df, ddf, x0, tolg, kmax, tolH=1e-12, verbose=False):
    """Newton-Raphson para f'(x) = 0 a partir de x0.

    Para quando |f'(x_{k+1})| < tolg (flag 0) ou ao fim de kmax iterações
    (flag 1, com aviso). Erro se |f''(x_k)| < tolH (curvatura ~ 0).
    Newton procura zeros de f', não mínimos: no fim, classificar o ponto com
    f''(x) (> 0 mínimo, < 0 máximo); essa avaliação extra não é feita aqui.
    f só é usada para leitura (coluna f(x_k) e fx) e não conta em nfev;
    com f=None nada disso é calculado.

    history tem uma linha por k = 0, ..., nit-1:
    [k, x_k, f(x_k), f'(x_k), f''(x_k), x_{k+1}, f'(x_{k+1})].
    """
    temf = f is not None
    x = x0
    hist = []
    flag = 1
    ngev = 0; nhev = 0
    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    g = df(x)
    ngev += 1
    for k in range(1, kmax + 1):
        H = ddf(x)
        nhev += 1
        if abs(H) < tolH:
            raise ArithmeticError("curvatura demasiado pequena")
        xk, gk = x, g
        x = x - g / H                 # passo de Newton
        g = df(x)                     # 1 avaliação de f'
        ngev += 1
        fk = f(xk) if temf else np.nan          # só para leitura
        hist.append([k - 1, xk, fk, gk, H, x, g])
        if abs(g) < tolg:
            flag = 0
            break
    # ----------------------------------------------------------------------

    fx = f(x) if temf else np.nan               # só para leitura
    hist = np.array(hist, float).reshape(-1, 7)
    nit = hist.shape[0]
    if flag == 0:
        message = "|f'(x)| < tolg em %d iterações" % nit
    else:
        message = "atingiu o número máximo de iterações (kmax = %d) sem |f'| < tolg" % kmax
        warnings.warn("newton_1d: " + message)

    if verbose:
        print("%3s %13s %13s %14s %12s %13s %14s" % COLS)
        for row in hist:
            print("%3d %13.7f %13.6f %+14.6e %12.5f %13.7f %+14.6e" % tuple(row))

    return NewtonResult(x=x, fx=fx, nit=nit, ngev=ngev, nhev=nhev,
                        history=hist, flag=flag, message=message)
