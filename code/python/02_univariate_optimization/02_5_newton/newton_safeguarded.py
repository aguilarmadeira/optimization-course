"""Newton salvaguardado (bisseção + Newton): procura um zero de f' em [a, b].

Otimização — deck 2.5 (Newton-Raphson), página «Newton salvaguardado».
Reproduz o exemplo dos slides: ver ex02_5_newton_safeguarded.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. tolerância no comprimento do intervalo).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

COLS = ("k", "x_k", "x_Newton", "tipo", "x_k+1", "a", "b", "df(x_k+1)")
TIPOS = {0: "Newton", 1: "bisseção"}


@dataclass
class SafeguardedResult:
    """Resultado de `newton_safeguarded` (os mesmos campos que `info` no MATLAB)."""
    x: float                  # último iterando
    fx: float                 # f(x) (nan se f is None) -- só para leitura
    nit: int                  # iterações
    ngev: int                 # avaliações de f': 2 (f'(a), f'(b)) + 1 (f'(x0)) + nit
    nhev: int                 # avaliações de f'': nit
    nfev: int = 0             # o método não usa f
    nbis: int = 0             # passos de bisseção
    a: float = np.nan         # intervalo final [a, b]
    b: float = np.nan
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 8)))
    cols: tuple = COLS
    flag: int = 0             # 0: |f'| < tolg; 1: atingiu kmax
    message: str = ""


def newton_safeguarded(f, df, ddf, a, b, x0, tolg, kmax, verbose=False):
    """Newton salvaguardado para f'(x) = 0 em [a, b], com f'(a) < 0 < f'(b).

    Em cada iteração: (1) tenta o passo de Newton x_N = x_k - f'(x_k)/f''(x_k);
    (2) se x_N sair de (a, b), ou se não der progresso suficiente
    (|x_N - x_k| > metade do passo de duas iterações antes, o teste de
    rtsafe), faz um passo de bisseção, x_{k+1} = (a + b)/2; (3) atualiza o
    intervalo com o sinal de f'(x_{k+1}): negativo -> a = x_{k+1}, senão
    b = x_{k+1}.  Pára quando |f'(x_{k+1})| < tolg (flag 0) ou ao fim de
    kmax iterações (flag 1, com aviso).

    f só é usada para leitura (fx) e não conta; pode ser None.
    Contagens: n_g = 3 + nit (f'(a), f'(b) para verificar o enquadramento,
    f'(x0) e uma por iteração), n_H = nit.
    history tem uma linha por k = 0, ..., nit-1:
    [k, x_k, x_N, tipo (0 Newton, 1 bisseção), x_{k+1}, a, b, f'(x_{k+1})].
    """
    ga = df(a); gb = df(b)
    ngev = 2
    if not (ga < 0 < gb):
        raise ValueError("newton_safeguarded: é preciso f'(a) < 0 < f'(b)")
    if not (a <= x0 <= b):
        raise ValueError("newton_safeguarded: x0 tem de estar em [a, b]")
    x = x0
    hist = []
    flag = 1
    nhev = 0; nbis = 0
    dxold = b - a; dx = dxold
    # --- núcleo (o dos slides, com as contagens) --------------------------
    g = df(x)
    ngev += 1
    for k in range(1, kmax + 1):
        H = ddf(x)
        nhev += 1
        xN = x - g / H if H != 0 else np.inf     # (1) passo de Newton
        if not (a < xN < b) or abs(xN - x) > abs(dxold) / 2:
            xn = (a + b) / 2                     # (2) bisseção
            tipo = 1; nbis += 1
        else:
            xn = xN
            tipo = 0
        dxold = dx; dx = xn - x
        gn = df(xn)                              # 1 avaliação de f'
        ngev += 1
        if gn < 0:                               # (3) novo intervalo
            a = xn
        else:
            b = xn
        hist.append([k - 1, x, xN, tipo, xn, a, b, gn])
        x = xn; g = gn
        if abs(g) < tolg:
            flag = 0
            break
    # ----------------------------------------------------------------------

    fx = f(x) if f is not None else np.nan       # só para leitura
    hist = np.array(hist, float).reshape(-1, 8)
    nit = hist.shape[0]
    if flag == 0:
        message = "|f'(x)| < tolg em %d iterações (%d de bisseção)" % (nit, nbis)
    else:
        message = "atingiu o número máximo de iterações (kmax = %d) sem |f'| < tolg" % kmax
        warnings.warn("newton_safeguarded: " + message)

    if verbose:
        print("%3s %12s %12s %9s %12s %24s %11s" % ("k", "x_k", "x_Newton", "tipo", "x_k+1", "[a, b]", "df(x_k+1)"))
        for r in hist:
            print("%3d %12.4e %12.4e %9s %12.4e  [%10.3e; %10.3e] %+11.2e"
                  % (r[0], r[1], r[2], TIPOS[int(r[3])], r[4], r[5], r[6], r[7]))

    return SafeguardedResult(x=x, fx=fx, nit=nit, ngev=ngev, nhev=nhev, nbis=nbis,
                             a=a, b=b, history=hist, flag=flag, message=message)
