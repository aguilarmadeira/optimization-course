"""Método da secção áurea: min f(x) em [a, b], com N reduções.

Otimização — deck 2.2 (Método da secção áurea).
Reproduz o exemplo dos slides: ver ex02_2_golden_section.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não verifica a unimodalidade de f).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np

COLS = ("k", "a", "b", "x1", "x2", "f(x1)", "f(x2)", "L_k")


@dataclass
class GoldenResult:
    """Resultado de `golden_section` (os mesmos campos que `info` no MATLAB)."""
    x: float                  # ponto médio do intervalo final
    fx: float                 # f(x)
    nit: int                  # reduções feitas (= N)
    nfev: int                 # avaliações de f: N + 3 (inclui f(x))
    nfev_loc: int             # avaliações de f para localizar x*: N + 2
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 8)))
    cols: tuple = COLS
    flag: int = 0
    message: str = ""


def golden_section(f, a, b, N, verbose=False):
    """Secção áurea em [a, b] com N reduções do intervalo.

    Parâmetros
    ----------
    f : função de uma variável, unimodal em [a, b]
    a, b : extremos do intervalo inicial, a < b
    N : número de reduções (L_N = tau**N * (b - a)); para garantir
        L_N <= tolx: N = ceil(log(tolx/(b - a)) / log(tau))
    verbose : se True, imprime a tabela das iterações

    Devolve
    -------
    GoldenResult com x (ponto médio final, |x - x*| <= L_N/2), fx = f(x),
    nfev = N + 3 (N + 2 para localizar x*, mais f(x)), nfev_loc = N + 2,
    nit = N e history com uma linha por k = 0, ..., N:
    [k, a, b, x1, x2, f(x1), f(x2), L_k].
    """
    if not a < b:
        raise ValueError("golden_section: é preciso a < b.")
    nfev = 0
    hist = np.zeros((N + 1, 8))

    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    tau = (np.sqrt(5) - 1) / 2
    L = b - a
    x1 = b - tau * L; f1 = f(x1)          # x1 < x2
    x2 = a + tau * L; f2 = f(x2)
    nfev += 2
    hist[0] = [0, a, b, x1, x2, f1, f2, L]
    for k in range(1, N + 1):
        if f1 < f2:                       # mínimo em [a, x2]
            b = x2; x2 = x1; f2 = f1
            L = b - a; x1 = b - tau * L; f1 = f(x1)
        else:                             # mínimo em [x1, b]
            a = x1; x1 = x2; f1 = f2
            L = b - a; x2 = a + tau * L; f2 = f(x2)
        nfev += 1                         # uma avaliação nova por redução
        hist[k] = [k, a, b, x1, x2, f1, f2, L]
    x = (a + b) / 2
    # ----------------------------------------------------------------------

    nfev_loc = nfev                       # N + 2
    fx = f(x)                             # avaliação no ponto médio final
    nfev += 1                             # N + 3

    if verbose:
        print("%3s %9s %9s %9s %9s %10s %10s %9s" % COLS)
        for row in hist:
            print("%3d %9.4f %9.4f %9.4f %9.4f %10.4f %10.4f %9.4f" % tuple(row))

    return GoldenResult(x=x, fx=fx, nit=N, nfev=nfev, nfev_loc=nfev_loc,
                        history=hist, flag=0,
                        message="fez as %d reduções pedidas: L_N = %.4g" % (N, L))
