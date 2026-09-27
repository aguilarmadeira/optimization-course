"""Pesquisa em grelha: min f(x1, x2) em [a1,b1] x [a2,b2], m valores por variável.

Otimização — deck 3.2.2 (Pesquisa em grelha).
Reproduz os exemplos dos slides: ver ex03_2_2_grid_search.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. só serve para n = 2, como no deck).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class GridSearchResult:
    """Resultado de `grid_search` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # melhor nó da grelha
    fx: float                 # f(x)
    nfev: int                 # avaliações de f: m^2
    nit: int = 0              # não há iterações: uma só passagem pela grelha
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    h: np.ndarray = field(default_factory=lambda: np.zeros(2))   # espaçamento
    X1: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    X2: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    F: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    cols: tuple = ("x1", "x2", "f")
    flag: int = 0
    message: str = ""


def grid_search(f, a, b, m, verbose=False):
    """Pesquisa em grelha m x m em [a1, b1] x [a2, b2] (n = 2, como no deck).

    Parâmetros
    ----------
    f : função de x = np.array([x1, x2])
    a, b : limites inferior e superior, a = (a1, a2), b = (b1, b2)
    m : número de valores por variável (m^2 avaliações)
    verbose : se True, imprime os 5 melhores nós

    Devolve
    -------
    GridSearchResult com x (melhor nó), fx, nfev = m^2, nit = 0, h (espaçamento
    (b - a)/(m - 1)), X1, X2, F (a grelha do meshgrid e os valores de f) e
    history com todos os nós, do melhor para o pior: [x1, x2, f].
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size != 2 or b.size != 2:
        raise ValueError("grid_search: esta versão é para n = 2.")
    if not (np.all(a < b) and m >= 2):
        raise ValueError("grid_search: é preciso a < b e m >= 2.")

    # --- núcleo (o mesmo dos slides) ----------------------------------------
    X1, X2 = np.meshgrid(np.linspace(a[0], b[0], m),
                         np.linspace(a[1], b[1], m))
    F = np.vectorize(lambda u, v: f(np.array([u, v])), otypes=[float])(X1, X2)
    i = np.argmin(F); fmin = F.flat[i]; x = np.array([X1.flat[i], X2.flat[i]])
    # ----------------------------------------------------------------------

    j = np.argsort(F.ravel(), kind="stable")
    hist = np.c_[X1.ravel()[j], X2.ravel()[j], F.ravel()[j]]
    if verbose:
        print("%4s %10s %10s %12s" % ("i", "x1", "x2", "f"))
        for k, row in enumerate(hist[:5], 1):
            print("%4d %10.4f %10.4f %12.4f" % (k, *row))

    return GridSearchResult(x=x, fx=fmin, nfev=F.size, h=(b - a) / (m - 1),
                            X1=X1, X2=X2, F=F, history=hist, flag=0,
                            message="avaliou a grelha %d x %d (%d avaliações)" % (m, m, F.size))
