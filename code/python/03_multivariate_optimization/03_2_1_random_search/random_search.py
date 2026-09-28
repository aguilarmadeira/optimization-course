"""Pesquisa aleatória pura: min f(x) em D = [a, b], com N amostras uniformes.

Otimização — deck 3.2.1 (Método de pesquisa aleatória).
Reproduz os exemplos dos slides: ver ex03_2_1_random_search.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não trata valores de f que sejam NaN).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class RandomSearchResult:
    """Resultado de `random_search` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # melhor ponto encontrado (x_best)
    fx: float                 # f(x_best)
    nit: int                  # amostras geradas (= N)
    nfev: int                 # avaliações de f (= N)
    nhit: float = np.inf      # 1.ª avaliação com f < ftarget (inf se nunca)
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0
    message: str = ""


def random_search(f, a, b, N, s=None, ftarget=-np.inf, verbose=False):
    """Pesquisa aleatória pura em D = [a, b] (caixa em R^n) com N amostras.

    Parâmetros
    ----------
    f : função de x (vetor numpy de n componentes)
    a, b : limites inferior e superior de D (vetores de n componentes)
    N : número de amostras (n_f = N)
    s : semente (inteiro) ou um numpy.random.Generator já criado; os pontos
        são gerados com numpy.random.default_rng(s), um de cada vez,
        x = a + (b - a) * rng.random(n)
    ftarget : alvo para info.nhit (por omissão -inf: sem alvo)
    verbose : se True, imprime as melhorias do melhor ponto

    Devolve
    -------
    RandomSearchResult com x = x_best, fx = f_best, nfev = nit = N,
    nhit (1.ª avaliação com f < ftarget; inf se nunca) e history com uma
    linha por melhoria: [i, x_best, f_best] (o melhor até agora nunca piora).
    """
    a = np.asarray(a, float); b = np.asarray(b, float); n = a.size
    rng = s if isinstance(s, np.random.Generator) else np.random.default_rng(s)
    nfev = 0; nhit = np.inf; hist = []

    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    fbest = np.inf; xbest = None
    for i in range(1, N + 1):
        x = a + (b - a) * rng.random(n)   # em [a, b]
        fx = f(x)
        nfev += 1
        if fx < ftarget and nhit == np.inf:
            nhit = nfev
        if fx < fbest:
            xbest = x
            fbest = fx
            hist.append([i, *xbest, fbest])
    # ----------------------------------------------------------------------

    cols = ("i",) + tuple("x%d" % (j + 1) for j in range(n)) + ("f_best",)
    hist = np.array(hist)
    if verbose:
        print(("%6s" + " %10s" * (n + 1)) % cols)
        for row in hist:
            print(("%6d" + " %10.4f" * (n + 1)) % tuple(row))

    return RandomSearchResult(x=xbest, fx=fbest, nit=N, nfev=nfev, nhit=nhit,
                              history=hist, cols=cols, flag=0,
                              message="gerou as %d amostras pedidas" % N)
