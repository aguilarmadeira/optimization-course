"""Pesquisa aleatória localizada: m amostras em x + r[-1,1]^n; r encolhe sem progresso.

Otimização — deck 3.2.1 (Método de pesquisa aleatória, variante localizada).
Reproduz os exemplos dos slides: ver ex03_2_1_random_search.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não há limites para as variáveis).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class LocalRandomSearchResult:
    """Resultado de `local_random_search` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # melhor ponto (o ponto atual)
    fx: float                 # f(x)
    nit: int                  # iterações feitas (k)
    nfev: int                 # avaliações de f: 1 + m*nit (inclui f(x0))
    nhit: float = np.inf      # 1.ª avaliação com f < ftarget (inf se nunca)
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0
    message: str = ""


def local_random_search(f, x0, r0, m=20, gamma=0.9, tolx=1e-8, kmax=5000,
                        ftarget=-np.inf, s=None, verbose=False):
    """Pesquisa aleatória localizada a partir de x0, com raio inicial r0.

    Parâmetros
    ----------
    f : função de x (vetor numpy de n componentes)
    x0 : ponto inicial
    r0 : raio inicial (as amostras são x + r*xi, xi uniforme em [-1,1]^n)
    m : amostras por iteração (por omissão 20)
    gamma : fator de redução do raio quando não há progresso (0.9)
    tolx : para quando r <= tolx (1e-8)
    kmax : número máximo de iterações (5000)
    ftarget : alvo para nhit (por omissão -inf: sem alvo)
    s : semente (inteiro) ou um numpy.random.Generator já criado; em cada
        iteração, X = x + r*(2*rng.random((m, n)) - 1)
    verbose : se True, imprime a tabela das iterações

    Devolve
    -------
    LocalRandomSearchResult com x, fx, nit, nfev = 1 + m*nit (inclui f(x0)),
    nhit (1.ª avaliação com f < ftarget, contando f(x0) como a 1.ª; é o
    custo até atingir o alvo usado na comparação final do capítulo) e
    history com uma linha por k = 0, ..., nit: [k, n_f, r, f(x), x].
    """
    x = np.array(x0, float); n = x.size
    rng = s if isinstance(s, np.random.Generator) else np.random.default_rng(s)

    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    fx = f(x); r = r0; k = 0
    nfev = 1
    nhit = 1 if fx < ftarget else np.inf
    hist = [[k, nfev, r, fx, *x]]
    while True:
        X = x + r * (2 * rng.random((m, n)) - 1)   # m pontos em x + [-r,r]^n
        v = np.zeros(m)
        for j in range(m):
            v[j] = f(X[j])
            if v[j] < ftarget and nhit == np.inf:
                nhit = nfev + j + 1
        nfev += m
        i = int(np.argmin(v))
        if v[i] < fx:                     # o melhor deles melhora f(x): aceitar
            x = X[i]; fx = v[i]
        else:                             # sem progresso: encolher
            r = gamma * r
        k += 1
        hist.append([k, nfev, r, fx, *x])
        if r <= tolx or k >= kmax:
            break
    # ----------------------------------------------------------------------

    flag = 0 if r <= tolx else 1
    msg = ("r = %.3g <= tolx" % r) if flag == 0 else ("atingiu kmax = %d iterações" % kmax)
    cols = ("k", "n_f", "r", "f(x)") + tuple("x%d" % (j + 1) for j in range(n))
    hist = np.array(hist)
    if verbose:
        print(("%5s %7s %10s %12s" + " %10s" * n) % cols)
        for row in hist:
            print(("%5d %7d %10.3e %12.4e" + " %10.6f" * n) % tuple(row))

    return LocalRandomSearchResult(x=x, fx=fx, nit=k, nfev=nfev, nhit=nhit,
                                   history=hist, cols=cols, flag=flag, message=msg)
