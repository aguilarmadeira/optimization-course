"""Simulated annealing: min f(x) em X = [lb, ub], com arrefecimento geométrico.

Um único ponto x; em cada iteração propõe-se um vizinho x' = x + sigma*zeta,
zeta ~ N(0, I), projetado em X, e aceita-se pelo critério de Metropolis:
sempre se Delta f <= 0; se Delta f > 0, com probabilidade exp(-Delta f / T).
Depois T <- c*T. Devolve-se o melhor ponto visitado.

Otimização — deck 5.2 (Simulated annealing).
Reproduz os exemplos dos slides: ver ex05_2_simulated_annealing.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não adapta sigma nem T ao longo da corrida,
não tem reaquecimentos nem critério de paragem além de Nmax).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class SAResult:
    """Resultado de `simulated_annealing` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # melhor ponto visitado (x_best)
    fx: float                 # f(x_best)
    nit: int                  # iterações feitas (= Nmax)
    nfev: int                 # avaliações de f: Nmax + 1 (inclui f(x0))
    nacc: int                 # propostas aceites
    T: float                  # temperatura no fim: T0 * c**Nmax
    xk: np.ndarray = field(default_factory=lambda: np.zeros(0))   # ponto corrente final
    fk: float = np.nan        # f no ponto corrente final
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0
    message: str = ""


def simulated_annealing(f, x0, T0, c, sigma, Nmax, lb=-np.inf, ub=np.inf, rng=None):
    """Simulated annealing com arrefecimento geométrico (assinatura do slide).

    Parâmetros
    ----------
    f : função objetivo de x (vetor numpy)
    x0 : ponto inicial (escalar ou vetor)
    T0 : temperatura inicial, T0 > 0 (na escala de f)
    c : fator de arrefecimento, 0 < c < 1 (T <- c*T em cada iteração)
    sigma : desvio-padrão da perturbação gaussiana (passo)
    Nmax : número de iterações (uma avaliação de f por iteração)
    lb, ub : limites de X (escalares ou vetores); o vizinho é projetado em X
    rng : gerador numpy (np.random.Generator) ou semente inteira; None: novo
        gerador. Para repetir uma corrida, passar np.random.default_rng(s).

    Devolve
    -------
    SAResult com x (melhor ponto visitado), fx = f(x), nfev = Nmax + 1 (conta
    f(x0)), nacc (propostas aceites), T (final), xk e fk (ponto corrente final)
    e history com uma linha por k = 0, 1, ..., Nmax:
    [k, T, aceite, f(x'), f(x_k), f_best, x'_1..x'_n, x_k,1..x_k,n]
    (T é a temperatura usada na iteração k, T0*c**(k-1); x' é a proposta e
    x_k o ponto corrente depois da iteração k; na linha 0, x' = x_0 e T = T0).
    """
    if rng is None or np.isscalar(rng):            # semente -> gerador
        rng = np.random.default_rng(rng)
    x = np.atleast_1d(np.asarray(x0, float)).copy()
    n = x.size
    hist = np.zeros((Nmax + 1, 6 + 2 * n))

    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    fx = f(x); nfev = 1
    xbest = x.copy(); fbest = fx; T = T0; nacc = 0
    hist[0] = np.concatenate(([0, T0, 1, fx, fx, fbest], x, x))
    for k in range(1, Nmax + 1):
        xn = x + sigma * rng.normal(size=x.shape)          # vizinho
        xn = np.minimum(np.maximum(xn, lb), ub)            # projetar em X
        fn = f(xn); df = fn - fx
        nfev += 1                                          # uma avaliação por iteração
        if df <= 0 or rng.random() < np.exp(-df / T):      # Metropolis
            x = xn; fx = fn; aceite = 1; nacc += 1
            if fx < fbest:
                xbest = x.copy(); fbest = fx
        else:
            aceite = 0
        hist[k] = np.concatenate(([k, T, aceite, fn, fx, fbest], xn, x))
        T = c * T                                          # arrefecer
    # ----------------------------------------------------------------------

    cols = ("k", "T", "aceite", "f(x')", "f(x_k)", "f_best") + \
        tuple("x'_%d" % (i + 1) for i in range(n)) + tuple("x_k,%d" % (i + 1) for i in range(n))
    return SAResult(x=xbest, fx=fbest, nit=Nmax, nfev=nfev, nacc=nacc, T=T, xk=x, fk=fx,
                    history=hist, cols=cols, flag=0,
                    message="fez as %d iterações pedidas: %d aceites, T final = %.3g"
                            % (Nmax, nacc, T))
