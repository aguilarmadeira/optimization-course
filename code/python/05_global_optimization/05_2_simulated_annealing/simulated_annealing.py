"""Simulated annealing como nas aulas (Deb, 2012, cap. 6): min f(x) em X = [lb, ub].

Passo 1: escolher x0, a temperatura inicial T0 (alta), o número n de pontos
         aceites a cada temperatura, o fator alpha (0 < alpha < 1, nas aulas
         entre 0.5 e 0.99) e a temperatura mínima Ts; fazer t = 0.
Passo 2: gerar um vizinho x' = N(x_t) (aqui x' = x_t + sigma*zeta,
         zeta ~ N(0, I), projetado em X).
Passo 3: Delta f = f(x') - f(x_t). Se Delta f < 0, aceitar (x_{t+1} = x',
         t = t + 1); senão, gerar r ~ U(0, 1): se r < exp(-Delta f / T),
         aceitar; senão, voltar ao passo 2.
Passo 4: se T < Ts, terminar; senão, se t mod n = 0, fazer T = alpha*T.
         Voltar ao passo 2.

O contador t só avança quando se aceita um ponto: a temperatura baixa depois
de n pontos ACEITES, e a T baixa (quase tudo rejeitado) cada patamar custa
muitas avaliações. Por isso há também um limite nfmax de avaliações (flag 1).
O algoritmo das aulas devolve o ponto final x_t; aqui devolve-se também o
melhor ponto visitado (xbest), que na prática convém guardar.

Otimização — deck 5.2 (Simulated annealing).
Reproduz os exemplos dos slides: ver ex05_2_simulated_annealing.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não adapta sigma nem T ao longo da corrida,
não tem reaquecimentos).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class SAResult:
    """Resultado de `simulated_annealing` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # ponto final x_t (o que o algoritmo das aulas devolve)
    fx: float                 # f(x_t)
    xbest: np.ndarray         # melhor ponto visitado
    fbest: float              # f(xbest)
    nit: int                  # pontos aceites (t final)
    nfev: int                 # avaliações de f (inclui f(x0))
    nacc: int                 # = nit
    T: float                  # temperatura no fim
    nT: int = 0               # reduções de temperatura feitas
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0             # 0: T < Ts; 1: atingiu nfmax
    message: str = ""


def simulated_annealing(f, x0, T0, alpha, n, Ts, sigma, lb=-np.inf, ub=np.inf,
                        rng=None, nfmax=100000):
    """Simulated annealing das aulas (patamares de n pontos aceites).

    Parâmetros
    ----------
    f : função objetivo de x (vetor numpy)
    x0 : ponto inicial (escalar ou vetor)
    T0 : temperatura inicial, T0 > 0 (na escala de f)
    alpha : fator de redução da temperatura, 0 < alpha < 1 (T = alpha*T)
    n : número de pontos aceites a cada temperatura
    Ts : temperatura mínima: termina quando T < Ts (verificado após cada
         aceitação, como no passo 4)
    sigma : desvio-padrão da perturbação gaussiana (vizinhança)
    lb, ub : limites de X (escalares ou vetores); o vizinho é projetado em X
    rng : gerador numpy (np.random.Generator) ou semente inteira; None: novo
        gerador. Para repetir uma corrida, passar np.random.default_rng(s).
    nfmax : limite de avaliações de f (salvaguarda; flag 1 se for atingido)

    Devolve
    -------
    SAResult com x = ponto final, xbest = melhor visitado, nfev (conta f(x0)),
    nit = nacc = t, T final, nT, e history com uma linha por avaliação:
    [n_f, t, T, aceite, f(x'), f(x_t), f_best, x'_1..x'_n, x_t,1..x_t,n]
    (a linha 0 é o ponto inicial; T é a temperatura usada nessa proposta).
    """
    if rng is None or np.isscalar(rng):            # semente -> gerador
        rng = np.random.default_rng(rng)
    x = np.atleast_1d(np.asarray(x0, float)).copy()
    nv = x.size
    hist = []

    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    fx = f(x); nfev = 1                                      # passo 1
    xbest = x.copy(); fbest = fx; T = T0; t = 0; nT = 0; flag = 1
    hist.append(np.concatenate(([1, 0, T0, 1, fx, fx, fbest], x, x)))
    while nfev < nfmax:
        xn = x + sigma * rng.normal(size=x.shape)            # passo 2: vizinho
        xn = np.minimum(np.maximum(xn, lb), ub)              # projetar em X
        fn = f(xn); df = fn - fx; nfev += 1
        Tk = T
        if df < 0 or rng.random() < np.exp(-df / T):         # passo 3: Metropolis
            x = xn; fx = fn; t += 1; aceite = 1
            if fx < fbest:
                xbest = x.copy(); fbest = fx
        else:
            aceite = 0                                       # rejeitado: passo 2
        hist.append(np.concatenate(([nfev, t, Tk, aceite, fn, fx, fbest], xn, x)))
        if aceite:                                           # passo 4
            if T < Ts:
                flag = 0
                break
            if t % n == 0:
                T = alpha * T; nT += 1
    # ----------------------------------------------------------------------

    cols = ("n_f", "t", "T", "aceite", "f(x')", "f(x_t)", "f_best") + \
        tuple("x'_%d" % (i + 1) for i in range(nv)) + tuple("x_t,%d" % (i + 1) for i in range(nv))
    if flag == 0:
        msg = "T = %.3g < Ts: %d aceites, %d avaliações" % (T, t, nfev)
    else:
        msg = "atingiu nfmax = %d avaliações (T = %.3g, %d aceites)" % (nfmax, T, t)
    return SAResult(x=x, fx=fx, xbest=xbest, fbest=fbest, nit=t, nfev=nfev, nacc=t, T=T,
                    nT=nT, history=np.array(hist), cols=cols, flag=flag, message=msg)
