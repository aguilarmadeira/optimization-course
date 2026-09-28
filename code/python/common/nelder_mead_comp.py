"""Nelder–Mead da comparação de 4.3.2 (minimizador interno dos exemplos do cap. 4.3).

Otimização — decks 4.3.1 e 4.3.2 (penalização exterior e barreira).
Utilitário comum (code/python/common/), usado por penalty_exterior (4.3.1),
barrier_log (4.3.2) e pelos dois exemplos. É a VARIANTE do Nelder–Mead usada
na comparação do 4.3: a contração exterior só é aceite com «<» (o
`nelder_mead` do deck 3.2.3 aceita com «<=»). Existe como função à parte, e
não como opção de `nelder_mead`, para que o `nelder_mead` do deck 3.2.3
fique exatamente o do pseudocódigo dos slides.

Variante de `nelder_mead` (deck 3.2.3, pasta
03_multivariate_optimization/03_2_3_nelder_mead/), para reproduzir EXATAMENTE
o Nelder–Mead didático do script da comparação exterior vs. interior
(`make_figs_4comp.py`), que produziu os números do slide «Exterior vs.
interior» (844, 743 chamadas a P, n_f = 717, 26 fora de X; 864 de (1,1)).

O algoritmo é o do deck 3.2.3 — coeficientes (alpha, gamma, rho, sigma) =
(1, 2, 1/2, 1/2), simplex inicial x0, x0 + h_i e_i com h_i = 0.05*max(1,|x0_i|),
mesma paragem — com duas diferenças, ambas do script da comparação:
  1. contração exterior aceite só se f(xoc) < f(xr) (desigualdade ESTRITA;
     o pseudocódigo do deck 3.2.3 aceita com <=). Perto da convergência,
     com tolx = 1e-10, há empates exatos f(xoc) == f(xr) e a regra de
     desempate muda o número de avaliações;
  2. os pontos calculam-se como no script: xe = xc + 2(xc - xw),
     xoc = xc + (xc - xw)/2 (em aritmética exata é o mesmo que
     xc + gamma(xr - xc) e xc + rho(xr - xc); os arredondamentos diferem).
Tolerâncias por omissão: as da comparação (tolx = 1e-10, tolf = 1e-12, kmax = 2000).
Aceita f = inf (pontos fora de X na barreira): esses vértices nunca são aceites
quando há alternativa finita.

Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não deteta simplex degenerado nem reinicia).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class NelderMeadCompResult:
    """Resultado de `nelder_mead_comp` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # melhor vértice do simplex final
    fx: float                 # f(x)
    nit: int                  # iterações feitas
    nfev: int                 # avaliações de f (inclui as n+1 do simplex inicial)
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0
    message: str = ""


def nelder_mead_comp(f, x0, tolx=1e-10, tolf=1e-12, kmax=2000):
    """Nelder–Mead da comparação de 4.3.2, a partir de um ponto x0.

    Parâmetros
    ----------
    f : função de x (vetor numpy de n componentes); pode devolver inf
    x0 : ponto inicial -> simplex x0, x0 + h_i e_i, h_i = 0.05*max(1, |x0_i|)
    tolx, tolf : para quando f(x_w) - f(x_b) <= tolf e max_i ||x_i - x_b|| <= tolx
    kmax : número máximo de iterações

    Devolve
    -------
    NelderMeadCompResult com x, fx, nit, nfev (inclui as n+1 avaliações do
    simplex inicial) e history com uma linha por k = 0, ..., nit:
    [k, n_f, f(x_b), f(x_w), x_b].
    """
    x0 = np.array(x0, float).ravel()
    n = x0.size
    h = 0.05 * np.maximum(1, np.abs(x0))
    X = np.vstack([x0, x0 + np.diag(h)])
    nfev = 0

    # --- núcleo (deck 3.2.3, com as duas diferenças do script) ------------
    F = np.array([f(x) for x in X], float); nfev += n + 1
    hist = [[0, nfev, F.min(), F.max(), *X[int(np.argmin(F))]]]
    flag = 1
    for k in range(1, kmax + 1):
        i = np.argsort(F, kind="stable"); X = X[i]; F = F[i]   # ordenar
        if F[n] - F[0] <= tolf and np.max(np.linalg.norm(X - X[0], axis=1)) <= tolx:
            flag = 0
            break
        xw = X[n].copy()
        xc = X[:n].mean(axis=0)                     # centroide sem o pior
        xr = xc + (xc - xw); fr = f(xr); nfev += 1  # reflexão
        if fr < F[0]:                               # expandir
            xe = xc + 2 * (xc - xw); fe = f(xe); nfev += 1
            if fe < fr: X[n], F[n] = xe, fe
            else: X[n], F[n] = xr, fr
        elif fr < F[n - 1]:                         # aceitar
            X[n], F[n] = xr, fr
        else:                                       # contrair
            if fr < F[n]:
                xt = xc + 0.5 * (xc - xw)           # contração exterior
            else:
                xt = xc - 0.5 * (xc - xw)           # contração interior
            ft = f(xt); nfev += 1
            if ft < min(fr, F[n]):                  # estrito (script)
                X[n], F[n] = xt, ft
            else:                                   # encolher
                for j in range(1, n + 1):
                    X[j] = X[0] + 0.5 * (X[j] - X[0]); F[j] = f(X[j])
                nfev += n
        hist.append([k, nfev, F.min(), F.max(), *X[int(np.argmin(F))]])
    # ----------------------------------------------------------------------

    i = np.argsort(F, kind="stable"); X = X[i]; F = F[i]
    msg = ("f(x_w) - f(x_b) <= tolf e diâmetro <= tolx" if flag == 0
           else "atingiu kmax = %d iterações" % kmax)
    cols = ("k", "n_f", "f(x_b)", "f(x_w)") + tuple("x_b%d" % (j + 1) for j in range(n))
    return NelderMeadCompResult(x=X[0].copy(), fx=F[0], nit=len(hist) - 1, nfev=nfev,
                                history=np.array(hist), cols=cols, flag=flag, message=msg)
