"""Gradientes conjugados não lineares (Fletcher-Reeves) com reinícios.

Otimização — deck 3.3.3 (Gradientes conjugados).
Reproduz os exemplos dos slides: ver ex03_3_3_conjugate_gradients.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. pesquisa em linha com condições de Wolfe forte).

Importa `line_search` de code/python/common/. Os exemplos ex*.py tratam do caminho
(import uc_setup); fora deles, é preciso primeiro
    sys.path.insert(0, "<repo>/code/python"); import uc_setup

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np

from line_search import line_search


def _dot(a, b):
    """Produto interno a^T b como soma dos produtos (a mesma ordem de operações
    que g'*g no MATLAB/Octave, para que as duas versões deem os mesmos números)."""
    return float(np.sum(a * b))


@dataclass
class CGResult:
    """Resultado de `conj_grad` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # última iterada x_k
    fx: float                 # f(x_k)
    nit: int                  # iterações feitas
    nfev: int                 # avaliações de f: pesquisas em linha + f nas iteradas
    ngev: int                 # gradientes: nit + 1
    nhev: int = 0             # o método não usa a Hessiana
    nfev_ls: int = 0          # avaliações de f só nas pesquisas em linha
    nrestart: int = 0         # reinícios feitos (d = -g depois de x_0)
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0             # 0: ||g|| <= tolg; 1: atingiu kmax
    message: str = ""


def conj_grad(f, grad, x0, tolg=1e-6, kmax=1000, restart=None, ls="brent",
              lstol=None, verbose=False):
    """Fletcher-Reeves: d_0 = -g_0, d_{k+1} = -g_{k+1} + beta_{k+1} d_k,
    beta_{k+1} = ||g_{k+1}||^2 / ||g_k||^2.

    Notação das aulas e dos slides: S_k = direção (aqui d), lambda_k =
    comprimento do passo (aqui alpha), grad f_k = gradiente (aqui g).

    Parâmetros
    ----------
    f : função de R^n em R
    grad : gradiente de f (função que devolve um vetor)
    x0 : ponto inicial
    tolg : pára quando ||g_k|| <= tolg (flag 0)
    kmax : número máximo de iterações (flag 1)
    restart : reinício periódico, d = -g quando mod(k+1, restart) = 0; por
              omissão restart = n (teste 1 do slide); 0 = sem reinício
              periódico.  O reinício quando d deixa de ser de descida
              (g^T d >= 0, teste 2) está sempre ativo.
    ls, lstol : pesquisa em linha, "brent" (tol 1e-10, a das figuras) ou
              "fminbnd" (TolX 1e-4) -- ver line_search.py
    verbose : se True, imprime a tabela das iterações

    Contagens: n_g = nit + 1; n_f = avaliações das pesquisas em linha (cada
    uma volta a avaliar phi(0) = f(x_k) no modo "brent") + f em cada
    iterada (x_0 incluído).
    history tem uma linha por k = 0, ..., nit:
    [k, x_k (n colunas), f(x_k), ||g_k||, alpha_{k-1}, beta_k, d_k (n colunas)]
    (beta_k e d_k da última linha não são calculados: nan).
    """
    x = np.array(x0, float).ravel()
    n = x.size
    if restart is None:
        restart = n
    nfev = ngev = nfev_ls = nrestart = 0
    fx = f(x); nfev += 1                          # f na iterada x_0
    rows = []
    # --- núcleo (o dos slides, com as contagens) ---------------------------
    g = np.asarray(grad(x), float).ravel(); ngev += 1
    d = -g                                        # 1.a direção
    rows.append([0, *x, fx, np.linalg.norm(g), np.nan, np.nan, *d])
    k = 0
    flag = 0 if np.linalg.norm(g) <= tolg else 1
    while flag == 1 and k < kmax:
        phi = lambda a: f(x + a * d)              # pesquisa em linha
        alpha, nls, _ = line_search(phi, ls, lstol)
        nfev += nls; nfev_ls += nls
        x = x + alpha * d
        fx = f(x); nfev += 1                      # f na nova iterada
        gn = np.asarray(grad(x), float).ravel(); ngev += 1
        rows.append([k + 1, *x, fx, np.linalg.norm(gn), alpha, np.nan, *np.full(n, np.nan)])
        if np.linalg.norm(gn) <= tolg:
            flag = 0; k += 1; break
        beta = _dot(gn, gn) / _dot(g, g)          # Fletcher-Reeves
        d = -gn + beta * d
        if (restart > 0 and (k + 1) % restart == 0) or _dot(gn, d) >= 0:
            d = -gn; beta = 0.0; nrestart += 1    # reinício
        rows[-1][n + 4] = beta
        rows[-1][n + 5:] = d
        g = gn; k += 1
    # ----------------------------------------------------------------------

    hist = np.array(rows, float)
    cols = ("k",) + tuple("x%d" % (i + 1) for i in range(n)) + \
           ("f(x_k)", "||g_k||", "alpha_k-1", "beta_k") + tuple("d%d" % (i + 1) for i in range(n))
    if flag == 0:
        message = "||g|| <= tolg em %d iterações" % k
    else:
        message = "atingiu o número máximo de iterações (kmax = %d) sem ||g|| <= tolg" % kmax

    if verbose:
        print("%4s %22s %14s %11s %9s %8s %22s" % ("k", "x_k", "f(x_k)", "||g_k||", "alpha", "beta", "d_k"))
        for r in hist:
            xs = "; ".join("%9.4f" % v for v in r[1:n + 1])
            ds = "; ".join("%9.4f" % v for v in r[n + 5:])
            print("%4d  (%s) %14.6f %11.2e %9.4f %8.4f  (%s)" % (r[0], xs, r[n + 1], r[n + 2], r[n + 3], r[n + 4], ds))

    return CGResult(x=x, fx=fx, nit=k, nfev=nfev, ngev=ngev, nfev_ls=nfev_ls,
                    nrestart=nrestart, history=hist, cols=cols, flag=flag, message=message)
