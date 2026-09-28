"""Método do gradiente (steepest descent, método de Cauchy) com pesquisa em linha.

Otimização — deck 3.3.1 (Método do gradiente).
Reproduz os exemplos dos slides: ver ex03_3_1_steepest_descent.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. pesquisa em linha inexata com condições de Wolfe).

Importa `line_search` de code/python/common/. Os exemplos ex*.py tratam do caminho
(import uc_setup); fora deles, é preciso primeiro
    sys.path.insert(0, "<repo>/code/python"); import uc_setup

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np

from grad_fd import grad_fd
from line_search import line_search


@dataclass
class SDResult:
    """Resultado de `steepest_descent` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # último iterando x_k
    fx: float                 # f(x_k)
    nit: int                  # iterações (passos) feitas
    nfev: int                 # avaliações de f: pesquisas em linha + f nos iterandos (+ 2n por gradiente numérico)
    ngev: int                 # gradientes analíticos (0 se numérico)
    nhev: int = 0             # o método não usa a Hessiana
    nfev_ls: int = 0          # avaliações de f só nas pesquisas em linha
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0             # 0: ||g|| <= tolg; 1: atingiu kmax; 2: passo relativo <= tolx
    message: str = ""


def steepest_descent(f, grad, x0, tolg=1e-6, kmax=1000, tolx=0.0, ls="brent",
                     lstol=None, verbose=False):
    """Método do gradiente: x_{k+1} = x_k - alpha_k g_k, alpha_k por pesquisa em linha.

    Parâmetros
    ----------
    f : função de R^n em R
    grad : gradiente de f (função que devolve um vetor); None usa grad_fd
           (diferenças centrais: 2n avaliações de f por gradiente, em nfev)
    x0 : ponto inicial
    tolg : pára quando ||g_k|| <= tolg (flag 0)
    kmax : número máximo de iterações (flag 1; o gradiente em x_kmax não é
           calculado, por isso n_g = kmax)
    tolx : pára quando ||x_{k+1} - x_k|| / max(1, ||x_k||) <= tolx (flag 2);
           por omissão 0 (só pára se o ponto não mudar)
    ls, lstol : pesquisa em linha, "brent" (tol 1e-10, a das figuras) ou
           "fminbnd" (TolX 1e-4) -- ver line_search.py
    verbose : se True, imprime a tabela das iterações

    Contagens: n_g = nit + 1 (um gradiente por iteração mais o de paragem);
    n_f = avaliações das pesquisas em linha (cada uma volta a avaliar
    phi(0) = f(x_k) no modo "brent") + f em cada iterando (x_0 incluído).
    history tem uma linha por k = 0, ..., nit:
    [k, x_k (n colunas), f(x_k), ||g_k||, alpha_{k-1}]  (nan se não calculado).
    """
    x = np.array(x0, float).ravel()
    n = x.size
    nfev = 0; ngev = 0; nfev_ls = 0
    fx = f(x); nfev += 1                     # f no iterando x_0
    rows = [[0, *x, fx, np.nan, np.nan]]
    flag = 1
    k = 0
    # --- núcleo (o dos slides, com as contagens) ---------------------------
    while True:
        if k >= kmax:
            flag = 1; break
        if grad is None:
            g, nf = grad_fd(f, x); nfev += nf     # 2n avaliações de f
        else:
            g = np.asarray(grad(x), float).ravel(); ngev += 1
        rows[k][n + 2] = np.linalg.norm(g)
        if np.linalg.norm(g) <= tolg:
            flag = 0; break
        d = -g                                    # direção
        phi = lambda a: f(x + a * d)              # pesquisa em linha
        alpha, nls, _ = line_search(phi, ls, lstol)
        nfev += nls; nfev_ls += nls
        xn = x + alpha * d
        fx = f(xn); nfev += 1                     # f no novo iterando
        passo = np.linalg.norm(xn - x) / max(1.0, np.linalg.norm(x))
        x = xn; k += 1
        rows.append([k, *x, fx, np.nan, alpha])
        if passo <= tolx:
            flag = 2; break
    # ----------------------------------------------------------------------

    hist = np.array(rows, float)
    cols = ("k",) + tuple("x%d" % (i + 1) for i in range(n)) + ("f(x_k)", "||g_k||", "alpha_k-1")
    if flag == 0:
        message = "||g|| <= tolg em %d iterações" % k
    elif flag == 1:
        message = "atingiu o número máximo de iterações (kmax = %d) sem ||g|| <= tolg" % kmax
    else:
        message = "passo relativo <= tolx em %d iterações" % k

    if verbose:
        print("%4s %22s %14s %12s %9s" % ("k", "x_k", "f(x_k)", "||g_k||", "alpha"))
        for r in hist:
            xs = "; ".join("%9.4f" % v for v in r[1:n + 1])
            print("%4d  (%s) %14.5f %12.4e %9.4f" % (r[0], xs, r[n + 1], r[n + 2], r[n + 3]))

    return SDResult(x=x, fx=fx, nit=k, nfev=nfev, ngev=ngev, nfev_ls=nfev_ls,
                    history=hist, cols=cols, flag=flag, message=message)
