"""Método de Newton em R^n, puro ou amortecido: H_k d_k = -g_k.

Otimização — deck 3.3.2 (Método de Newton).
Reproduz os exemplos dos slides: ver ex03_3_2_newton.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. pesquisa em linha de Armijo/Wolfe a começar em
alpha = 1; aqui o passo amortecido é o minimizante da linha).

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
    que g'*d no MATLAB/Octave, para que as duas versões deem os mesmos números)."""
    return float(np.sum(a * b))


DIRECOES = {0: "Newton", 1: "Newton com H + mu I", 2: "-g (salvaguarda)"}


@dataclass
class NewtonNDResult:
    """Resultado de `newton_nd` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # último iterando x_k
    fx: float                 # f(x_k) (nan se f is None)
    nit: int                  # iterações feitas
    nfev: int                 # avaliações de f: f nos iterandos (+ pesquisas em linha, se amortecido)
    ngev: int                 # gradientes: nit + 1 (o último é o do teste de paragem)
    nhev: int                 # Hessianas: nit
    nfev_ls: int = 0          # avaliações de f só nas pesquisas em linha
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0             # 0: ||g|| <= tolg; 1: atingiu kmax
    message: str = ""


def _definida_positiva(H):
    """Teste de Cholesky: True se H é (numericamente) definida positiva."""
    try:
        np.linalg.cholesky(H)
        return True
    except np.linalg.LinAlgError:
        return False


def newton_nd(f, grad, hess, x0, tolg=1e-8, kmax=50, damped=False, modify=False,
              ls="brent", lstol=None, verbose=False):
    """Newton em R^n a partir de x0.

    Parâmetros
    ----------
    f : função de R^n em R (None só no Newton puro: então n_f = 0)
    grad, hess : gradiente (vetor) e Hessiana (matriz n x n) de f
    x0 : ponto inicial
    tolg : pára quando ||g_k|| <= tolg (flag 0)
    kmax : número máximo de iterações (flag 1, com aviso)
    damped : False (por omissão) = Newton puro, passo completo alpha = 1;
             True = amortecido: se g^T d >= 0 usa d = -g (salvaguarda) e
             alpha por pesquisa em linha
    modify : (só amortecido) se True e H_k não for definida positiva
             (Cholesky), usa H_k + mu I com mu = 1e-3 max(1, ||H_k||_F),
             10 mu, 100 mu, ... até ser definida positiva
    ls, lstol : pesquisa em linha do amortecido ("brent", tol 1e-10, a das
             figuras; ou "fminbnd") -- ver line_search.py
    verbose : se True, imprime a tabela das iterações

    Contagens: n_g = nit + 1, n_H = nit; n_f = f em cada iterando (x_0
    incluído) + as avaliações das pesquisas em linha (amortecido).
    history tem uma linha por k = 0, ..., nit:
    [k, x_k (n colunas), f(x_k), ||g_k||, alpha_{k-1}, lambda_min(H_{k-1}),
     direção_{k-1} (0 Newton, 1 H + mu I, 2 -g)]   (nan se não se aplica).
    O cálculo de lambda_min é só para leitura (não é usado pelo método).
    """
    if damped and f is None:
        raise ValueError("newton_nd: o Newton amortecido precisa de f")
    temf = f is not None
    x = np.array(x0, float).ravel()
    n = x.size
    nfev = ngev = nhev = nfev_ls = 0
    if temf:
        fx = f(x); nfev += 1                       # f no iterando x_0
    else:
        fx = np.nan
    rows = [[0, *x, fx, np.nan, np.nan, np.nan, np.nan]]
    flag = 1
    k = 0
    # --- núcleo (o dos slides, com as contagens) ---------------------------
    while True:
        if k >= kmax:
            flag = 1; break
        g = np.asarray(grad(x), float).ravel(); ngev += 1
        rows[k][n + 2] = np.linalg.norm(g)
        if np.linalg.norm(g) <= tolg:
            flag = 0; break
        H = np.asarray(hess(x), float); nhev += 1
        lmin = np.linalg.eigvalsh(H).min()        # só para leitura
        tipo = 0
        if damped and modify and not _definida_positiva(H):
            mu = 1e-3 * max(1.0, np.linalg.norm(H))
            while not _definida_positiva(H):
                H = H + mu * np.eye(n); mu *= 10
            tipo = 1
        d = -np.linalg.solve(H, g)                # H d = -g (sem inverter H)
        if damped:
            if _dot(g, d) >= 0:
                d = -g; tipo = 2                  # salvaguarda: não é de descida
            phi = lambda a: f(x + a * d)          # pesquisa em linha
            alpha, nls, _ = line_search(phi, ls, lstol)
            nfev += nls; nfev_ls += nls
        else:
            alpha = 1.0                           # passo completo
        x = x + alpha * d
        if temf:
            fx = f(x); nfev += 1                  # f no novo iterando
        k += 1
        rows.append([k, *x, fx, np.nan, alpha, lmin, tipo])
    # ----------------------------------------------------------------------

    hist = np.array(rows, float)
    cols = ("k",) + tuple("x%d" % (i + 1) for i in range(n)) + \
           ("f(x_k)", "||g_k||", "alpha_k-1", "lmin(H_k-1)", "direcao_k-1")
    if flag == 0:
        message = "||g|| <= tolg em %d iterações" % k
    else:
        message = "atingiu o número máximo de iterações (kmax = %d) sem ||g|| <= tolg" % kmax
        import warnings
        warnings.warn("newton_nd: " + message)

    if verbose:
        print("%4s %24s %16s %11s %9s  %s" % ("k", "x_k", "f(x_k)", "||g_k||", "alpha", "direção"))
        for r in hist:
            xs = "; ".join("%9.4f" % v for v in r[1:n + 1])
            dirs = DIRECOES.get(int(r[n + 5]), "") if not np.isnan(r[n + 5]) else ""
            print("%4d  (%s) %16.6e %11.2e %9.4f  %s" % (r[0], xs, r[n + 1], r[n + 2], r[n + 3], dirs))

    return NewtonNDResult(x=x, fx=fx, nit=k, nfev=nfev, ngev=ngev, nhev=nhev,
                          nfev_ls=nfev_ls, history=hist, cols=cols, flag=flag,
                          message=message)
