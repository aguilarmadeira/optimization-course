"""Método do simplex de Nelder–Mead: min f(x) em R^n, só com valores de f.

Otimização — deck 3.2.3 (Método do simplex de Nelder–Mead).
Notação das aulas: x_l = melhor vértice, x_g = o seguinte ao pior, x_h = pior;
x_c = centroide dos vértices exceto x_h; reflexão x_r = 2 x_c - x_h;
expansão (1+gamma) x_c - gamma x_h; contração interior (1-beta) x_c + beta x_h;
contração exterior (1+beta) x_c - beta x_h.  Por omissão gamma = 2, beta = 0.5.

Duas formas de parar:
  * eps (a das aulas): análise do erro Q = sqrt( sum_i (f(x_i) - f(x_c))^2 / (n+1) ),
    calculada em cada iteração com o novo simplex e o centroide dessa
    iteração (f(x_c) conta como uma avaliação); para quando Q <= eps;
  * tolx, tolf (Nelder–Mead padrão, como o fminsearch): para quando
    f(x_h) - f(x_l) <= tolf e max_i ||x_i - x_l|| <= tolx.
Regras nos casos-limite (regra = "padrao", por omissão): a expansão só fica
se f(x_e) < f(x_r) (senão fica x_r); uma contração só fica se melhorar
(exterior: f <= f(x_r); interior: f < f(x_h)); senão encolhe-se o simplex
para metade em direção a x_l.  Com regra = "deb" (Deb, 2012, sec. 3.3.2) o
ponto novo substitui sempre x_h e não há encolhimento.
Reproduz os exemplos dos slides: ver ex03_2_3_nelder_mead.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não deteta simplex degenerado nem reinicia).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class NelderMeadResult:
    """Resultado de `nelder_mead` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # melhor vértice do simplex final (x_l)
    fx: float                 # f(x)
    nit: int                  # iterações feitas
    nfev: int                 # avaliações de f (inclui as n+1 do simplex inicial)
    nhit: float = np.inf      # 1.ª avaliação com f < ftarget (inf se nunca)
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    ops: list = field(default_factory=list)   # operação de cada iteração
    cols: tuple = ()
    simplex: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    Q: float = np.nan         # último valor de Q (só com eps)
    flag: int = 0
    message: str = ""


def nelder_mead(f, X0, tolx=1e-6, tolf=1e-6, kmax=500, ftarget=-np.inf, verbose=False,
                gamma=2.0, beta=0.5, eps=None, regra="padrao"):
    """Simplex de Nelder–Mead a partir de um ponto ou de um simplex inicial.

    Parâmetros
    ----------
    f : função de x (vetor numpy de n componentes)
    X0 : ponto inicial x0 (n componentes) -> simplex x0, x0 + h_i e_i, com
         h_i = 0.05*max(1, |x0_i|); ou matriz (n+1) x n com os vértices do
         simplex inicial (nos exemplos do deck: arestas 1.2 no Himmelblau e
         0.5 no Rosenbrock)
    tolx, tolf : paragem do Nelder–Mead padrão (só se eps for None):
         f(x_h) - f(x_l) <= tolf e max_i ||x_i - x_l|| <= tolx; 1e-6, 1e-6
    kmax : número máximo de iterações (500)
    ftarget : alvo para nhit (por omissão -inf: sem alvo)
    verbose : se True, imprime a tabela das iterações
    gamma, beta : parâmetros de expansão e de contração (2 e 0.5)
    eps : se for dado, pára quando Q <= eps (a análise do erro das aulas)
    regra : "padrao" (por omissão) ou "deb" (ver o cabeçalho do módulo)

    Devolve
    -------
    NelderMeadResult com x (= x_l), fx, nit, nfev (inclui as n+1 avaliações
    do simplex inicial; 1-2 por iteração, mais 1 para f(x_c) quando se usa
    eps, e n mais com encolhimento), nhit, ops (operação de cada iteração),
    simplex final, Q e history com uma linha por k = 0, ..., nit:
    [k, n_f, f(x_l), f(x_h), Q, x_l] (Q = nan sem eps).
    """
    if regra not in ("padrao", "deb"):
        raise ValueError('nelder_mead: regra tem de ser "padrao" ou "deb".')
    sigma = 0.5                                       # encolhimento para metade
    X0 = np.array(X0, float)
    if X0.ndim == 1:                                  # simplex inicial a partir de x0
        x0 = X0; n = x0.size; h = 0.05 * np.maximum(1, np.abs(x0))
        X = np.vstack([x0, x0 + np.diag(h)])
    else:
        X = X0.copy(); n = X.shape[1]
        if X.shape[0] != n + 1:
            raise ValueError("nelder_mead: o simplex tem de ter n+1 vértices.")

    cnt = {"nfev": 0, "nhit": np.inf}

    def aval(x):                                      # f com as contagens
        v = f(x); cnt["nfev"] += 1
        if v < ftarget and cnt["nhit"] == np.inf:
            cnt["nhit"] = cnt["nfev"]
        return v

    def linha(k, F, Q):
        i = int(np.argmin(F)); return [k, cnt["nfev"], F[i], F.max(), Q, *X[i]]

    # --- núcleo (o algoritmo do deck, com as contagens) --------------------
    F = np.array([aval(x) for x in X])                # passo 1: n+1 avaliações
    hist = [linha(0, F, np.nan)]; ops = ["inicial"]
    flag = 1; Q = np.nan
    for k in range(1, kmax + 1):
        i = np.argsort(F, kind="stable"); X = X[i]; F = F[i]   # passo 2: ordenar
        # x_l = X[0], x_g = X[n-1], x_h = X[n]
        if eps is None and F[n] - F[0] <= tolf and np.max(np.linalg.norm(X - X[0], axis=1)) <= tolx:
            flag = 0                                  # paragem padrão (antes da iteração k)
            break
        xh = X[n].copy()
        xc = X[:n].mean(axis=0)                       # centroide sem o pior
        if eps is not None:
            fc = aval(xc)                             # f(x_c), para Q
        xr = xc + (xc - xh); fr = aval(xr); op = "reflexão"   # passo 3: x_r = 2 x_c - x_h
        encolher = False
        if fr < F[0]:                                 # expansão
            xe = xc + gamma * (xr - xc); fe = aval(xe)   # = (1+gamma) x_c - gamma x_h
            if regra == "deb" or fe < fr: X[n], F[n], op = xe, fe, "expansão"
            else: X[n], F[n] = xr, fr
        elif fr < F[n - 1]:                           # aceitar a reflexão
            X[n], F[n] = xr, fr
        elif fr < F[n]:                               # contração exterior
            xo = xc + beta * (xr - xc); fo = aval(xo); op = "contr. ext."   # = (1+beta) x_c - beta x_h
            if regra == "deb" or fo <= fr: X[n], F[n] = xo, fo
            else: encolher = True
        else:                                         # contração interior
            xi = xc - beta * (xc - xh); fi = aval(xi); op = "contr. int."   # = (1-beta) x_c + beta x_h
            if regra == "deb" or fi < F[n]: X[n], F[n] = xi, fi
            else: encolher = True
        if encolher:                                  # x_i = x_l + sigma (x_i - x_l), i != l
            op = "encolhimento"
            for j in range(1, n + 1):
                X[j] = X[0] + sigma * (X[j] - X[0]); F[j] = aval(X[j])
        if eps is not None:                           # passo 4: análise do erro
            Q = float(np.sqrt(np.mean((F - fc) ** 2)))
        hist.append(linha(k, F, Q)); ops.append(op)
        if eps is not None and Q <= eps:
            flag = 0
            break
    # ----------------------------------------------------------------------

    i = np.argsort(F, kind="stable"); X = X[i]; F = F[i]
    nit = len(hist) - 1
    if flag == 0:
        msg = "Q <= eps" if eps is not None else "f(x_h) - f(x_l) <= tolf e diâmetro <= tolx"
    else:
        msg = "atingiu kmax = %d iterações" % kmax
    cols = ("k", "n_f", "f(x_l)", "f(x_h)", "Q") + tuple("x_l%d" % (j + 1) for j in range(n))
    hist = np.array(hist)
    if verbose:
        print(("%4s %5s %12s %12s %10s" + " %10s" * n + "  %s") % (cols + ("operação",)))
        for row, o in zip(hist, ops):
            print(("%4d %5d %12.4f %12.4f %10.4f" + " %10.4f" * n + "  %s") % (*row, o))

    return NelderMeadResult(x=X[0].copy(), fx=F[0], nit=nit, nfev=cnt["nfev"],
                            nhit=cnt["nhit"], history=hist, ops=ops, cols=cols,
                            simplex=X, Q=Q, flag=flag, message=msg)
