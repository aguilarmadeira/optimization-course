"""Método de Nelder–Mead (simplex): min f(x) em R^n, só com valores de f.

Otimização — deck 3.2.3 (Método de Nelder–Mead).
Coeficientes (alpha, gamma, rho, sigma) = (1, 2, 1/2, 1/2) e regras de
aceitação/encolhimento do pseudocódigo do deck («O algoritmo»).
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
    x: np.ndarray             # melhor vértice do simplex final
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
    flag: int = 0
    message: str = ""


def nelder_mead(f, X0, tolx=1e-6, tolf=1e-6, kmax=500, ftarget=-np.inf, verbose=False):
    """Nelder–Mead a partir de um ponto ou de um simplex inicial.

    Parâmetros
    ----------
    f : função de x (vetor numpy de n componentes)
    X0 : ponto inicial x0 (n componentes) -> simplex x0, x0 + h_i e_i, com
         h_i = 0.05*max(1, |x0_i|) (o do deck); ou matriz (n+1) x n com os
         vértices do simplex inicial (nos exemplos do deck: arestas 1.2 no
         Himmelblau e 0.5 no Rosenbrock)
    tolx, tolf : para quando f(x_w) - f(x_b) <= tolf e o simplex tem
         diâmetro <= tolx (medido como max_i ||x_i - x_b||); 1e-6, 1e-6
    kmax : número máximo de iterações (500)
    ftarget : alvo para nhit (por omissão -inf: sem alvo)
    verbose : se True, imprime a tabela das iterações

    Devolve
    -------
    NelderMeadResult com x (melhor vértice), fx, nit, nfev (inclui as n+1
    avaliações do simplex inicial; 1-2 por iteração, n+2 com encolhimento),
    nhit (1.ª avaliação com f < ftarget), ops (operação de cada iteração),
    simplex final e history com uma linha por k = 0, ..., nit:
    [k, n_f, f(x_b), f(x_w), x_b].
    """
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5     # coeficientes originais
    X0 = np.array(X0, float)
    if X0.ndim == 1:                                  # simplex inicial do deck
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

    def linha(k, F):
        i = int(np.argmin(F)); return [k, cnt["nfev"], F[i], F.max(), *X[i]]

    # --- núcleo (o pseudocódigo do deck, com as contagens) ----------------
    F = np.array([aval(x) for x in X])                # n+1 avaliações
    hist = [linha(0, F)]; ops = ["inicial"]
    flag = 1
    for k in range(1, kmax + 1):
        i = np.argsort(F, kind="stable"); X = X[i]; F = F[i]   # ordenar
        # x_b = X[0], x_s = X[n-1], x_w = X[n]
        if F[n] - F[0] <= tolf and np.max(np.linalg.norm(X - X[0], axis=1)) <= tolx:
            flag = 0                                  # parou (antes de fazer a iteração k)
            break
        xw = X[n].copy()
        xc = X[:n].mean(axis=0)                       # centroide sem o pior
        xr = xc + alpha * (xc - xw); fr = aval(xr); op = "reflexão"
        encolher = False
        if fr < F[0]:                                 # expandir
            xe = xc + gamma * (xr - xc); fe = aval(xe)
            if fe < fr: X[n], F[n], op = xe, fe, "expansão"
            else: X[n], F[n] = xr, fr
        elif fr < F[n - 1]:                           # aceitar
            X[n], F[n] = xr, fr
        elif fr < F[n]:                               # contração exterior
            xoc = xc + rho * (xr - xc); foc = aval(xoc); op = "contr. ext."
            if foc <= fr: X[n], F[n] = xoc, foc
            else: encolher = True
        else:                                         # contração interior
            xic = xc - rho * (xc - xw); fic = aval(xic); op = "contr. int."
            if fic < F[n]: X[n], F[n] = xic, fic
            else: encolher = True
        if encolher:                                  # x_i <- x_b + sigma (x_i - x_b), i != b
            op = "encolhimento"
            for j in range(1, n + 1):
                X[j] = X[0] + sigma * (X[j] - X[0]); F[j] = aval(X[j])
        hist.append(linha(k, F)); ops.append(op)
    # ----------------------------------------------------------------------

    i = np.argsort(F, kind="stable"); X = X[i]; F = F[i]
    nit = len(hist) - 1
    msg = ("f(x_w) - f(x_b) <= tolf e diâmetro <= tolx" if flag == 0
           else "atingiu kmax = %d iterações" % kmax)
    cols = ("k", "n_f", "f(x_b)", "f(x_w)") + tuple("x_b%d" % (j + 1) for j in range(n))
    hist = np.array(hist)
    if verbose:
        print(("%4s %5s %12s %12s" + " %10s" * n + "  %s") % (cols + ("operação",)))
        for row, o in zip(hist, ops):
            print(("%4d %5d %12.4f %12.4f" + " %10.4f" * n + "  %s") % (*row, o))

    return NelderMeadResult(x=X[0].copy(), fx=F[0], nit=nit, nfev=cnt["nfev"],
                            nhit=cnt["nhit"], history=hist, ops=ops, cols=cols,
                            simplex=X, flag=flag, message=msg)
