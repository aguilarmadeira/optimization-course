"""Agregação de objetivos: soma ponderada e Tchebycheff ponderado.

weighted_sum:          para cada vetor de pesos w,  min sum_i w_i f_i(x)
weighted_tchebycheff:  para cada w,  min theta  s.a.  theta - w_i (f_i(x) - z_i^id) >= 0
                       (a reformulação sem o max, da nota do slide)

Cada problema escalar é resolvido com scipy.optimize.minimize(method='SLSQP'),
como no slide (a UC não reimplementa o SLSQP); o resto usa só numpy.

Otimização — deck 6.3 (Métodos de agregação de objetivos).
Reproduz o exemplo dos slides: ver ex06_3_aggregation.py
Implementação didática — cada problema escalar é resolvido localmente; só um
minimizante global da soma ponderada tem garantia de ser eficiente.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize


@dataclass
class AggResult:
    """Resultado de `weighted_sum` / `weighted_tchebycheff` (campos de `info` no MATLAB)."""
    x: np.ndarray             # um ponto por linha (nw x n)
    fx: np.ndarray            # F(x) de cada ponto (nw x m)
    nit: int                  # soma das iterações do SLSQP
    nfev: int                 # chamadas a F (soma ponderada: = soma de r.nfev; inclui diferenças finitas)
    ngev: int = 0             # os gradientes são numéricos, dentro do SLSQP
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    status: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    flag: int = 0
    message: str = ""


def _pontos(X0, n):
    X0 = np.asarray(X0, float)
    return X0.reshape(-1, n)


def _tabela(cols, hist, m, fmt_w):
    print("".join(" %7s" % c for c in cols))
    for h in hist:
        print("".join(fmt_w % v for v in h[:m]) + "".join(" %7.4f" % v for v in h[m:-1])
              + " %7d" % h[-1])


def _saida(cols, X, F, hist, status, nit, nfev, nw):
    falhou = bool(np.any(status != 0))
    msg = ("SLSQP: algum problema não terminou bem (ver status)" if falhou
           else "SLSQP: %d problemas escalares resolvidos" % nw)
    return AggResult(x=X, fx=F, nit=nit, nfev=nfev, history=hist, cols=cols,
                     status=status, flag=int(falhou), message=msg)


def weighted_sum(F, X0, bounds, W, cons=(), verbose=False):
    """Soma ponderada: para cada linha w de W, min w @ F(x) s.a. cons, bounds.

    Parâmetros
    ----------
    F : função que devolve o vetor [f1(x), ..., fm(x)]
    X0 : pontos iniciais, um por linha:
         - um só ponto: arranque a quente (a solução para um w é o ponto
           inicial do w seguinte), como no slide;
         - vários: para cada w parte de todos e fica com o de menor soma
           ponderada (sem arranque a quente); útil em problemas não convexos
    bounds : limites das variáveis, [(l1, u1), ...]
    W : pesos, um vetor por linha (nw x m), w_i >= 0, sum_i w_i = 1
    cons : restrições no formato do SLSQP ({'type': 'ineq', 'fun': g}, g >= 0)
    verbose : se True, imprime a tabela dos pesos

    Devolve
    -------
    AggResult com x, fx, nit, nfev (chamadas a F = soma de r.nfev) e history
    com uma linha por w: [w, x, F(x), n_f(w)].
    """
    W = np.atleast_2d(np.asarray(W, float)); nw, m = W.shape
    n = len(bounds); X0 = _pontos(X0, n)
    X = np.zeros((nw, n)); Fx = np.zeros((nw, m))
    hist = np.zeros((nw, 2 * m + n + 1)); status = np.zeros(nw, int)
    cont = [0]

    def Fc(z):
        cont[0] += 1
        return np.asarray(F(z), float)

    nit = 0; x0 = X0[0]
    for k, w in enumerate(W):
        antes = cont[0]
        melhor = None
        for s in range(len(X0)):
            if len(X0) > 1:
                x0 = X0[s]
            r = minimize(lambda z: float(w @ Fc(z)), x0, bounds=bounds,
                         constraints=list(cons), method="SLSQP")
            nit += r.nit
            if melhor is None or r.fun < melhor.fun:
                melhor = r
        x0 = melhor.x                                   # a quente (só com um ponto inicial)
        X[k] = melhor.x; Fx[k] = F(melhor.x); status[k] = melhor.status
        hist[k] = [*w, *melhor.x, *Fx[k], cont[0] - antes]

    cols = (tuple("w%d" % (i + 1) for i in range(m)) + tuple("x%d" % (i + 1) for i in range(n))
            + tuple("f%d" % (i + 1) for i in range(m)) + ("nf",))
    if verbose:
        _tabela(cols, hist, m, " %7.2f")
    return _saida(cols, X, Fx, hist, status, nit, cont[0], nw)


def weighted_tchebycheff(F, X0, bounds, W, zid, cons=(), verbose=False):
    """Tchebycheff ponderado: min max_i w_i (f_i(x) - z_i^id), pela reformulação

        min theta   s.a.   theta - w_i (f_i(x) - z_i^id) >= 0,  i = 1..m,  cons,  bounds

    em z = [x, theta] (theta livre). Parâmetros como em weighted_sum, mais o
    ponto ideal zid; theta0 = max_i w_i (f_i(x0) - z_i^id).

    Devolve AggResult; history = [w, x, F(x), theta, n_f(w)] (n_f = chamadas a
    F, que aqui aparece nas restrições).
    """
    W = np.atleast_2d(np.asarray(W, float)); nw, m = W.shape
    n = len(bounds); X0 = _pontos(X0, n); zid = np.asarray(zid, float)
    X = np.zeros((nw, n)); Fx = np.zeros((nw, m))
    hist = np.zeros((nw, 2 * m + n + 2)); status = np.zeros(nw, int)
    cont = [0]

    def Fc(z):
        cont[0] += 1
        return np.asarray(F(z), float)

    bz = list(bounds) + [(None, None)]
    cons_z = [{"type": c["type"], "fun": (lambda z, g=c["fun"]: g(z[:n]))} for c in cons]
    nit = 0; x0 = X0[0]
    for k, w in enumerate(W):
        gz = {"type": "ineq", "fun": lambda z, w=w: z[n] - w * (Fc(z[:n]) - zid)}
        antes = cont[0]
        melhor = None
        for s in range(len(X0)):
            if len(X0) > 1:
                x0 = X0[s]
            z0 = np.r_[x0, np.max(w * (np.asarray(F(x0)) - zid))]
            r = minimize(lambda z: z[n], z0, bounds=bz, constraints=[gz] + cons_z,
                         method="SLSQP")
            nit += r.nit
            if melhor is None or r.fun < melhor.fun:
                melhor = r
        x0 = melhor.x[:n]                               # a quente
        X[k] = melhor.x[:n]; Fx[k] = F(melhor.x[:n]); status[k] = melhor.status
        hist[k] = [*w, *melhor.x[:n], *Fx[k], melhor.x[n], cont[0] - antes]

    cols = (tuple("w%d" % (i + 1) for i in range(m)) + tuple("x%d" % (i + 1) for i in range(n))
            + tuple("f%d" % (i + 1) for i in range(m)) + ("theta", "nf"))
    if verbose:
        _tabela(cols, hist, m, " %7.3f")
    return _saida(cols, X, Fx, hist, status, nit, cont[0], nw)
