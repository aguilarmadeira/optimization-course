"""Método do epsilon-constrangimento (dois objetivos).

Para cada valor e de E resolve o subproblema do cap. 4

     min f2(x)   s.a.   g(x) = e - f1(x) >= 0,   bounds

com scipy.optimize.minimize(method='SLSQP') e arranque a quente (a solução
para um e é o ponto inicial do seguinte), e guarda o multiplicador u de
f1 <= e (a taxa de troca, u = -df2*/de).

Otimização — deck 6.4 (Método do epsilon-constrangimento).
Reproduz o exemplo dos slides: ver ex06_4_epsilon_constraint.py
O subproblema é resolvido pelo SLSQP do SciPy, como no slide (a UC não
reimplementa o SLSQP); o resto usa só numpy.
Implementação didática — resolve cada subproblema localmente: num problema
não convexo as garantias do slide exigem subproblemas resolvidos globalmente
(p. ex. com vários pontos iniciais).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import NonlinearConstraint, minimize


@dataclass
class EpsResult:
    """Resultado de `epsilon_constraint` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # um ponto por linha (len(E) x n)
    fx: np.ndarray            # [f1(x), f2(x)] de cada ponto
    u: np.ndarray             # multiplicadores de f1 <= e
    nit: int                  # soma das iterações do SLSQP
    nfev: int                 # soma de r.nfev = avaliações de f2 (inclui as das diferenças finitas)
    nfev1: int                # chamadas a f1 (restrição e o seu gradiente numérico)
    ncalls: int               # nfev + nfev1 (chamadas escalares)
    ngev: int = 0             # os gradientes são numéricos, dentro do SLSQP
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 8)))
    cols: tuple = ()
    status: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    flag: int = 0
    message: str = ""


def _multiplicador(r, f1, f2, e, bounds):
    """u de f1 <= e: r.multipliers do SLSQP (SciPy recente, 1.16+); senão,
    r.v de trust-constr a partir do ponto obtido (como no slide, 4.2).
    Estas chamadas extra não entram nas contagens."""
    m = getattr(r, "multipliers", None)
    if m is not None:
        return float(np.ravel(m)[0])
    rt = minimize(f2, r.x, method="trust-constr", bounds=bounds,
                  constraints=NonlinearConstraint(f1, -np.inf, e))
    return float(np.ravel(rt.v[0])[0])


def epsilon_constraint(f1, f2, x0, bounds, E, verbose=False, grafico=False):
    """epsilon-constrangimento: min f2 s.a. f1 <= e, para cada e de E.

    Parâmetros
    ----------
    f1, f2 : funções de x com valor escalar; f1 é o objetivo limitado, f2 o
        objetivo minimizado (como no slide)
    x0 : ponto inicial do primeiro subproblema
    bounds : limites das variáveis, [(l1, u1), (l2, u2), ...]
    E : valores de e (p. ex. np.linspace(0.05, 0.95, 10)); e tem as unidades
        de f1 (não é uma tolerância numérica)
    verbose : se True, imprime a tabela do varrimento
    grafico : se True, desenha os pontos obtidos (precisa de matplotlib)

    Devolve
    -------
    EpsResult com x, fx, u, nit e as contagens: nfev = soma de r.nfev
    (avaliações de f2), nfev1 (chamadas a f1), ncalls = nfev + nfev1, e
    history com uma linha por e: [e, x, f1, f2, u, nf2(e), nf1(e)].
    """
    E = np.ravel(np.asarray(E, float)); ne = len(E)
    x0 = np.asarray(x0, float); n = len(x0)
    X = np.zeros((ne, n)); F = np.zeros((ne, 2)); U = np.zeros(ne)
    hist = np.zeros((ne, n + 6)); status = np.zeros(ne, int)
    cont = {"f1": 0, "f2": 0}               # contadores diretos

    def f1c(z):
        cont["f1"] += 1
        return f1(z)

    def f2c(z):
        cont["f2"] += 1
        return f2(z)

    nfev = 0; nit = 0
    # --- núcleo (o do slide, com as contagens) ------------------------------
    for k, e in enumerate(E):
        antes = (cont["f1"], cont["f2"])
        g = {"type": "ineq", "fun": lambda z, e=e: e - f1c(z)}   # forma da UC: g >= 0
        r = minimize(f2c, x0, bounds=bounds, constraints=[g], method="SLSQP")
        u = _multiplicador(r, f1, f2, e, bounds)                  # taxa de troca
        nfev += r.nfev                                            # avaliações de f2
        nit += r.nit
        x0 = r.x                                                  # a quente
        # ----------------------------------------------------------------------
        nf1k = cont["f1"] - antes[0]; nf2k = cont["f2"] - antes[1]
        X[k] = r.x; F[k] = [f1(r.x), f2(r.x)]; U[k] = u; status[k] = r.status
        hist[k] = [e, *r.x, *F[k], u, nf2k, nf1k]

    cols = ("e",) + tuple("x%d" % (i + 1) for i in range(n)) + ("f1", "f2", "u", "nf2", "nf1")
    falhou = bool(np.any(status != 0))
    msg = ("SLSQP: algum subproblema não terminou bem (ver status)" if falhou
           else "SLSQP: %d subproblemas resolvidos" % ne)

    if verbose:
        print("%6s" % "e" + "".join(" %8s" % c for c in cols[1:n + 1])
              + " %8s %8s %8s %5s %5s" % cols[n + 1:])
        for h in hist:
            print("%6.2f" % h[0] + "".join(" %8.4f" % v for v in h[1:n + 1])
                  + " %8.4f %8.4f %8.4f %5d %5d" % tuple(h[n + 1:]))

    if grafico:
        import matplotlib.pyplot as plt
        plt.plot(F[:, 0], F[:, 1], "o")
        plt.xlabel("$f_1$"); plt.ylabel("$f_2$")
        plt.title("$\\varepsilon$-constrangimento: um ponto por $\\varepsilon$")
        plt.show()

    return EpsResult(x=X, fx=F, u=U, nit=nit, nfev=nfev, nfev1=cont["f1"],
                     ncalls=cont["f1"] + cont["f2"], history=hist, cols=cols,
                     status=status, flag=int(falhou), message=msg)
