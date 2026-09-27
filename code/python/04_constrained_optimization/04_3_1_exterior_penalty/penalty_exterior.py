"""Método da função de penalização exterior (EPFM).

min f(x)  s.a.  g_j(x) >= 0 (j = 1..J),  h_l(x) = 0 (l = 1..K):
minimiza-se P(x, R) = f(x) + R [ sum_j <g_j(x)>^2 + sum_l h_l(x)^2 ],
<g> = min(g, 0), para R crescente (R <- c R, c > 1), com arranque a quente.

Otimização — deck 4.3.1 (Método da função de penalização exterior).
Reproduz os exemplos dos slides: ver ex04_3_1_exterior_penalty.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não põe as restrições em escalas comparáveis
nem adapta c ou a tolerância interna ao longo dos ciclos).

Importa `nelder_mead_comp` de code/python/common/. Os exemplos ex*.py tratam do caminho
(import uc_setup); fora deles, é preciso primeiro
    sys.path.insert(0, "<repo>/code/python"); import uc_setup

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np

from nelder_mead_comp import nelder_mead_comp


@dataclass
class PenaltyResult:
    """Resultado de `penalty_exterior` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # x^(t) do último ciclo
    fx: float                 # f(x)
    nit: int                  # ciclos exteriores t feitos
    nfev: int                 # avaliações de f: chamadas a P nos ciclos + 1 (f(x) final)
    R: float                  # R usado no último ciclo
    u: np.ndarray = field(default_factory=lambda: np.zeros(0))       # -2R<g_j(x)>
    lam: np.ndarray = field(default_factory=lambda: np.zeros(0))     # -2R h_l(x)  (info.lambda)
    ngev: int = 0             # o minimizador interno por omissão só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0
    message: str = ""


def _vetor(v):
    return np.atleast_1d(np.asarray(v, float)).ravel()


def penalty_exterior(f, g, h, x0, R0, c, tolviol, tolx, tmax, interno=None, verbose=False):
    """Penalização exterior com arranque a quente (assinatura do slide).

    Parâmetros
    ----------
    f : função objetivo de x (vetor numpy)
    g : x -> valor ou vetor com as desigualdades g_j(x) >= 0 (lambda x: [] se não houver)
    h : x -> valor ou vetor com as igualdades h_l(x) = 0 (lambda x: [] se não houver)
    x0 : ponto inicial x^(0) (não precisa de ser admissível)
    R0 : R^(0) > 0;  c : fator de aumento, c > 1
    tolviol, tolx : para no ciclo t se viol(x^(t)) <= tolviol e
        ||x^(t) - x^(t-1)|| / max(1, ||x^(t-1)||) <= tolx
    tmax : número máximo de ciclos exteriores
    interno : minimizador sem restrições do cap. 3, interno(P, x) -> objeto com
        .x e .nfev (chamadas a P); por omissão, o Nelder–Mead da UC
        (nelder_mead_comp, tolx = 1e-10, tolf = 1e-12, kmax = 2000)
    verbose : se True, imprime a tabela dos ciclos

    Devolve
    -------
    PenaltyResult com x, fx = f(x), nit (ciclos), nfev (todas as chamadas a f:
    chamadas a P nos ciclos + a avaliação final f(x); cada chamada a P avalia
    f uma vez), R (o do último ciclo), u = -2R<g(x)>, lam = -2R h(x), e
    history com uma linha por t = 0, 1, ..., nit (t = 0: ponto inicial):
    [t, R, n_f(ciclo), n_f(acum.), viol, x_1..x_n, u_1..u_J, lambda_1..lambda_K]
    (na linha t, R é o valor usado para obter x^(t); na linha 0, R, u e lambda
    ficam nan).
    """
    if interno is None:
        interno = lambda P, x: nelder_mead_comp(P, x)
    br = lambda z: np.minimum(z, 0)                     # operador <.>
    x = np.array(x0, float).ravel()
    n = x.size
    J = _vetor(g(x)).size
    K = _vetor(h(x)).size
    viol = lambda z: max(np.max(np.abs(br(_vetor(g(z)))), initial=0.0),
                         np.max(np.abs(_vetor(h(z))), initial=0.0))
    nfev = 0
    hist = [[0, np.nan, 0, 0, viol(x), *x] + [np.nan] * (J + K)]
    R = float(R0); Rult = R; flag = 1

    # --- núcleo (o dos slides, com as contagens) --------------------------
    for t in range(1, tmax + 1):
        P = lambda z, R=R: f(z) + R * (np.sum(br(_vetor(g(z)))**2) + np.sum(_vetor(h(z))**2))
        r = interno(P, x)                               # a partir de x
        xn = np.asarray(r.x, float).ravel()
        nfev += r.nfev
        v = viol(xn)
        Rult = R
        hist.append([t, R, r.nfev, nfev, v, *xn,
                     *(-2 * R * br(_vetor(g(xn)))), *(-2 * R * _vetor(h(xn)))])
        if v <= tolviol and np.linalg.norm(xn - x) / max(1, np.linalg.norm(x)) <= tolx:
            x = xn; flag = 0
            break
        x = xn; R = c * R                               # aumentar R
    # ----------------------------------------------------------------------

    fx = f(x); nfev += 1                                # avaliação final
    u = -2 * Rult * br(_vetor(g(x)))                    # u_R = -2R<g>
    lam = -2 * Rult * _vetor(h(x))                      # lambda_R = -2R h
    cols = (("t", "R", "n_f(ciclo)", "n_f(acum)", "viol")
            + tuple("x%d" % (i + 1) for i in range(n))
            + tuple("u%d" % (j + 1) for j in range(J))
            + tuple("lambda%d" % (l + 1) for l in range(K)))
    hist = np.array(hist, float)
    nit = len(hist) - 1
    if flag == 0:
        msg = "viol <= tolviol e passo relativo <= tolx no ciclo %d (R = %g)" % (nit, Rult)
    else:
        msg = "atingiu tmax = %d ciclos (R = %g)" % (tmax, Rult)

    if verbose:
        print(("%3s %9s %10s %9s %10s" + " %10s" * (n + J + K)) % cols)
        for row in hist:
            Rs = "%9s" % "-" if np.isnan(row[1]) else "%9.1e" % row[1]
            rest = " ".join("%10s" % "-" if np.isnan(v) else "%10.6f" % v for v in row[5:])
            print("%3d %s %10d %9d %10.2e %s" % (row[0], Rs, row[2], row[3], row[4], rest))

    return PenaltyResult(x=x, fx=fx, nit=nit, nfev=nfev, R=Rult, u=u, lam=lam,
                         history=hist, cols=cols, flag=flag, message=msg)
