"""Método da barreira logarítmica (penalização interior).

min f(x)  s.a.  g_j(x) >= 0 (j = 1..J), a partir de x^(0) com todos os g_j > 0:
minimiza-se P(x, R) = f(x) - R sum_j ln g_j(x) (P = inf se algum g_j <= 0)
para R decrescente (R <- c R, 0 < c < 1), com arranque a quente.

Otimização — deck 4.3.2 (Método da função de penalização interior).
Reproduz os exemplos dos slides: ver ex04_3_2_barrier.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não tem Fase I: exige um ponto estritamente
admissível; não trata igualdades).

Importa `nelder_mead_comp` de code/python/common/. Os exemplos ex*.py tratam do caminho
(import uc_setup); fora deles, é preciso primeiro
    sys.path.insert(0, "<repo>/code/python"); import uc_setup

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np

from nelder_mead_comp import nelder_mead_comp


def _vetor(v):
    return np.atleast_1d(np.asarray(v, float)).ravel()


def pbar(f, g, x, R, contador=None):
    """P(x, R) = f(x) - R sum_j ln g_j(x); inf (sem avaliar f) se algum g_j(x) <= 0.

    contador : (opcional) dicionário; cada chamada com x fora de X soma 1 a
    contador["fora"] (essas chamadas avaliam g mas não f).
    """
    gx = _vetor(g(x))
    if np.any(gx <= 0):                    # nunca avaliar ln de um número <= 0
        if contador is not None:
            contador["fora"] += 1
        return np.inf
    return f(x) - R * np.sum(np.log(gx))


@dataclass
class BarrierResult:
    """Resultado de `barrier_log` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # x^(t) do último ciclo (estritamente admissível)
    fx: float                 # f(x)
    nit: int                  # ciclos exteriores t feitos
    nfev: int                 # chamadas a P nos ciclos + 1 (f(x) final): majorante de n_f
    nPev: int                 # chamadas a P nos ciclos (sem a avaliação final)
    nfora: int                # chamadas a P com x fora de X (avaliam g, não f)
    nfev_efetivo: int         # avaliações de f de facto: nPev - nfora + 1
    R: float                  # R usado no último ciclo
    u: np.ndarray = field(default_factory=lambda: np.zeros(0))   # u_j = R/g_j(x)
    ngev: int = 0             # o minimizador interno por omissão só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    flag: int = 0
    message: str = ""


def barrier_log(f, g, x0, R0, c, tolcomp, tolx, tmax, interno=None, verbose=False):
    """Barreira logarítmica com arranque a quente (assinatura do slide).

    Parâmetros
    ----------
    f : função objetivo de x (vetor numpy)
    g : x -> valor ou vetor com as desigualdades g_j(x) >= 0
    x0 : ponto inicial ESTRITAMENTE admissível (todos os g_j(x0) > 0)
    R0 : R^(0) > 0;  c : fator de redução, 0 < c < 1
    tolcomp, tolx : para no ciclo t se J R <= tolcomp (J R = sum_j u_j g_j,
        complementaridade perturbada) e ||x^(t) - x^(t-1)|| / max(1, ||x^(t-1)||) <= tolx
    tmax : número máximo de ciclos exteriores
    interno : minimizador sem restrições do cap. 3, interno(P, x) -> objeto com
        .x e .nfev (chamadas a P); tem de aceitar P = inf. Por omissão, o
        Nelder–Mead da UC (nelder_mead_comp, tolx = 1e-10, tolf = 1e-12, kmax = 2000)
    verbose : se True, imprime a tabela dos ciclos

    Devolve
    -------
    BarrierResult com x, fx = f(x), nit (ciclos), nfev = chamadas a P + 1
    (como no slide: majorante de n_f), nPev (chamadas a P), nfora (chamadas a
    P fora de X, que avaliam g mas não f), nfev_efetivo = nPev - nfora + 1,
    R (o do último ciclo), u = R/g(x), e history com uma linha por
    t = 0, 1, ..., nit (t = 0: ponto inicial):
    [t, R, chamadas a P (ciclo), fora de X (ciclo), chamadas a P (acum.), x_1..x_n, u_1..u_J]
    (na linha t, R é o valor usado para obter x^(t); na linha 0, R e u ficam nan).
    """
    if interno is None:
        interno = lambda P, x: nelder_mead_comp(P, x)
    x = np.array(x0, float).ravel()
    n = x.size
    g0 = _vetor(g(x)); J = g0.size
    if np.any(g0 <= 0):
        raise ValueError("barrier_log: x0 tem de ser estritamente admissível (g_j(x0) > 0).")
    cont = {"fora": 0}
    nPev = 0
    hist = [[0, np.nan, 0, 0, 0, *x] + [np.nan] * J]
    R = float(R0); Rult = R; flag = 1

    # --- núcleo (o dos slides, com as contagens) --------------------------
    for t in range(1, tmax + 1):
        P = lambda z, R=R: pbar(f, g, z, R, cont)      # inf fora de X
        fora0 = cont["fora"]
        r = interno(P, x)                              # a partir de x
        xn = np.asarray(r.x, float).ravel()
        nPev += r.nfev
        u = R / _vetor(g(xn))                          # u_j = R/g_j
        Rult = R
        hist.append([t, R, r.nfev, cont["fora"] - fora0, nPev, *xn, *u])
        if u.size * R <= tolcomp and np.linalg.norm(xn - x) / max(1, np.linalg.norm(x)) <= tolx:
            x = xn; flag = 0                           # J*R = sum(u.*g)
            break
        x = xn; R = c * R                              # reduzir R
    # ----------------------------------------------------------------------

    fx = f(x)                                          # avaliação final
    nfev = nPev + 1
    nfora = cont["fora"]
    u = Rult / _vetor(g(x))
    cols = (("t", "R", "P(ciclo)", "fora(ciclo)", "P(acum)")
            + tuple("x%d" % (i + 1) for i in range(n))
            + tuple("u%d" % (j + 1) for j in range(J)))
    hist = np.array(hist, float)
    nit = len(hist) - 1
    if flag == 0:
        msg = "J*R <= tolcomp e passo relativo <= tolx no ciclo %d (R = %g)" % (nit, Rult)
    else:
        msg = "atingiu tmax = %d ciclos (R = %g)" % (tmax, Rult)

    if verbose:
        print(("%3s %9s %9s %11s %8s" + " %10s" * (n + J)) % cols)
        for row in hist:
            Rs = "%9s" % "-" if np.isnan(row[1]) else "%9.1e" % row[1]
            rest = " ".join("%10s" % "-" if np.isnan(v) else "%10.6f" % v for v in row[5:])
            print("%3d %s %9d %11d %8d %s" % (row[0], Rs, row[2], row[3], row[4], rest))

    return BarrierResult(x=x, fx=fx, nit=nit, nfev=nfev, nPev=nPev, nfora=nfora,
                         nfev_efetivo=nPev - nfora + 1, R=Rult, u=u, history=hist,
                         cols=cols, flag=flag, message=msg)
