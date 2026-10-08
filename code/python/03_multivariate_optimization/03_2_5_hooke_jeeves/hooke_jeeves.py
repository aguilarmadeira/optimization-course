"""Método de Hooke–Jeeves (pesquisa em padrão), sem derivadas.

Alterna o movimento exploratório (exploratory.py, eixo a eixo: avalia +-Delta_j
e fica com o melhor dos três pontos, como em Deb, 2012) com o movimento de
padrão X^P_{k+1} = X_k + (X_k - X_{k-1}) = 2 X_k - X_{k-1} (a = 2), seguido de
nova exploração em torno de X^P_{k+1}.

Duas variantes (parâmetro variante):
  * "deb" (por omissão; a das aulas, Deb, 2012, sec. 3.3.3):
      passo 2  explorar em torno de X_k; se melhorar, X_{k+1} = resultado, passo 4;
               senão, passo 3;
      passo 3  se ||Delta|| < eps, parar; senão Delta = Delta/R e voltar ao passo 2;
      passo 4  padrão X^P_{k+1} = 2 X_k - X_{k-1};
      passo 5  explorar em torno de X^P_{k+1}; o resultado é X_{k+1};
      passo 6  se f(X_{k+1}) < f(X_k), voltar ao passo 4; senão, passo 3
               (com X_k como ponto base).
    Aqui T é eps (um escalar) e R o fator de redução (2 por omissão).
  * "1961" (Hooke e Jeeves, 1961): quando o padrão falha, recua-se para X_k e
    explora-se de novo com o mesmo Delta; Delta só se divide por 2 quando a
    exploração em torno da base falha; pára quando todos os Delta_j < T_j.
Em ambas, Delta só diminui (não é reposto em Delta_0 a cada sucesso).

Otimização — deck 3.2.5 (Método de Hooke–Jeeves).
Reproduz os exemplos dos slides: ver ex03_2_5_hooke_jeeves.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não escala as variáveis: P0 e T devem ser
escolhidos relativamente à escala de cada variável). Parar a passo
finito não certifica estacionariedade.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from dataclasses import dataclass, field

import numpy as np

from exploratory import exploratory


@dataclass
class HookeJeevesResult:
    """Resultado de `hooke_jeeves` (os mesmos campos que `info` no MATLAB)."""
    x: np.ndarray             # ponto base final
    fx: float                 # f(x) (já avaliado; não há avaliação extra no fim)
    nit: int                  # número de movimentos (linhas de history)
    nfev: int                 # 1 + explorações (2n cada) + 1 por ponto tentativo
    ngev: int = 0             # o método só usa valores de f
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    cols: tuple = ()
    ftrace: np.ndarray = field(default_factory=lambda: np.zeros(0))
    P: np.ndarray = field(default_factory=lambda: np.zeros(0))
    flag: int = 0
    message: str = ""


def hooke_jeeves(f, x0, a, P0, T, kmax=10000, verbose=False, menor=None, variante="deb", R=2.0):
    """Hooke–Jeeves a partir de x0, com fator a, passos iniciais P0 (Delta_0) e tolerância T.

    Parâmetros
    ----------
    f : função de um vetor x (array de dimensão n)
    x0 : ponto inicial
    a : fator do movimento de padrão (tipicamente 2); None = 2
    P0 : passos iniciais Delta_0, um por coordenada (escalar = igual em todas)
    T : variante "deb": eps (escalar), pára quando ||Delta|| < eps;
        variante "1961": tolerâncias, uma por coordenada (escalar = igual em
        todas), pára quando todos os Delta_j < T_j
    kmax : número máximo de movimentos (explorações a partir da base +
        movimentos de padrão)
    verbose : se True, imprime uma linha por movimento
    menor : (opcional) função menor(u, v) que diz se o valor u é melhor do que
        v; por omissão, u < v (ver exploratory; deck 4.3.3, regra de
        admissibilidade, com f a devolver o par (v(x), f(x)))
    variante : "deb" (por omissão, a das aulas) ou "1961" (ver o cabeçalho)
    R : fator de redução de Delta na variante "deb" (2)

    Devolve
    -------
    HookeJeevesResult com x, fx, nit, nfev = 1 (f(x0)) + 2n por exploração
    + 1 por ponto tentativo (cada padrão gasta 1 + 2n);
    history com uma linha por movimento:
    [m, tipo, xb(1:n), xt(1:n), f(xt), xe(1:n), f(xe), f_ref, P(1:n), sucesso, nfev]
    (tipo = 0: exploração em torno de xb, xt e f(xt) = NaN, f_ref = f(xb);
    tipo = 1: padrão a partir da base xb, com ponto tentativo xt e exploração
    em torno de xt até xe', nas colunas xe, f(xe), f_ref = f(xe) anterior;
    sucesso = 1 se f(xe) < f_ref; nfev acumulado); ftrace com os valores de f
    pela ordem em que foram avaliados; flag = 0 se parou pelo critério de
    paragem, 1 se atingiu kmax.
    """
    if a is None:
        a = 2.0
    if variante not in ("deb", "1961"):
        raise ValueError('hooke_jeeves: variante tem de ser "deb" ou "1961".')
    deb = variante == "deb"
    if menor is None:
        menor = lambda u, v: u < v        # a ordem habitual
    reg = lambda y: y[-1] if np.ndim(y) else y   # valor registado
    xb = np.array(x0, dtype=float).ravel()
    n = xb.size
    P = np.array(P0, dtype=float).ravel() * np.ones(n)
    if deb:
        if np.size(T) != 1:
            raise ValueError('hooke_jeeves: na variante "deb", T é eps (um escalar).')
        eps = float(np.ravel(T)[0])
        T = np.full(n, eps)
    else:
        T = np.array(T, dtype=float).ravel() * np.ones(n)
    if np.any(P <= 0) or np.any(T <= 0):
        raise ValueError("hooke_jeeves: P0 e T têm de ser positivos.")

    fb = f(xb); nfev = 1                  # f(x0) conta
    ftrace = [reg(fb)]
    hist = []
    m = 0
    flag = 1
    xe, fe = xb.copy(), fb
    nan = [np.nan] * n

    while m < kmax:                       # passo 2: exploração em torno da base xb (= X_k)
        r = exploratory(f, xb, fb, P, menor)
        xe, fe = r.x, r.fx
        nfev += r.nfev; ftrace.extend(r.ftrace)
        m += 1
        hist.append([m, 0, *xb, *nan, np.nan, *xe, reg(fe), reg(fb), *P, menor(fe, fb), nfev])
        if not menor(fe, fb):             # falhou (fe >= fb): passo 3
            if deb:
                if np.linalg.norm(P) < eps:
                    flag = 0
                    break
                P = P / R
            else:
                P = P / 2
                if np.all(P < T):
                    flag = 0
                    break
        else:
            while m < kmax:               # movimentos de padrão
                xt = xb + a * (xe - xb); ft = f(xt)          # padrão: ponto tentativo
                nfev += 1; ftrace.append(reg(ft))
                r = exploratory(f, xt, ft, P, menor)         # exploração em torno de xt
                xn, fn = r.x, r.fx
                nfev += r.nfev; ftrace.extend(r.ftrace)
                m += 1
                hist.append([m, 1, *xb, *xt, reg(ft), *xn, reg(fn), reg(fe), *P, menor(fn, fe), nfev])
                if not menor(fn, fe):     # rejeitado (fn >= fe): recuar para xe (= X_k)
                    xb, fb = xe, fe
                    if deb:               # passo 3 logo a seguir
                        if np.linalg.norm(P) < eps:
                            flag = 0
                        else:
                            P = P / R
                    break
                else:                     # aceite: continuar o padrão
                    xb, fb = xe, fe
                    xe, fe = xn, fn
            if flag == 0:
                break
    if menor(fe, fb):                     # saiu por kmax a meio de um padrão aceite
        xb, fb = xe, fe
    x, fx = xb, fb

    cols = (("m", "tipo") + tuple("xb%d" % (j + 1) for j in range(n))
            + tuple("xt%d" % (j + 1) for j in range(n)) + ("f(xt)",)
            + tuple("xe%d" % (j + 1) for j in range(n)) + ("f(xe)", "f_ref")
            + tuple("P%d" % (j + 1) for j in range(n)) + ("sucesso", "nfev"))
    hist = np.array(hist, dtype=float).reshape(-1, 4 * n + 7)

    if flag == 0:
        if deb:
            message = "||Delta|| = %.3g < eps = %.3g ao fim de %d movimentos" % (np.linalg.norm(P), eps, m)
        else:
            message = "todos os P_j < T_j (max P_j = %.3g) ao fim de %d movimentos" % (P.max(), m)
    else:
        message = "atingiu kmax = %d movimentos (max P_j = %.3g)" % (kmax, P.max())

    if verbose:
        print("%4s %-9s" % ("m", "movimento")
              + "".join("%9s" % c for c in cols[2:2 + n])
              + "".join("%9s" % c for c in cols[2 + n:2 + 2 * n]) + "%10s" % "f(xt)"
              + "".join("%9s" % c for c in cols[3 + 2 * n:3 + 3 * n]) + "%10s" % "f(xe)"
              + "".join("%9s" % c for c in cols[5 + 3 * n:5 + 4 * n])
              + "  %-10s %5s" % ("resultado", "n_f"))
        for h in hist:
            tipo = int(h[1])
            s = "%4d " % h[0] + ("explor.  " if tipo == 0 else "padrão   ")
            s += "".join("%9.4f" % v for v in h[2:2 + n])
            if tipo == 0:
                s += "%9s" % "-" * n + "%10s" % "-"
            else:
                s += "".join("%9.4f" % v for v in h[2 + n:2 + 2 * n]) + "%10.4f" % h[2 + 2 * n]
            s += "".join("%9.4f" % v for v in h[3 + 2 * n:3 + 3 * n]) + "%10.4f" % h[3 + 3 * n]
            s += "".join("%9.4f" % v for v in h[5 + 3 * n:5 + 4 * n])
            if tipo == 0:
                res = "melhora" if h[5 + 4 * n] else ("falha" if deb else "falha:P/2")
            else:
                res = "aceite" if h[5 + 4 * n] else "rejeitado"
            print(s + "  %-10s %5d" % (res, h[6 + 4 * n]))

    return HookeJeevesResult(x=x, fx=fx, nit=m, nfev=nfev, history=hist, cols=cols,
                             ftrace=np.array(ftrace), P=P, flag=flag, message=message)
