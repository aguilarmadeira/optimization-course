"""Interpolação quadrática sucessiva (método de Powell, 1D).

Otimização — deck 2.3 (Interpolação quadrática sucessiva).
Reproduz o exemplo dos slides: ver ex02_3_powell.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. limitar o passo, pontos repetidos).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

COLS = ("k", "x1", "x2", "x3", "f1", "f2", "f3", "a1", "a2", "xb", "f(xb)")


@dataclass
class PowellResult:
    """Resultado de `powell_quadratic` (os mesmos campos que `info` no MATLAB)."""
    x: float                  # o melhor dos quatro pontos da última iteração
    fx: float
    nit: int                  # iterações (cada uma calcula um xb)
    nfev: int                 # 3 iniciais + 1 por iteração + 1 por reflexão
    nrefl: int                # reflexões (salvaguarda para a2 <= 0)
    ngev: int = 0
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 11)))
    cols: tuple = COLS
    flag: int = 0             # 0: critério satisfeito; 1: atingiu kmax
    message: str = ""


def powell_quadratic(f, x1, Delta, tolx, tolf, kmax=100, verbose=False):
    """Interpolação quadrática sucessiva a partir de x1, com passo Delta.

    Para quando |F_min - f(xb)| <= tolf e |X_min - xb| <= tolx (as duas).
    Salvaguardas: se a2 <= 0 (parábola sem mínimo), reflete-se o pior ponto
    no melhor, 2*X_min - x_pior, e repete-se; depois de cada iteração
    guardam-se o melhor ponto e os dois que o enquadram (se não existirem,
    os três melhores).

    Devolve um PowellResult; history tem uma linha por iteração k = 0, 1, ...:
    [k, x1, x2, x3, f1, f2, f3, a1, a2, xb, f(xb)].
    """
    # --- arranque: 3 avaliações ------------------------------------------
    x2 = x1 + Delta
    f1 = f(x1); f2 = f(x2)
    if f1 > f2:
        x3 = x1 + 2 * Delta           # f desce: o terceiro ponto vai para a frente
    else:
        x3 = x1 - Delta               # f sobe: o terceiro ponto vai para trás
    f3 = f(x3)
    nfev = 3
    X = np.array([x1, x2, x3], float)
    F = np.array([f1, f2, f3], float)

    hist = []
    k = 0; nrefl = 0; flag = 1
    xb = fb = None
    while k < kmax:
        # --- núcleo de uma iteração (o mesmo dos slides) -------------------
        i = np.argsort(X, kind="stable"); X = X[i]; F = F[i]   # x1 < x2 < x3
        imin = int(np.argmin(F)); Fmin = F[imin]; Xmin = X[imin]
        a1 = (F[1] - F[0]) / (X[1] - X[0])
        a2 = ((F[2] - F[0]) / (X[2] - X[0]) - a1) / (X[2] - X[1])
        if a2 <= 0:                   # parábola sem mínimo:
            # refletir o pior ponto no melhor e repetir
            iw = int(np.argmax(F))
            X[iw] = 2 * Xmin - X[iw]; F[iw] = f(X[iw])
            nfev += 1; nrefl += 1
            if nrefl > kmax:
                flag = 1
                break
            continue
        xb = (X[0] + X[1]) / 2 - a1 / (2 * a2)
        fb = f(xb)                    # 1 avaliação
        nfev += 1
        # ------------------------------------------------------------------
        hist.append([k, *X, *F, a1, a2, xb, fb])
        k += 1

        # parar se xb e fb mudam pouco
        if abs(Fmin - fb) <= tolf and abs(Xmin - xb) <= tolx:
            flag = 0
            break

        # senão, guardar o melhor ponto e os dois que o enquadram
        Xa = np.append(X, xb); Fa = np.append(F, fb)
        i = np.argsort(Xa, kind="stable"); Xa = Xa[i]; Fa = Fa[i]
        ib = int(np.argmin(Fa))
        if 0 < ib < 3:
            X = Xa[ib - 1:ib + 2].copy(); F = Fa[ib - 1:ib + 2].copy()
        else:                         # sem enquadramento: os três melhores
            j = np.argsort(Fa, kind="stable")[:3]
            X = Xa[j]; F = Fa[j]

    # o melhor dos pontos disponíveis (na paragem: os quatro da última iteração)
    if k > 0:
        Xa = np.append(X, xb); Fa = np.append(F, fb)
    else:
        Xa = X; Fa = F
    ib = int(np.argmin(Fa)); x = float(Xa[ib]); fx = float(Fa[ib])

    hist = np.array(hist, float).reshape(-1, 11)
    if flag == 0:
        message = ("|F_min - f(xb)| <= tolf e |X_min - xb| <= tolx em %d iterações" % k)
    else:
        message = "atingiu o número máximo de iterações (kmax = %d)" % kmax
        warnings.warn("powell_quadratic: " + message)

    if verbose:
        print("%3s %9s %9s %9s %10s %10s %10s %9s %8s %9s %10s" % COLS)
        for row in hist:
            print("%3d %9.4f %9.4f %9.4f %10.4f %10.4f %10.4f %9.3f %8.3f %9.4f %10.4f"
                  % tuple(row))

    return PowellResult(x=x, fx=fx, nit=k, nfev=nfev, nrefl=nrefl,
                        history=hist, flag=flag, message=message)
