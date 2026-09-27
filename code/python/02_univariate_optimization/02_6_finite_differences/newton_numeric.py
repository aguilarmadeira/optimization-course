"""Newton 1D com f' e f'' por diferenças centradas (só usa valores de f).

Otimização — deck 2.6 (Derivadas numéricas; Newton numérico).
Reproduz o exemplo dos slides: ver ex02_6_finite_differences.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. ddf <= 0, passo amortecido).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

COLS = ("k", "x_k", "h", "df_num", "ddf_num", "x_k+1")


@dataclass
class NewtonNumericResult:
    """Resultado de `newton_numeric` (os mesmos campos que `info` no MATLAB)."""
    x: float                  # último iterando x_{k+1}
    fx: float                 # f(x)
    nit: int                  # passos de Newton
    nfev: int                 # 3 por iteração + 1 (f(x) no fim)
    ngev: int = 0             # derivadas numéricas: contam em nfev
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 6)))
    cols: tuple = COLS
    flag: int = 0             # 0: |df| < tolg; 1: atingiu kmax
    message: str = ""


def newton_numeric(f, x0, tolg, kmax, hrel=0.01, verbose=False):
    """Newton numérico a partir de x0.

    Em cada iteração: h = max(hrel*|x_k|, 1e-4) (com hrel = 0.01 é a regra da
    UC); com f(x_k - h), f(x_k), f(x_k + h) calcula df ~ f'(x_k) e
    ddf ~ f''(x_k) e dá o passo x_{k+1} = x_k - df/ddf. Para quando
    |df| < tolg (df calculado em x_k) e devolve x_{k+1}; aviso se atingir kmax.

    Com h fixo o método estagna perto de x*: converge para o zero de f'_num,
    que não é o zero de f' (erro de truncatura ~ h^2 |f'''|/6).

    history tem uma linha por k = 0, ..., nit-1: [k, x_k, h, df, ddf, x_{k+1}].
    """
    x = x0
    nfev = 0
    hist = []
    flag = 1
    for k in range(1, kmax + 1):
        # --- núcleo (o mesmo dos slides) ----------------------------------
        h = max(hrel * abs(x), 1e-4)
        fm = f(x - h); f0 = f(x); fp = f(x + h)
        df = (fp - fm) / (2 * h)
        ddf = (fp - 2 * f0 + fm) / h**2
        if ddf == 0:
            raise ArithmeticError("f''_num = 0: aumentar h")
        xk = x
        x = x - df / ddf              # passo de Newton
        # ------------------------------------------------------------------
        nfev += 3
        hist.append([k - 1, xk, h, df, ddf, x])
        if abs(df) < tolg:
            flag = 0
            break
    fx = f(x)
    nfev += 1

    hist = np.array(hist, float).reshape(-1, 6)
    nit = hist.shape[0]
    if flag == 0:
        message = "|f'_num(x_k)| < tolg em k = %d; devolve x_%d" % (nit - 1, nit)
    else:
        message = ("atingiu o número máximo de iterações (kmax = %d) sem |f'_num| < tolg" % kmax)
        warnings.warn("newton_numeric: " + message)

    if verbose:
        print("%3s %11s %8s %14s %12s %11s" % COLS)
        for row in hist:
            print("%3d %11.7f %8.4f %+14.6e %12.6f %11.7f" % tuple(row))

    return NewtonNumericResult(x=x, fx=fx, nit=nit, nfev=nfev, history=hist,
                               flag=flag, message=message)
