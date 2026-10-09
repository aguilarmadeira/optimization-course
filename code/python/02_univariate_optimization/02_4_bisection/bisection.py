"""Método da bisseção sobre f': min f(x) em [a, b], como nas aulas (Deb).

Passo 1: escolher a < b com f'(a) < 0 < f'(b) e eps > 0; fazer x1 = a, x2 = b.
Passo 2: fazer z = (x1 + x2)/2 e calcular f'(z).
Passo 3: se |f'(z)| <= eps, terminar (x* ~ z); senão, se f'(z) < 0, fazer
         x1 = z; se f'(z) > 0, fazer x2 = z; voltar ao passo 2.

Variante (opcional, tolx > 0): parar também quando x2 - x1 <= tolx e devolver
o ponto médio, que garante |x - x*| <= tolx/2. Com tolg = 0 é a bisseção com
um número fixo de reduções, N = ceil(log2((b - a)/tolx)).

Otimização — deck 2.4 (Bisseção).
Reproduz o exemplo dos slides: ver ex02_4_bisection.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. não verifica a unicidade do zero de f').

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

COLS = ("k", "x1", "x2", "z", "f(z)", "df(z)")


@dataclass
class BisectionResult:
    """Resultado de `bisection` (os mesmos campos que `info` no MATLAB)."""
    x: float                  # z da última iteração (ou o ponto médio, na variante)
    fx: float                 # f(x) (nan se f is None) -- só para leitura
    nit: int                  # iterações (avaliações de f'(z))
    ngev: int                 # avaliações de f': nit + 2
    L: float                  # x2 - x1 no fim
    nfev: int = 0             # o método só usa f'
    nhev: int = 0
    history: np.ndarray = field(default_factory=lambda: np.zeros((0, 6)))
    cols: tuple = COLS
    flag: int = 0             # 0: |f'(z)| <= tolg; 1: atingiu kmax; 2: x2 - x1 <= tolx
    message: str = ""


def bisection(f, df, a, b, tolg, kmax=100, tolx=0.0, verbose=False):
    """Bisseção sobre f' em [a, b]; exige df(a) < 0 < df(b).

    Termina quando |f'(z)| <= tolg (flag 0, x = z), como nas aulas e no Deb.
    Com tolx > 0, termina também quando x2 - x1 <= tolx (flag 2, x = ponto
    médio, |x - x*| <= tolx/2). Ao fim de kmax iterações, para com aviso
    (flag 1, x = último z).

    f só é usada para leitura (coluna f(z) da tabela e fx); o método não a
    usa e essas avaliações não entram em nfev (n_f = 0, como nos slides).
    Com f=None nada disso é calculado.

    history tem uma linha por iteração k = 1, ..., nit:
    [k, x1, x2, z, f(z), f'(z)] (x1, x2 no início da iteração).
    """
    if not a < b:
        raise ValueError("bisection: é preciso a < b.")
    temf = f is not None

    # --- verificar o enquadramento: 2 avaliações de f' --------------------
    da = df(a); db = df(b)
    ngev = 2
    if not (da < 0 and db > 0):
        raise ValueError("enquadramento")

    hist = []
    flag = 1
    x1, x2 = a, b
    z = (x1 + x2) / 2
    # --- núcleo (o mesmo dos slides, com as contagens) --------------------
    for k in range(1, kmax + 1):
        z = (x1 + x2) / 2; d = df(z)                 # passo 2
        ngev += 1
        fk = f(z) if temf else np.nan                # só para leitura
        hist.append([k, x1, x2, z, fk, d])
        if abs(d) <= tolg:                           # passo 3: terminou
            flag = 0
            break
        if d < 0:
            x1 = z                                   # mínimo à direita
        else:
            x2 = z                                   # mínimo à esquerda
        if tolx > 0 and x2 - x1 <= tolx:             # variante
            flag = 2
            break
    # ----------------------------------------------------------------------

    x = (x1 + x2) / 2 if flag == 2 else z
    fx = f(x) if temf else np.nan                    # só para leitura
    hist = np.array(hist, float).reshape(-1, 6)
    nit = hist.shape[0]
    L = x2 - x1
    if flag == 0:
        message = "|f'(z)| <= tolg na iteração %d" % nit
    elif flag == 2:
        message = "x2 - x1 = %.2g <= tolx ao fim de %d iterações (ponto médio)" % (L, nit)
    else:
        message = "atingiu o número máximo de iterações (kmax = %d)" % kmax
        warnings.warn("bisection: " + message)

    if verbose:
        print("%3s %9s %9s %9s %10s %10s" % COLS)
        for row in hist:
            print("%3d %9.4f %9.4f %9.4f %10.4f %+10.4f" % tuple(row))

    return BisectionResult(x=x, fx=fx, nit=nit, ngev=ngev, L=L,
                           history=hist, flag=flag, message=message)
