"""Derivadas numéricas f'(x) e f''(x) por diferenças centradas.

Com os mesmos três valores f(x-h), f(x), f(x+h) (3 avaliações de f):
    df  = (f(x+h) - f(x-h)) / (2h)            erro O(h^2)
    ddf = (f(x+h) - 2 f(x) + f(x-h)) / h^2    erro O(h^2)

h grande: erro de truncatura; h pequeno: erro de arredondamento
(cancelamento). A regra de 1 % é uma regra de trabalho, não um passo ótimo
universal.

Otimização — deck 2.6 (Derivadas numéricas).
Reproduz o exemplo dos slides: ver ex02_6_finite_differences.py
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""


def num_derivs(f, x, h=None):
    """Devolve (df, ddf, info) com as aproximações centradas de f'(x), f''(x).

    h por omissão: a regra da UC, h = 0.01|x| se |x| > 0.01, 1e-4 caso
    contrário (= max(0.01|x|, 1e-4)).

    info (dicionário, os mesmos campos que no MATLAB):
      nfev = 3, ngev = nhev = 0 (as derivadas numéricas contam em nfev),
      h, fm, f0, fp (= f(x-h), f(x), f(x+h)),
      df_fwd = (f(x+h) - f(x))/h  (progressiva, erro O(h)),
      df_bwd = (f(x) - f(x-h))/h  (regressiva, erro O(h)),
      estas duas sem custo extra (reutilizam os mesmos valores).
    """
    if h is None:
        h = max(0.01 * abs(x), 1e-4)

    # --- núcleo (o mesmo dos slides) --------------------------------------
    fm = f(x - h); f0 = f(x); fp = f(x + h)
    df = (fp - fm) / (2 * h)
    ddf = (fp - 2 * f0 + fm) / h**2
    # ----------------------------------------------------------------------

    info = dict(nfev=3, ngev=0, nhev=0, h=h, fm=fm, f0=f0, fp=fp,
                df_fwd=(fp - f0) / h, df_bwd=(f0 - fm) / h)
    return df, ddf, info
