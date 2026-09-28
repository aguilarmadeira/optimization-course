"""Gradiente numérico por diferenças centrais, coordenada a coordenada.

    df/dx_i ~ (f(x + h_i e_i) - f(x - h_i e_i)) / (2 h_i)      erro O(h_i^2)

Custo: 2n avaliações de f (conta em n_f; n_g = 0).

Otimização — decks 3.3 (Derivadas numéricas em n dimensões) e 3.3.1.
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import numpy as np


def grad_fd(f, x, h=None):
    """Gradiente de f em x por diferenças centrais.

    Parâmetros
    ----------
    f : função de R^n em R
    x : ponto (vetor com n componentes)
    h : (opcional) passo, escalar ou vetor com n componentes; por omissão a
        regra de 2.6 em cada coordenada: h_i = max(0.01 |x_i|, 1e-4)

    Devolve
    -------
    (g, nfev): o gradiente aproximado e o número de avaliações de f (2n).
    Mesmo compromisso de 2.6: h grande, erro de truncatura; h pequeno,
    erro de arredondamento.
    """
    x = np.asarray(x, float).ravel()
    n = x.size
    if h is None:
        h = np.maximum(0.01 * np.abs(x), 1e-4)
    h = np.broadcast_to(np.asarray(h, float), (n,))
    g = np.zeros(n)
    nfev = 0
    # --- núcleo (o dos slides) ---------------------------------------------
    for i in range(n):
        e = np.zeros(n); e[i] = h[i]
        g[i] = (f(x + e) - f(x - e)) / (2 * h[i])
        nfev += 2
    # ----------------------------------------------------------------------
    return g, nfev
