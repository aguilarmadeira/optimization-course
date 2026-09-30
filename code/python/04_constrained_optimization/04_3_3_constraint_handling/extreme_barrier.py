"""Barreira extrema: +inf fora da região admissível, sem avaliar f.

    extreme_barrier(f, g, x)     restrições explícitas g(x) >= 0 (escalar ou
                                 vetor; convenção do cap. 4): f(x) se todas
                                 as g_j(x) >= 0; senão
                                 +inf e f NÃO é avaliada.
    extreme_barrier(F, None, x)  restrição oculta (não há g): F(x), ou +inf
                                 se a avaliação falha (NaN, inf ou exceção).

Uso: F = lambda x: extreme_barrier(f, g, x); depois, qualquer método que só
compare valores (p. ex. hooke_jeeves, 3.2.5) minimiza F. É o que fazem o DMS
e o MultiGLODS: f = +inf nos pontos inadmissíveis, que não são avaliados.

Otimização — deck 4.3.3 (restrições na prática).
J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import numpy as np


def extreme_barrier(f, g, x):
    if g is not None:
        return np.inf if np.any(np.asarray(g(x)) < 0) else f(x)
    try:
        y = f(x)
    except Exception:              # o simulador falhou com erro
        return np.inf
    return np.inf if np.isnan(y) else y   # ... ou devolveu NaN
