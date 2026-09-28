"""Função de Rosenbrock, com gradiente e Hessiana só quando são pedidos.

    f(x) = (1 - x1)^2 + 100 (x2 - x1^2)^2, mínimo f = 0 em (1, 1).

É o equivalente Python do fg.m (que usa nargout): aqui diz-se quantas saídas
se querem com nout (1: f; 2: f, g; 3: f, g, H).  Quem só precisa de f (p. ex.
a pesquisa em linha) não paga o resto.

Otimização — deck 3.3 (Métodos com derivadas) e comparação do cap. 3.
J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import numpy as np


def fg(x, nout=1):
    """f(x); com nout = 2 devolve (f, g); com nout = 3 devolve (f, g, H)."""
    f = (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2
    if nout == 1:
        return f
    g = np.array([-2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0]**2),     # gradiente só se for pedido
                  200 * (x[1] - x[0]**2)])
    if nout == 2:
        return f, g
    H = np.array([[2 - 400 * (x[1] - 3 * x[0]**2), -400 * x[0]],       # Hessiana só se for pedida
                  [-400 * x[0], 200]])
    return f, g, H
