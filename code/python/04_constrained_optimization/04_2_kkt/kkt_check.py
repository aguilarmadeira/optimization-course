"""Verifica as condições KKT num ponto, na convenção da UC.

Otimização — deck 4.2 (Condições KKT). Usada em ex04_2_kkt.py.
Implementação didática: só verifica as condições de 1.ª ordem; não verifica
a qualificação das restrições (LICQ) nem classifica o ponto.

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
from types import SimpleNamespace

import numpy as np


def kkt_check(df, g=None, Jg=None, u=None, h=None, Jh=None, lam=None, tol=1e-8):
    """Condições KKT de min f s.a. g_j(x) >= 0, h_l(x) = 0, com L = f - u'g - lambda'h.

    Tudo avaliado no ponto x:
    df : gradiente de f (n,)
    g, Jg, u : valores g_j(x) (J,), jacobiana (J, n) com linha j = grad g_j,
               multiplicadores (J,);  None se não houver desigualdades
    h, Jh, lam : o mesmo para as igualdades (K,), (K, n), (K,)
    tol : tolerância

    Devolve (ok, res): ok é True se as quatro condições se cumprem (<= tol);
    res tem os resíduos
      estac  estacionariedade  ||df - Jg'u - Jh'lambda||_inf
      admis  admissibilidade   max(max(-g_j, 0), |h_l|)
      compl  complementaridade max |u_j g_j|
      sinal  sinal             max(-u_j, 0)   (u_j >= 0; lambda livre)
      ativas índices j (a começar em 1, como nos slides) com |g_j| <= tol
      ok     o mesmo que ok
    """
    df = np.ravel(np.asarray(df, float)); n = df.size
    if g is None:
        g, Jg, u = np.zeros(0), np.zeros((0, n)), np.zeros(0)
    if h is None:
        h, Jh, lam = np.zeros(0), np.zeros((0, n)), np.zeros(0)
    g = np.ravel(np.asarray(g, float)); u = np.ravel(np.asarray(u, float))
    h = np.ravel(np.asarray(h, float)); lam = np.ravel(np.asarray(lam, float))
    Jg = np.asarray(Jg, float).reshape(g.size, n)
    Jh = np.asarray(Jh, float).reshape(h.size, n)

    r = df - Jg.T @ u - Jh.T @ lam
    res = SimpleNamespace(
        estac=float(np.linalg.norm(r, np.inf)),
        admis=float(np.max(np.concatenate([[0.0], -g, np.abs(h)]))),
        compl=float(np.max(np.concatenate([[0.0], np.abs(u * g)]))),
        sinal=float(np.max(np.concatenate([[0.0], -u]))),
        ativas=[j + 1 for j in range(g.size) if abs(g[j]) <= tol])
    ok = max(res.estac, res.admis, res.compl, res.sinal) <= tol
    res.ok = ok
    return ok, res
