"""Pesquisa em linha para os métodos de 3.3: alpha_k ~ argmin_{alpha>0} phi(alpha).

Dois métodos, ambos só com valores de phi (cap. 2):

* "brent" (por omissão) — a pesquisa de elevada precisão usada nas figuras e
  na tabela comparativa do capítulo 3: enquadramento a partir de (0, 1e-3)
  (passos em razão áurea + interpolação parabólica) seguido do método de
  Brent (secção áurea + interpolação parabólica, 2.2-2.3) com tolerância
  relativa tol = 1e-10.  É uma tradução linha a linha de
  scipy.optimize.bracket + Brent (minimize_scalar, method="brent"), por isso
  dá exatamente os mesmos pontos e as mesmas contagens que os scripts das
  figuras (cerca de 24 avaliações por pesquisa na quadrática do 3.3.1).
* "fminbnd" — o algoritmo do fminbnd do MATLAB (Brent em [0, amax], secção
  áurea + interpolação parabólica, TolX = 1e-4 por omissão), com o
  alargamento do slide «No computador»: se alpha ~ amax, repetir com
  10*amax.  Na quadrática do 3.3.1 gasta 6 avaliações por pesquisa.

Otimização — decks 3.3.1, 3.3.2 e 3.3.3 (pesquisa em linha; cap. 2).
Utilitário comum (code/python/common/): usado por steepest_descent (3.3.1),
newton_nd (3.3.2), conj_grad (3.3.3) e pela comparação do cap. 3; existe só
aqui. Os exemplos encontram-no através de `import uc_setup`.
Implementação didática — algumas salvaguardas de software profissional
não estão incluídas (p. ex. condições de Wolfe).

J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).
"""
import math


def line_search(phi, method="brent", tol=None, amax=1.0):
    """Minimiza phi(alpha) para alpha >= 0.

    Parâmetros
    ----------
    phi : função de uma variável, phi(alpha) = f(x_k + alpha d_k)
    method : "brent" (por omissão) ou "fminbnd"
    tol : "brent": tolerância relativa de Brent (por omissão 1e-10);
          "fminbnd": TolX (por omissão 1e-4)
    amax : "fminbnd": extremo direito inicial do intervalo [0, amax]

    Devolve
    -------
    (alpha, nfev, fvals): o passo alpha >= 0, o número de avaliações de phi
    e a lista dos valores de phi pela ordem em que foram calculados.
    Nota: "brent" volta a avaliar phi(0) = f(x_k) (é a primeira avaliação
    do enquadramento) e essa avaliação conta.
    """
    fvals = []

    def F(t):
        v = phi(t)
        fvals.append(v)
        return v

    if method == "brent":
        alpha = _brent(F, 1e-10 if tol is None else tol)
        alpha = max(alpha, 0.0)
    elif method == "fminbnd":
        tolx = 1e-4 if tol is None else tol
        b = amax
        alpha = _fminbnd(F, 0.0, b, tolx)
        while alpha > 0.99 * b:            # mínimo no extremo: alargar o intervalo
            b = 10 * b
            alpha = _fminbnd(F, 0.0, b, tolx)
    else:
        raise ValueError("line_search: método desconhecido '%s'" % method)
    return alpha, len(fvals), fvals


def _bracket(F, xa=0.0, xb=1e-3, grow_limit=110.0, maxiter=1000):
    """Enquadramento (xa, xb, xc) com F(xb) <= F(xa), F(xc) (tradução de scipy.optimize.bracket)."""
    gold = 1.618034
    verysmall = 1e-21
    fa = F(xa)
    fb = F(xb)
    if fa < fb:                            # trocar para que fa > fb
        xa, xb = xb, xa
        fa, fb = fb, fa
    xc = xb + gold * (xb - xa)
    fc = F(xc)
    it = 0
    while fc < fb:
        tmp1 = (xb - xa) * (fb - fc)
        tmp2 = (xb - xc) * (fb - fa)
        val = tmp2 - tmp1
        denom = 2.0 * verysmall if abs(val) < verysmall else 2.0 * val
        w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom   # parábola
        wlim = xb + grow_limit * (xc - xb)
        if it > maxiter:
            raise RuntimeError("line_search: não encontrou enquadramento")
        it += 1
        if (w - xc) * (xb - w) > 0.0:
            fw = F(w)
            if fw < fc:
                xa = xb; xb = w; fa = fb; fb = fw
                break
            elif fw > fb:
                xc = w; fc = fw
                break
            w = xc + gold * (xc - xb)
            fw = F(w)
        elif (w - wlim) * (wlim - xc) >= 0.0:
            w = wlim
            fw = F(w)
        elif (w - wlim) * (xc - w) > 0.0:
            fw = F(w)
            if fw < fc:
                xb = xc; xc = w
                w = xc + gold * (xc - xb)
                fb = fc; fc = fw
                fw = F(w)
        else:
            w = xc + gold * (xc - xb)
            fw = F(w)
        xa = xb; xb = xc; xc = w
        fa = fb; fb = fc; fc = fw
    ok = ((fb < fc and fb <= fa) or (fb < fa and fb <= fc)) and \
         (xa < xb < xc or xc < xb < xa) and \
         all(math.isfinite(v) for v in (xa, xb, xc))
    if not ok:
        raise RuntimeError("line_search: não encontrou enquadramento válido")
    return xa, xb, xc, fa, fb, fc


def _brent(F, tol, maxiter=500):
    """Método de Brent a partir do enquadramento (tradução de scipy.optimize.Brent)."""
    xa, xb, xc, fa, fb, fc = _bracket(F)
    mintol = 1.0e-11
    cg = 0.3819660
    x = w = v = xb
    fw = fv = fx = fb
    a, b = (xa, xc) if xa < xc else (xc, xa)
    deltax = 0.0
    rat = 0.0
    it = 0
    while it < maxiter:
        tol1 = tol * abs(x) + mintol
        tol2 = 2.0 * tol1
        xmid = 0.5 * (a + b)
        if abs(x - xmid) < (tol2 - 0.5 * (b - a)):          # convergência
            break
        if abs(deltax) <= tol1:
            deltax = (a - x) if x >= xmid else (b - x)       # secção áurea
            rat = cg * deltax
        else:                                                # parábola
            tmp1 = (x - w) * (fx - fv)
            tmp2 = (x - v) * (fx - fw)
            p = (x - v) * tmp2 - (x - w) * tmp1
            tmp2 = 2.0 * (tmp2 - tmp1)
            if tmp2 > 0.0:
                p = -p
            tmp2 = abs(tmp2)
            dx_temp = deltax
            deltax = rat
            if (p > tmp2 * (a - x)) and (p < tmp2 * (b - x)) and \
                    (abs(p) < abs(0.5 * tmp2 * dx_temp)):
                rat = p * 1.0 / tmp2
                u = x + rat
                if (u - a) < tol2 or (b - u) < tol2:
                    rat = tol1 if xmid - x >= 0 else -tol1
            else:
                deltax = (a - x) if x >= xmid else (b - x)
                rat = cg * deltax
        if abs(rat) < tol1:                                  # andar pelo menos tol1
            u = x + tol1 if rat >= 0 else x - tol1
        else:
            u = x + rat
        fu = F(u)
        if fu > fx:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v = w; w = u; fv = fw; fw = fu
            elif fu <= fv or v == x or v == w:
                v = u; fv = fu
        else:
            if u >= x:
                a = x
            else:
                b = x
            v = w; w = x; x = u
            fv = fw; fw = fx; fx = fu
        it += 1
    return x


def _fminbnd(F, a, b, tolx, maxfun=500):
    """Algoritmo do fminbnd do MATLAB em [a, b] (secção áurea + parábola)."""
    seps = math.sqrt(2.220446049250313e-16)
    c = 0.5 * (3.0 - math.sqrt(5.0))
    v = a + c * (b - a)
    w = xf = v
    e = rat = 0.0
    fx = F(xf)
    num = 1
    fv = fw = fx
    xm = 0.5 * (a + b)
    tol1 = seps * abs(xf) + tolx / 3.0
    tol2 = 2.0 * tol1
    while abs(xf - xm) > (tol2 - 0.5 * (b - a)):
        golden = True
        if abs(e) > tol1:                                    # tentar a parábola
            golden = False
            r = (xf - w) * (fx - fv)
            q = (xf - v) * (fx - fw)
            p = (xf - v) * q - (xf - w) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q = abs(q)
            r = e
            e = rat
            if abs(p) < abs(0.5 * q * r) and p > q * (a - xf) and p < q * (b - xf):
                rat = p / q
                x = xf + rat
                if (x - a) < tol2 or (b - x) < tol2:
                    rat = tol1 * (math.copysign(1.0, xm - xf) if xm != xf else 1.0)
            else:
                golden = True
        if golden:                                           # secção áurea
            e = (a - xf) if xf >= xm else (b - xf)
            rat = c * e
        si = math.copysign(1.0, rat) if rat != 0 else 1.0
        x = xf + si * max(abs(rat), tol1)
        fu = F(x)
        num += 1
        if fu <= fx:
            if x >= xf:
                a = xf
            else:
                b = xf
            v, fv = w, fw
            w, fw = xf, fx
            xf, fx = x, fu
        else:
            if x < xf:
                a = x
            else:
                b = x
            if fu <= fw or w == xf:
                v, fv = w, fw
                w, fw = x, fu
            elif fu <= fv or v == xf or v == w:
                v, fv = x, fu
        xm = 0.5 * (a + b)
        tol1 = seps * abs(xf) + tolx / 3.0
        tol2 = 2.0 * tol1
        if num >= maxfun:
            break
    return xf
