function [alpha, nfev, fvals] = LineSearch(phi, metodo, tol, amax)
%LINESEARCH  Pesquisa em linha para os métodos de 3.3: alpha ~ argmin phi(alpha), alpha >= 0.
%
%   [alpha, nfev, fvals] = LineSearch(phi)
%   [alpha, nfev, fvals] = LineSearch(phi, metodo, tol, amax)
%
%   Entradas
%     phi      função de uma variável (handle), phi(alpha) = f(x_k + alpha*d_k)
%     metodo   (opcional) 'brent' (por omissão) ou 'fminbnd'
%       'brent'   a pesquisa de elevada precisão das figuras e da tabela
%                 comparativa do capítulo 3: enquadramento a partir de
%                 (0, 1e-3) (passos em razão áurea + interpolação parabólica)
%                 e depois o método de Brent (secção áurea + interpolação
%                 parabólica, 2.2-2.3) com tolerância relativa tol = 1e-10.
%                 É uma tradução linha a linha de scipy.optimize.bracket +
%                 Brent (minimize_scalar, method='brent'): dá os mesmos
%                 pontos e as mesmas contagens que os scripts das figuras
%                 (cerca de 24 avaliações por pesquisa na quadrática do 3.3.1).
%       'fminbnd' o algoritmo do fminbnd do MATLAB (Brent em [0, amax],
%                 TolX = 1e-4 por omissão), com o alargamento do slide
%                 «No computador»: se alpha ~ amax, repetir com 10*amax.
%                 Na quadrática do 3.3.1 gasta 6 avaliações por pesquisa.
%                 Não usa o fminbnd: corre igual em MATLAB e Octave.
%     tol      (opcional) 'brent': tolerância relativa (1e-10);
%              'fminbnd': TolX (1e-4)
%     amax     (opcional, 1) 'fminbnd': extremo direito inicial
%
%   Saídas
%     alpha    passo, alpha >= 0
%     nfev     número de avaliações de phi
%     fvals    valores de phi pela ordem em que foram calculados
%
%   Nota: 'brent' volta a avaliar phi(0) = f(x_k) (é a primeira avaliação do
%   enquadramento) e essa avaliação conta.
%
%   Otimização — decks 3.3.1, 3.3.2 e 3.3.3 (pesquisa em linha; cap. 2).
%   Utilitário comum (code/matlab/common/): usado por SteepestDescent (3.3.1),
%   NewtonND (3.3.2), ConjGrad (3.3.3) e pela comparação do cap. 3; existe só
%   aqui. Os exemplos acrescentam-no ao caminho com uc_setup.
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. condições de Wolfe).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 2 || isempty(metodo), metodo = 'brent'; end
if nargin < 3, tol = []; end
if nargin < 4 || isempty(amax), amax = 1; end

switch metodo
  case 'brent'
    if isempty(tol), tol = 1e-10; end
    [alpha, fvals] = ls_brent(phi, tol);
    alpha = max(alpha, 0);
  case 'fminbnd'
    if isempty(tol), tol = 1e-4; end
    b = amax;
    [alpha, fvals] = ls_fminbnd(phi, 0, b, tol, []);
    while alpha > 0.99*b              % mínimo no extremo: alargar o intervalo
      b = 10*b;
      [alpha, fvals] = ls_fminbnd(phi, 0, b, tol, fvals);
    end
  otherwise
    error('LineSearch: método desconhecido ''%s''', metodo);
end
nfev = numel(fvals);
end

% =========================================================================
function [xa, xb, xc, fa, fb, fc, fv] = ls_bracket(phi)
% Enquadramento (xa, xb, xc) com phi(xb) <= phi(xa), phi(xc)
% (tradução de scipy.optimize.bracket, a partir de xa = 0, xb = 1e-3).
gold = 1.618034;  verysmall = 1e-21;  grow_limit = 110.0;  maxiter = 1000;
xa = 0.0;  xb = 1e-3;
fa = phi(xa);  fb = phi(xb);  fv = [fa, fb];
if fa < fb                              % trocar para que fa > fb
  t = xa; xa = xb; xb = t;  t = fa; fa = fb; fb = t;
end
xc = xb + gold*(xb - xa);
fc = phi(xc);  fv(end + 1) = fc;
it = 0;
while fc < fb
  tmp1 = (xb - xa)*(fb - fc);
  tmp2 = (xb - xc)*(fb - fa);
  val = tmp2 - tmp1;
  if abs(val) < verysmall, denom = 2.0*verysmall; else, denom = 2.0*val; end
  w = xb - ((xb - xc)*tmp2 - (xb - xa)*tmp1)/denom;      % parábola
  wlim = xb + grow_limit*(xc - xb);
  if it > maxiter, error('LineSearch: não encontrou enquadramento'); end
  it = it + 1;
  if (w - xc)*(xb - w) > 0.0
    fw = phi(w);  fv(end + 1) = fw;
    if fw < fc
      xa = xb; xb = w; fa = fb; fb = fw;
      break
    elseif fw > fb
      xc = w; fc = fw;
      break
    end
    w = xc + gold*(xc - xb);
    fw = phi(w);  fv(end + 1) = fw;
  elseif (w - wlim)*(wlim - xc) >= 0.0
    w = wlim;
    fw = phi(w);  fv(end + 1) = fw;
  elseif (w - wlim)*(xc - w) > 0.0
    fw = phi(w);  fv(end + 1) = fw;
    if fw < fc
      xb = xc; xc = w;
      w = xc + gold*(xc - xb);
      fb = fc; fc = fw;
      fw = phi(w);  fv(end + 1) = fw;
    end
  else
    w = xc + gold*(xc - xb);
    fw = phi(w);  fv(end + 1) = fw;
  end
  xa = xb; xb = xc; xc = w;
  fa = fb; fb = fc; fc = fw;
end
ok = ((fb < fc && fb <= fa) || (fb < fa && fb <= fc)) && ...
     ((xa < xb && xb < xc) || (xc < xb && xb < xa)) && all(isfinite([xa xb xc]));
if ~ok, error('LineSearch: não encontrou enquadramento válido'); end
end

% =========================================================================
function [x, fv] = ls_brent(phi, tol)
% Método de Brent a partir do enquadramento (tradução de scipy.optimize.Brent).
[xa, xb, xc, fa, fb, fc, fv] = ls_bracket(phi); %#ok<ASGLU>
mintol = 1.0e-11;  cg = 0.3819660;  maxiter = 500;
x = xb; w = xb; v = xb;
fx = fb; fw = fb; fvv = fb;
if xa < xc, a = xa; b = xc; else, a = xc; b = xa; end
deltax = 0.0;  rat = 0.0;  it = 0;
while it < maxiter
  tol1 = tol*abs(x) + mintol;
  tol2 = 2.0*tol1;
  xmid = 0.5*(a + b);
  if abs(x - xmid) < (tol2 - 0.5*(b - a))           % convergência
    break
  end
  if abs(deltax) <= tol1
    if x >= xmid, deltax = a - x; else, deltax = b - x; end   % secção áurea
    rat = cg*deltax;
  else                                              % parábola
    tmp1 = (x - w)*(fx - fvv);
    tmp2 = (x - v)*(fx - fw);
    p = (x - v)*tmp2 - (x - w)*tmp1;
    tmp2 = 2.0*(tmp2 - tmp1);
    if tmp2 > 0.0, p = -p; end
    tmp2 = abs(tmp2);
    dx_temp = deltax;
    deltax = rat;
    if (p > tmp2*(a - x)) && (p < tmp2*(b - x)) && (abs(p) < abs(0.5*tmp2*dx_temp))
      rat = p*1.0/tmp2;
      u = x + rat;
      if (u - a) < tol2 || (b - u) < tol2
        if xmid - x >= 0, rat = tol1; else, rat = -tol1; end
      end
    else
      if x >= xmid, deltax = a - x; else, deltax = b - x; end
      rat = cg*deltax;
    end
  end
  if abs(rat) < tol1                                % andar pelo menos tol1
    if rat >= 0, u = x + tol1; else, u = x - tol1; end
  else
    u = x + rat;
  end
  fu = phi(u);  fv(end + 1) = fu; %#ok<AGROW>
  if fu > fx
    if u < x, a = u; else, b = u; end
    if fu <= fw || w == x
      v = w; w = u; fvv = fw; fw = fu;
    elseif fu <= fvv || v == x || v == w
      v = u; fvv = fu;
    end
  else
    if u >= x, a = x; else, b = x; end
    v = w; w = x; x = u;
    fvv = fw; fw = fx; fx = fu;
  end
  it = it + 1;
end
end

% =========================================================================
function [xf, fv] = ls_fminbnd(phi, a, b, tolx, fv)
% Algoritmo do fminbnd do MATLAB em [a, b] (secção áurea + parábola).
seps = sqrt(eps);  c = 0.5*(3.0 - sqrt(5.0));  maxfun = 500;
v = a + c*(b - a);  w = v;  xf = v;
e = 0.0;  rat = 0.0;
fx = phi(xf);  fv(end + 1) = fx;  num = 1;
fvv = fx;  fw = fx;
xm = 0.5*(a + b);
tol1 = seps*abs(xf) + tolx/3.0;  tol2 = 2.0*tol1;
while abs(xf - xm) > (tol2 - 0.5*(b - a))
  golden = true;
  if abs(e) > tol1                                  % tentar a parábola
    golden = false;
    r = (xf - w)*(fx - fvv);
    q = (xf - v)*(fx - fw);
    p = (xf - v)*q - (xf - w)*r;
    q = 2.0*(q - r);
    if q > 0.0, p = -p; end
    q = abs(q);
    r = e;
    e = rat;
    if abs(p) < abs(0.5*q*r) && p > q*(a - xf) && p < q*(b - xf)
      rat = p/q;
      x = xf + rat;
      if (x - a) < tol2 || (b - x) < tol2
        si = sign(xm - xf) + (xm == xf);
        rat = tol1*si;
      end
    else
      golden = true;
    end
  end
  if golden                                         % secção áurea
    if xf >= xm, e = a - xf; else, e = b - xf; end
    rat = c*e;
  end
  si = sign(rat) + (rat == 0);
  x = xf + si*max(abs(rat), tol1);
  fu = phi(x);  fv(end + 1) = fu; %#ok<AGROW>
  num = num + 1;
  if fu <= fx
    if x >= xf, a = xf; else, b = xf; end
    v = w; fvv = fw;
    w = xf; fw = fx;
    xf = x; fx = fu;
  else
    if x < xf, a = x; else, b = x; end
    if fu <= fw || w == xf
      v = w; fvv = fw;
      w = x; fw = fu;
    elseif fu <= fvv || v == xf || v == w
      v = x; fvv = fu;
    end
  end
  xm = 0.5*(a + b);
  tol1 = seps*abs(xf) + tolx/3.0;  tol2 = 2.0*tol1;
  if num >= maxfun, break, end
end
end
