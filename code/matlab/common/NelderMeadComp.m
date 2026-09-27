function [x, fx, info] = NelderMeadComp(f, x0, tolx, tolf, kmax)
%NELDERMEADCOMP  Nelder–Mead da comparação de 4.3.2 (minimizador interno dos exemplos de 4.3).
%
%   [x, fx, info] = NelderMeadComp(f, x0)
%   [x, fx, info] = NelderMeadComp(f, x0, tolx, tolf, kmax)
%
%   Utilitário comum (code/matlab/common/), usado por PenaltyExterior (4.3.1),
%   BarrierLog (4.3.2) e pelos dois exemplos. É a VARIANTE do Nelder–Mead
%   usada na comparação do 4.3: a contração exterior só é aceite com «<»
%   (o NelderMead do deck 3.2.3 aceita com «<=»). Existe como função à
%   parte, e não como opção de NelderMead, para que o NelderMead do deck
%   3.2.3 fique exatamente o do pseudocódigo dos slides.
%
%   Variante de NelderMead (deck 3.2.3, pasta
%   03_multivariate_optimization/03_2_3_nelder_mead/), para reproduzir
%   EXATAMENTE o Nelder–Mead didático do script da comparação exterior vs.
%   interior (make_figs_4comp.py), que produziu os números do slide
%   «Exterior vs. interior» (844, 743 chamadas a P, n_f = 717, 26 fora de X;
%   864 de (1,1)).
%
%   O algoritmo é o do deck 3.2.3 — coeficientes (alpha, gamma, rho, sigma)
%   = (1, 2, 1/2, 1/2), simplex inicial x0, x0 + h_i e_i com
%   h_i = 0.05*max(1,|x0_i|), mesma paragem — com duas diferenças, ambas do
%   script da comparação:
%     1. contração exterior aceite só se f(xoc) < f(xr) (desigualdade
%        ESTRITA; o pseudocódigo do deck 3.2.3 aceita com <=). Perto da
%        convergência, com tolx = 1e-10, há empates exatos f(xoc) == f(xr)
%        e a regra de desempate muda o número de avaliações;
%     2. os pontos calculam-se como no script: xe = xc + 2(xc - xw),
%        xoc = xc + (xc - xw)/2 (em aritmética exata é o mesmo que
%        xc + gamma(xr - xc) e xc + rho(xr - xc); os arredondamentos diferem).
%   Aceita f = Inf (pontos fora de X na barreira).
%
%   Entradas
%     f        função (handle) de x, vetor linha; pode devolver Inf
%     x0       ponto inicial (vetor)
%     tolx     (opcional, 1e-10) para quando f(xw) - f(xb) <= tolf e
%     tolf     (opcional, 1e-12)   max_i ||x_i - xb|| <= tolx
%     kmax     (opcional, 2000) número máximo de iterações
%
%   Saídas
%     x, fx    melhor vértice do simplex final e f(x)
%     info     .nfev (inclui as n+1 avaliações do simplex inicial), .ngev,
%              .nhev (0), .nit, .history ([k n_f f(xb) f(xw) xb], uma linha
%              por k = 0..nit), .cols, .flag (0 tolerâncias, 1 kmax), .message
%
%   Otimização — decks 4.3.1 e 4.3.2. Implementação didática — algumas
%   salvaguardas de software profissional não estão incluídas (p. ex. não
%   deteta simplex degenerado nem reinicia).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 3 || isempty(tolx), tolx = 1e-10; end
if nargin < 4 || isempty(tolf), tolf = 1e-12; end
if nargin < 5 || isempty(kmax), kmax = 2000; end

x0 = x0(:)';  n = numel(x0);  h = 0.05*max(1, abs(x0));
X = [x0; repmat(x0, n, 1) + diag(h)];
nfev = 0;

% --- núcleo (deck 3.2.3, com as duas diferenças do script) ---------------
F = zeros(n + 1, 1);
for i = 1:n + 1, F(i) = f(X(i, :)); end
nfev = nfev + n + 1;
[fb, ib] = min(F);
hist = [0, nfev, fb, max(F), X(ib, :)];
flag = 1;
for k = 1:kmax
  [F, i] = sort(F);  X = X(i, :);                   % ordenar
  if F(n+1) - F(1) <= tolf && max(sqrt(sum((X - repmat(X(1, :), n + 1, 1)).^2, 2))) <= tolx
    flag = 0;
    break
  end
  xw = X(n+1, :);
  xc = mean(X(1:n, :), 1);                          % centroide sem o pior
  xr = xc + (xc - xw);  fr = f(xr);  nfev = nfev + 1;   % reflexão
  if fr < F(1)                                      % expandir
    xe = xc + 2*(xc - xw);  fe = f(xe);  nfev = nfev + 1;
    if fe < fr, X(n+1, :) = xe;  F(n+1) = fe;  else, X(n+1, :) = xr;  F(n+1) = fr;  end
  elseif fr < F(n)                                  % aceitar
    X(n+1, :) = xr;  F(n+1) = fr;
  else                                              % contrair
    if fr < F(n+1)
      xt = xc + 0.5*(xc - xw);                      % contração exterior
    else
      xt = xc - 0.5*(xc - xw);                      % contração interior
    end
    ft = f(xt);  nfev = nfev + 1;
    if ft < min(fr, F(n+1))                         % estrito (script)
      X(n+1, :) = xt;  F(n+1) = ft;
    else                                            % encolher
      for j = 2:n + 1
        X(j, :) = X(1, :) + 0.5*(X(j, :) - X(1, :));  F(j) = f(X(j, :));
      end
      nfev = nfev + n;
    end
  end
  [fb, ib] = min(F);
  hist(end + 1, :) = [k, nfev, fb, max(F), X(ib, :)];
end
% -------------------------------------------------------------------------

[F, i] = sort(F);  X = X(i, :);
x = X(1, :);  fx = F(1);
cols = [{'k', 'n_f', 'f(x_b)', 'f(x_w)'}, cell(1, n)];
for j = 1:n, cols{4 + j} = sprintf('x_b%d', j); end
info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = size(hist, 1) - 1;
info.history = hist;
info.cols = cols;
info.flag = flag;
if flag == 0
  info.message = 'f(x_w) - f(x_b) <= tolf e diâmetro <= tolx';
else
  info.message = sprintf('atingiu kmax = %d iterações', kmax);
end
end
