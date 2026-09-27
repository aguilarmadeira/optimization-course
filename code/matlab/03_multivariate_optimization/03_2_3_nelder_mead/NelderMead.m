function [x, fx, info] = NelderMead(f, X0, tolx, tolf, kmax, ftarget, verbose)
%NELDERMEAD  Método de Nelder–Mead (simplex): min f(x) em R^n, só com valores de f.
%
%   [x, fx, info] = NelderMead(f, X0)
%   [x, fx, info] = NelderMead(f, X0, tolx, tolf, kmax, ftarget, verbose)
%
%   Coeficientes (alpha, gamma, rho, sigma) = (1, 2, 1/2, 1/2) e regras de
%   aceitação/encolhimento do pseudocódigo do deck («O algoritmo»):
%     f(xr) < f(xb)           expandir: xe = xc + gamma(xr - xc); fica o melhor de xe, xr
%     f(xr) < f(xs)           aceitar xr
%     f(xr) < f(xw)           contração exterior: xoc = xc + rho(xr - xc);
%                             aceitar se f(xoc) <= f(xr), senão encolher
%     f(xr) >= f(xw)          contração interior: xic = xc - rho(xc - xw);
%                             aceitar se f(xic) < f(xw), senão encolher
%     encolher                xi <- xb + sigma(xi - xb) para todo i ~= b
%
%   Entradas
%     f        função (handle) de x, vetor linha com n componentes
%     X0       ponto inicial x0 (vetor linha) -> simplex x0, x0 + h_i e_i,
%              com h_i = 0.05*max(1,|x0_i|) (o do deck); ou matriz (n+1) x n
%              com os vértices do simplex inicial, um por linha (nos
%              exemplos do deck: arestas 1.2 no Himmelblau e 0.5 no Rosenbrock)
%     tolx     (opcional, 1e-6) para quando f(xw) - f(xb) <= tolf e o simplex
%     tolf     (opcional, 1e-6)   tem diâmetro <= tolx, medido como
%                                 max_i ||x_i - xb||
%     kmax     (opcional, 500) número máximo de iterações
%     ftarget  (opcional, -Inf) alvo para info.nhit
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        melhor vértice do simplex final (xb)
%     fx       f(xb)
%     info     estrutura com
%       .nfev      avaliações de f, incluindo as n+1 do simplex inicial:
%                  1-2 por iteração, n+2 numa iteração com encolhimento
%       .nhit      primeira avaliação com f < ftarget (Inf se nunca)
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       iterações feitas
%       .history   uma linha por k = 0, 1, ..., nit (depois da iteração k):
%                  [k  n_f  f(xb)  f(xw)  xb]
%       .ops       operação de cada linha de .history ('inicial', 'reflexão',
%                  'expansão', 'contr. ext.', 'contr. int.', 'encolhimento')
%       .cols      nomes das colunas de .history
%       .simplex   vértices do simplex final, ordenados (o melhor primeiro)
%       .flag      0 (tolerâncias atingidas) ou 1 (atingiu kmax)
%       .message   mensagem de paragem
%
%   Otimização — deck 3.2.3 (Método de Nelder–Mead). O deck usa fminsearch;
%   esta é a implementação da UC, com o pseudocódigo do deck.
%   Reproduz os exemplos dos slides: ver ex03_2_3_nelder_mead.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não deteta simplex degenerado nem reinicia).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 3 || isempty(tolx), tolx = 1e-6; end
if nargin < 4 || isempty(tolf), tolf = 1e-6; end
if nargin < 5 || isempty(kmax), kmax = 500; end
if nargin < 6 || isempty(ftarget), ftarget = -Inf; end
if nargin < 7, verbose = false; end
alpha = 1;  gamma = 2;  rho = 1/2;  sigma = 1/2;    % coeficientes originais

if min(size(X0)) == 1                               % simplex inicial do deck
  x0 = X0(:)';  n = numel(x0);  h = 0.05*max(1, abs(x0));
  X = [x0; repmat(x0, n, 1) + diag(h)];
else
  X = X0;  n = size(X, 2);
  if size(X, 1) ~= n + 1, error('NelderMead: o simplex tem de ter n+1 vértices.'); end
end

nfev = 0;  nhit = Inf;

% --- núcleo (o pseudocódigo do deck, com as contagens) -------------------
F = zeros(n + 1, 1);
for i = 1:n + 1
  [F(i), nfev, nhit] = aval(f, X(i, :), nfev, nhit, ftarget);   % n+1 avaliações
end
[fb, ib] = min(F);
hist = [0, nfev, fb, max(F), X(ib, :)];  ops = {'inicial'};
flag = 1;
for k = 1:kmax
  [F, i] = sort(F);  X = X(i, :);                   % ordenar
  % xb = X(1,:), xs = X(n,:), xw = X(n+1,:)
  if F(n+1) - F(1) <= tolf && max(sqrt(sum((X - repmat(X(1, :), n + 1, 1)).^2, 2))) <= tolx
    flag = 0;                                       % parou (antes de fazer a iteração k)
    break
  end
  xw = X(n+1, :);
  xc = mean(X(1:n, :), 1);                          % centroide sem o pior
  xr = xc + alpha*(xc - xw);  [fr, nfev, nhit] = aval(f, xr, nfev, nhit, ftarget);
  op = 'reflexão';  encolher = false;
  if fr < F(1)                                      % expandir
    xe = xc + gamma*(xr - xc);  [fe, nfev, nhit] = aval(f, xe, nfev, nhit, ftarget);
    if fe < fr
      X(n+1, :) = xe;  F(n+1) = fe;  op = 'expansão';
    else
      X(n+1, :) = xr;  F(n+1) = fr;
    end
  elseif fr < F(n)                                  % aceitar
    X(n+1, :) = xr;  F(n+1) = fr;
  elseif fr < F(n+1)                                % contração exterior
    xoc = xc + rho*(xr - xc);  [foc, nfev, nhit] = aval(f, xoc, nfev, nhit, ftarget);
    op = 'contr. ext.';
    if foc <= fr, X(n+1, :) = xoc;  F(n+1) = foc;  else, encolher = true;  end
  else                                              % contração interior
    xic = xc - rho*(xc - xw);  [fic, nfev, nhit] = aval(f, xic, nfev, nhit, ftarget);
    op = 'contr. int.';
    if fic < F(n+1), X(n+1, :) = xic;  F(n+1) = fic;  else, encolher = true;  end
  end
  if encolher                                       % xi <- xb + sigma(xi - xb), i ~= b
    op = 'encolhimento';
    for j = 2:n + 1
      X(j, :) = X(1, :) + sigma*(X(j, :) - X(1, :));
      [F(j), nfev, nhit] = aval(f, X(j, :), nfev, nhit, ftarget);
    end
  end
  [fb, ib] = min(F);
  hist(end + 1, :) = [k, nfev, fb, max(F), X(ib, :)];  ops{end + 1} = op;
end
% -------------------------------------------------------------------------

[F, i] = sort(F);  X = X(i, :);
x = X(1, :);  fx = F(1);
cols = [{'k', 'n_f', 'f(x_b)', 'f(x_w)'}, cell(1, n)];
for j = 1:n, cols{4 + j} = sprintf('x_b%d', j); end

info.nfev = nfev;
info.nhit = nhit;
info.ngev = 0;
info.nhev = 0;
info.nit = size(hist, 1) - 1;
info.history = hist;
info.ops = ops;
info.cols = cols;
info.simplex = X;
info.flag = flag;
if flag == 0
  info.message = 'f(x_w) - f(x_b) <= tolf e diâmetro <= tolx';
else
  info.message = sprintf('atingiu kmax = %d iterações', kmax);
end

if verbose
  fprintf('%4s %5s %12s %12s', cols{1:4});  fprintf(' %10s', cols{5:end});
  fprintf('  operação\n');
  for i = 1:size(hist, 1)
    fprintf('%4d %5d %12.4f %12.4f', hist(i, 1:4));  fprintf(' %10.4f', hist(i, 5:end));
    fprintf('  %s\n', ops{i});
  end
end
end

% -------------------------------------------------------------------------
function [v, nfev, nhit] = aval(f, x, nfev, nhit, ftarget)
% avalia f em x e atualiza as contagens (n_f e a 1.ª avaliação abaixo do alvo)
v = f(x);  nfev = nfev + 1;
if v < ftarget && isinf(nhit), nhit = nfev; end
end
