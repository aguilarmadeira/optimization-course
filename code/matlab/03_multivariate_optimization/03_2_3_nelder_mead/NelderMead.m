function [x, fx, info] = NelderMead(f, X0, tolx, tolf, kmax, ftarget, verbose, opts)
%NELDERMEAD  Método do simplex de Nelder–Mead: min f(x) em R^n, só com valores de f.
%
%   [x, fx, info] = NelderMead(f, X0)
%   [x, fx, info] = NelderMead(f, X0, tolx, tolf, kmax, ftarget, verbose)
%   [x, fx, info] = NelderMead(f, X0, [], [], kmax, [], verbose, opts)
%
%   Notação das aulas: xl = melhor vértice, xg = o seguinte ao pior,
%   xh = pior; xc = centroide dos vértices exceto xh.
%     reflexão            xr = 2 xc - xh
%     f(xr) < f(xl)       expansão: (1+gamma) xc - gamma xh; fica se f < f(xr), senão xr
%     f(xr) < f(xg)       aceitar xr
%     f(xr) < f(xh)       contração exterior: (1+beta) xc - beta xh;
%                         fica se f <= f(xr), senão encolher
%     f(xr) >= f(xh)      contração interior: (1-beta) xc + beta xh;
%                         fica se f < f(xh), senão encolher
%     encolher            xi = xl + (xi - xl)/2 para todo i ~= l
%   Com opts.regra = 'deb' (Deb, 2012, sec. 3.3.2) o ponto novo substitui
%   sempre xh e não há encolhimento.
%
%   Entradas
%     f        função (handle) de x, vetor linha com n componentes
%     X0       ponto inicial x0 (vetor linha) -> simplex x0, x0 + h_i e_i,
%              com h_i = 0.05*max(1,|x0_i|); ou matriz (n+1) x n com os
%              vértices do simplex inicial, um por linha
%     tolx     (opcional, 1e-6) paragem do Nelder–Mead padrão (como o
%     tolf     (opcional, 1e-6)   fminsearch), usada se opts.eps não for dado:
%                                 f(xh) - f(xl) <= tolf e max_i ||x_i - xl|| <= tolx
%     kmax     (opcional, 500) número máximo de iterações
%     ftarget  (opcional, -Inf) alvo para info.nhit
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%     opts     (opcional) estrutura com
%       .gamma   parâmetro de expansão (2)
%       .beta    parâmetro de contração (0.5)
%       .eps     se for dado, pára quando Q <= eps, com a análise do erro das aulas
%                Q = sqrt( sum_i (f(x_i) - f(xc))^2 / (n+1) ), calculada em cada
%                iteração com o novo simplex e o centroide dessa iteração
%                (f(xc) conta como uma avaliação)
%       .regra   'padrao' (por omissão) ou 'deb'
%
%   Saídas
%     x        melhor vértice do simplex final (xl)
%     fx       f(xl)
%     info     estrutura com
%       .nfev      avaliações de f, incluindo as n+1 do simplex inicial
%       .nhit      primeira avaliação com f < ftarget (Inf se nunca)
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       iterações feitas
%       .history   uma linha por k = 0, 1, ..., nit (depois da iteração k):
%                  [k  n_f  f(xl)  f(xh)  Q  xl]   (Q = NaN sem opts.eps)
%       .ops       operação de cada linha de .history
%       .cols      nomes das colunas de .history
%       .simplex   vértices do simplex final, ordenados (o melhor primeiro)
%       .Q         último valor de Q (NaN sem opts.eps)
%       .flag      0 (critério de paragem atingido) ou 1 (atingiu kmax)
%       .message   mensagem de paragem
%
%   Otimização — deck 3.2.3 (Método do simplex de Nelder–Mead).
%   Reproduz os exemplos dos slides: ver ex03_2_3_nelder_mead.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não deteta simplex degenerado nem reinicia).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 3 || isempty(tolx), tolx = 1e-6; end
if nargin < 4 || isempty(tolf), tolf = 1e-6; end
if nargin < 5 || isempty(kmax), kmax = 500; end
if nargin < 6 || isempty(ftarget), ftarget = -Inf; end
if nargin < 7 || isempty(verbose), verbose = false; end
if nargin < 8, opts = struct(); end
gamma = 2;  beta = 1/2;  eps_Q = [];  regra = 'padrao';
if isfield(opts, 'gamma'), gamma = opts.gamma; end
if isfield(opts, 'beta'),  beta = opts.beta;  end
if isfield(opts, 'eps'),   eps_Q = opts.eps;  end
if isfield(opts, 'regra'), regra = opts.regra; end
if ~any(strcmp(regra, {'padrao', 'deb'})), error('NelderMead: opts.regra tem de ser ''padrao'' ou ''deb''.'); end
deb = strcmp(regra, 'deb');  usaQ = ~isempty(eps_Q);
sigma = 1/2;                                        % encolhimento para metade

if min(size(X0)) == 1                               % simplex inicial a partir de x0
  x0 = X0(:)';  n = numel(x0);  h = 0.05*max(1, abs(x0));
  X = [x0; repmat(x0, n, 1) + diag(h)];
else
  X = X0;  n = size(X, 2);
  if size(X, 1) ~= n + 1, error('NelderMead: o simplex tem de ter n+1 vértices.'); end
end

nfev = 0;  nhit = Inf;  Q = NaN;

% --- núcleo (o algoritmo do deck, com as contagens) ----------------------
F = zeros(n + 1, 1);
for i = 1:n + 1
  [F(i), nfev, nhit] = aval(f, X(i, :), nfev, nhit, ftarget);   % passo 1: n+1 avaliações
end
[fl, il] = min(F);
hist = [0, nfev, fl, max(F), NaN, X(il, :)];  ops = {'inicial'};
flag = 1;
for k = 1:kmax
  [F, i] = sort(F);  X = X(i, :);                   % passo 2: ordenar
  % xl = X(1,:), xg = X(n,:), xh = X(n+1,:)
  if ~usaQ && F(n+1) - F(1) <= tolf && max(sqrt(sum((X - repmat(X(1, :), n + 1, 1)).^2, 2))) <= tolx
    flag = 0;                                       % paragem padrão (antes da iteração k)
    break
  end
  xh = X(n+1, :);
  xc = mean(X(1:n, :), 1);                          % centroide sem o pior
  if usaQ, [fc, nfev, nhit] = aval(f, xc, nfev, nhit, ftarget); end   % f(xc), para Q
  xr = xc + (xc - xh);  [fr, nfev, nhit] = aval(f, xr, nfev, nhit, ftarget);   % passo 3: xr = 2xc - xh
  op = 'reflexão';  encolher = false;
  if fr < F(1)                                      % expansão: (1+gamma)xc - gamma xh
    xe = xc + gamma*(xr - xc);  [fe, nfev, nhit] = aval(f, xe, nfev, nhit, ftarget);
    if deb || fe < fr
      X(n+1, :) = xe;  F(n+1) = fe;  op = 'expansão';
    else
      X(n+1, :) = xr;  F(n+1) = fr;
    end
  elseif fr < F(n)                                  % aceitar a reflexão
    X(n+1, :) = xr;  F(n+1) = fr;
  elseif fr < F(n+1)                                % contração exterior: (1+beta)xc - beta xh
    xo = xc + beta*(xr - xc);  [fo, nfev, nhit] = aval(f, xo, nfev, nhit, ftarget);
    op = 'contr. ext.';
    if deb || fo <= fr, X(n+1, :) = xo;  F(n+1) = fo;  else, encolher = true;  end
  else                                              % contração interior: (1-beta)xc + beta xh
    xi = xc - beta*(xc - xh);  [fi, nfev, nhit] = aval(f, xi, nfev, nhit, ftarget);
    op = 'contr. int.';
    if deb || fi < F(n+1), X(n+1, :) = xi;  F(n+1) = fi;  else, encolher = true;  end
  end
  if encolher                                       % xi = xl + sigma(xi - xl), i ~= l
    op = 'encolhimento';
    for j = 2:n + 1
      X(j, :) = X(1, :) + sigma*(X(j, :) - X(1, :));
      [F(j), nfev, nhit] = aval(f, X(j, :), nfev, nhit, ftarget);
    end
  end
  if usaQ, Q = sqrt(mean((F - fc).^2)); end         % passo 4: análise do erro
  [fl, il] = min(F);
  hist(end + 1, :) = [k, nfev, fl, max(F), Q, X(il, :)];  ops{end + 1} = op;
  if usaQ && Q <= eps_Q
    flag = 0;
    break
  end
end
% -------------------------------------------------------------------------

[F, i] = sort(F);  X = X(i, :);
x = X(1, :);  fx = F(1);
cols = [{'k', 'n_f', 'f(x_l)', 'f(x_h)', 'Q'}, cell(1, n)];
for j = 1:n, cols{5 + j} = sprintf('x_l%d', j); end

info.nfev = nfev;
info.nhit = nhit;
info.ngev = 0;
info.nhev = 0;
info.nit = size(hist, 1) - 1;
info.history = hist;
info.ops = ops;
info.cols = cols;
info.simplex = X;
info.Q = Q;
info.flag = flag;
if flag == 0
  if usaQ, info.message = 'Q <= eps';
  else, info.message = 'f(x_h) - f(x_l) <= tolf e diâmetro <= tolx'; end
else
  info.message = sprintf('atingiu kmax = %d iterações', kmax);
end

if verbose
  fprintf('%4s %5s %12s %12s %10s', cols{1:5});  fprintf(' %10s', cols{6:end});
  fprintf('  operação\n');
  for i = 1:size(hist, 1)
    fprintf('%4d %5d %12.4f %12.4f %10.4f', hist(i, 1:5));  fprintf(' %10.4f', hist(i, 6:end));
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
