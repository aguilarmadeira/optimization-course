function [x, fx, info] = ConjGrad(f, grad, x0, tolg, kmax, opts)
%CONJGRAD  Gradientes conjugados não lineares (Fletcher-Reeves) com reinícios.
%
%   [x, fx, info] = ConjGrad(f, grad, x0, tolg, kmax)
%   [x, fx, info] = ConjGrad(f, grad, x0, tolg, kmax, opts)
%
%   d_0 = -g_0,  d_{k+1} = -g_{k+1} + beta_{k+1} d_k,
%   beta_{k+1} = ||g_{k+1}||^2 / ||g_k||^2  (Fletcher-Reeves).
%
%   Notação das aulas e dos slides: S_k = direção (aqui d), lambda_k =
%   comprimento do passo (aqui alpha), grad f_k = gradiente (aqui g).
%
%   Entradas
%     f      função (handle) de R^n em R
%     grad   gradiente de f (handle que devolve um vetor coluna)
%     x0     ponto inicial (vetor coluna)
%     tolg   pára quando ||g_k|| <= tolg (flag 0)
%     kmax   número máximo de iterações (flag 1)
%     opts   (opcional) estrutura com os campos (todos opcionais)
%       .restart  reinício periódico, d = -g quando mod(k+1, restart) == 0;
%                 por omissão restart = n (teste 1 do slide); 0 = sem
%                 reinício periódico. O reinício quando d deixa de ser de
%                 descida (g'*d >= 0, teste 2) está sempre ativo.
%       .ls       pesquisa em linha: 'brent' (por omissão; tol 1e-10, a das
%                 figuras) ou 'fminbnd' (TolX 1e-4) -- ver LineSearch
%       .lstol    tolerância da pesquisa em linha ([] = a de LineSearch)
%       .verbose  (false) se true, imprime a tabela das iterações
%
%   Saídas
%     x      última iterada x_k
%     fx     f(x_k)
%     info   estrutura com
%       .nfev      avaliações de f: pesquisas em linha (cada uma volta a
%                  avaliar phi(0) = f(x_k) no modo 'brent') + f em cada
%                  iterada (x_0 incluído)
%       .nfev_ls   avaliações de f só nas pesquisas em linha
%       .ngev      gradientes: nit + 1
%       .nhev      0 (o método não usa a Hessiana)
%       .nit       número de iterações
%       .nrestart  reinícios feitos (d = -g depois de x_0)
%       .history   uma linha por k = 0, 1, ..., nit:
%                  [k  x_k'  f(x_k)  ||g_k||  alpha_{k-1}  beta_k  d_k']
%                  (beta_k e d_k da última linha não são calculados: NaN)
%       .cols      nomes das colunas de .history
%       .flag      0 (||g|| <= tolg) ou 1 (atingiu kmax)
%       .message   mensagem de paragem
%
%   Precisa de LineSearch.m (code/matlab/common/). Os exemplos tratam do
%   caminho; fora deles, chamar antes addpath('<repo>/code/matlab'); uc_setup
%
%   Otimização — deck 3.3.3 (Gradientes conjugados).
%   Reproduz os exemplos dos slides: ver ex03_3_3_conjugate_gradients.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. pesquisa em linha com condições de Wolfe forte).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 6, opts = struct(); end
x = x0(:);
n = numel(x);
restart = campo(opts, 'restart', n);
ls = campo(opts, 'ls', 'brent');
lstol = campo(opts, 'lstol', []);
verbose = campo(opts, 'verbose', false);

nfev = 0;  ngev = 0;  nfev_ls = 0;  nrestart = 0;
fx = f(x);  nfev = nfev + 1;                   % f na iterada x_0
% --- núcleo (igual ao dos slides, com as contagens) ----------------------
g = grad(x);  g = g(:);  ngev = ngev + 1;
d = -g;                                        % 1.a direção
hist = [0, x', fx, norm(g), NaN, NaN, d'];
k = 0;
if norm(g) <= tolg, flag = 0; else, flag = 1; end
while flag == 1 && k < kmax
  phi = @(a) f(x + a*d);                       % pesquisa em linha
  [alpha, nls] = LineSearch(phi, ls, lstol);
  nfev = nfev + nls;  nfev_ls = nfev_ls + nls;
  x = x + alpha*d;
  fx = f(x);  nfev = nfev + 1;                 % f na nova iterada
  gn = grad(x);  gn = gn(:);  ngev = ngev + 1;
  hist(k + 2, :) = [k + 1, x', fx, norm(gn), alpha, NaN, NaN(1, n)];
  if norm(gn) <= tolg, flag = 0; k = k + 1; break, end
  beta = (gn'*gn)/(g'*g);                      % Fletcher-Reeves
  d = -gn + beta*d;
  if (restart > 0 && mod(k + 1, restart) == 0) || gn'*d >= 0
    d = -gn;  beta = 0;  nrestart = nrestart + 1;   % reinício
  end
  hist(k + 2, n + 5) = beta;
  hist(k + 2, n + 6:end) = d';
  g = gn;  k = k + 1;
end
% -------------------------------------------------------------------------

info.nfev = nfev;
info.nfev_ls = nfev_ls;
info.ngev = ngev;
info.nhev = 0;
info.nit = k;
info.nrestart = nrestart;
info.history = hist;
xc = cell(1, n);  dc = cell(1, n);
for i = 1:n, xc{i} = sprintf('x%d', i); dc{i} = sprintf('d%d', i); end
info.cols = [{'k'}, xc, {'f(x_k)', '||g_k||', 'alpha_k-1', 'beta_k'}, dc];
info.flag = flag;
if flag == 0
  info.message = sprintf('||g|| <= tolg em %d iterações', k);
else
  info.message = sprintf('atingiu o número máximo de iterações (kmax = %d) sem ||g|| <= tolg', kmax);
end

if verbose
  fmt = @(v) strjoin(arrayfun(@(t) sprintf('%9.4f', t), v, 'UniformOutput', false), '; ');
  fprintf('%4s %22s %14s %11s %9s %8s %22s\n', 'k', 'x_k', 'f(x_k)', '||g_k||', 'alpha', 'beta', 'd_k');
  for r = 1:size(hist, 1)
    fprintf('%4d  (%s) %14.6f %11.2e %9.4f %8.4f  (%s)\n', hist(r, 1), fmt(hist(r, 2:n + 1)), ...
            hist(r, n + 2:n + 5), fmt(hist(r, n + 6:end)));
  end
end
end

function v = campo(s, nome, omissao)
% valor do campo nome da estrutura s, ou o valor por omissão
if isfield(s, nome) && ~isempty(s.(nome)), v = s.(nome); else, v = omissao; end
end
