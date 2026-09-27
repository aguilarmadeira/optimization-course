function [x, fx, info] = SteepestDescent(f, grad, x0, tolg, kmax, opts)
%STEEPESTDESCENT  Método do gradiente (steepest descent, Cauchy) com pesquisa em linha.
%
%   [x, fx, info] = SteepestDescent(f, grad, x0, tolg, kmax)
%   [x, fx, info] = SteepestDescent(f, grad, x0, tolg, kmax, opts)
%
%   x_{k+1} = x_k - alpha_k g_k, com alpha_k da pesquisa em linha (cap. 2).
%
%   Entradas
%     f      função (handle) de R^n em R
%     grad   gradiente de f (handle que devolve um vetor); [] usa GradFD
%            (diferenças centrais: 2n avaliações de f por gradiente, em nfev)
%     x0     ponto inicial (vetor coluna)
%     tolg   pára quando ||g_k|| <= tolg (flag 0)
%     kmax   número máximo de iterações (flag 1; o gradiente em x_kmax não é
%            calculado, por isso n_g = kmax)
%     opts   (opcional) estrutura com os campos (todos opcionais)
%       .tolx     pára quando ||x_{k+1}-x_k||/max(1,||x_k||) <= tolx (flag 2);
%                 por omissão 0 (só pára se o ponto não mudar)
%       .ls       pesquisa em linha: 'brent' (por omissão; tol 1e-10, a das
%                 figuras) ou 'fminbnd' (TolX 1e-4) -- ver LineSearch
%       .lstol    tolerância da pesquisa em linha ([] = a de LineSearch)
%       .verbose  (false) se true, imprime a tabela das iterações
%
%   Saídas
%     x      último iterando x_k
%     fx     f(x_k)
%     info   estrutura com
%       .nfev     avaliações de f: pesquisas em linha (cada uma volta a
%                 avaliar phi(0) = f(x_k) no modo 'brent') + f em cada
%                 iterando (x_0 incluído) (+ 2n por gradiente numérico)
%       .nfev_ls  avaliações de f só nas pesquisas em linha
%       .ngev     gradientes: nit + 1 (0 se numérico)
%       .nhev     0 (o método não usa a Hessiana)
%       .nit      número de iterações (passos)
%       .history  uma linha por k = 0, 1, ..., nit:
%                 [k  x_k'  f(x_k)  ||g_k||  alpha_{k-1}]  (NaN se não calculado)
%       .cols     nomes das colunas de .history
%       .flag     0 (||g|| <= tolg), 1 (atingiu kmax), 2 (passo <= tolx)
%       .message  mensagem de paragem
%
%   Precisa de LineSearch.m (code/matlab/common/) e, se grad = [], de GradFD.m
%   (nesta pasta). Os exemplos tratam do caminho; fora deles, chamar antes
%   addpath('<repo>/code/matlab'); uc_setup
%
%   Otimização — deck 3.3.1 (Método do gradiente).
%   Reproduz os exemplos dos slides: ver ex03_3_1_steepest_descent.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. pesquisa em linha inexata com condições de Wolfe).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 6, opts = struct(); end
tolx = campo(opts, 'tolx', 0);
ls = campo(opts, 'ls', 'brent');
lstol = campo(opts, 'lstol', []);
verbose = campo(opts, 'verbose', false);

x = x0(:);
n = numel(x);
nfev = 0;  ngev = 0;  nfev_ls = 0;
fx = f(x);  nfev = nfev + 1;             % f no iterando x_0
hist = [0, x', fx, NaN, NaN];
flag = 1;
k = 0;
% --- núcleo (igual ao dos slides, com as contagens) ----------------------
while true
  if k >= kmax, flag = 1; break, end
  if isempty(grad)
    [g, nf] = GradFD(f, x);  nfev = nfev + nf;   % 2n avaliações de f
  else
    g = grad(x);  ngev = ngev + 1;               % n_g = 1
  end
  g = g(:);
  hist(k + 1, n + 3) = norm(g);
  if norm(g) <= tolg, flag = 0; break, end
  d = -g;                                        % direção
  phi = @(a) f(x + a*d);                         % pesquisa em linha
  [alpha, nls] = LineSearch(phi, ls, lstol);
  nfev = nfev + nls;  nfev_ls = nfev_ls + nls;
  xn = x + alpha*d;
  fx = f(xn);  nfev = nfev + 1;                  % f no novo iterando
  passo = norm(xn - x)/max(1, norm(x));
  x = xn;  k = k + 1;
  hist(k + 1, :) = [k, x', fx, NaN, alpha];
  if passo <= tolx, flag = 2; break, end
end
% -------------------------------------------------------------------------

info.nfev = nfev;
info.nfev_ls = nfev_ls;
info.ngev = ngev;
info.nhev = 0;
info.nit = k;
info.history = hist;
xc = cell(1, n);
for i = 1:n, xc{i} = sprintf('x%d', i); end
info.cols = [{'k'}, xc, {'f(x_k)', '||g_k||', 'alpha_k-1'}];
info.flag = flag;
switch flag
  case 0, info.message = sprintf('||g|| <= tolg em %d iterações', k);
  case 1, info.message = sprintf('atingiu o número máximo de iterações (kmax = %d) sem ||g|| <= tolg', kmax);
  otherwise, info.message = sprintf('passo relativo <= tolx em %d iterações', k);
end

if verbose
  fprintf('%4s %22s %14s %12s %9s\n', 'k', 'x_k', 'f(x_k)', '||g_k||', 'alpha');
  for r = 1:size(hist, 1)
    fprintf('%4d  (%s) %14.5f %12.4e %9.4f\n', hist(r, 1), ...
            strjoin(arrayfun(@(v) sprintf('%9.4f', v), hist(r, 2:n + 1), 'UniformOutput', false), '; '), ...
            hist(r, n + 2:n + 4));
  end
end
end

function v = campo(s, nome, omissao)
% valor do campo nome da estrutura s, ou o valor por omissão
if isfield(s, nome) && ~isempty(s.(nome)), v = s.(nome); else, v = omissao; end
end
