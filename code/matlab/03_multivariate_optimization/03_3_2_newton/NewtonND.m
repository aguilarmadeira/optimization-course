function [x, fx, info] = NewtonND(f, grad, hess, x0, tolg, kmax, opts)
%NEWTONND  Método de Newton em R^n, puro ou amortecido: H_k d_k = -g_k.
%
%   [x, fx, info] = NewtonND(f, grad, hess, x0, tolg, kmax)
%   [x, fx, info] = NewtonND(f, grad, hess, x0, tolg, kmax, opts)
%
%   Entradas
%     f       função (handle) de R^n em R ([] só no Newton puro: n_f = 0)
%     grad    gradiente de f (handle que devolve um vetor coluna)
%     hess    Hessiana de f (handle que devolve uma matriz n x n)
%     x0      ponto inicial (vetor coluna)
%     tolg    pára quando ||g_k|| <= tolg (flag 0)
%     kmax    número máximo de iterações (flag 1, com aviso)
%     opts    (opcional) estrutura com os campos (todos opcionais)
%       .damped   false (por omissão): Newton puro, passo completo alpha = 1;
%                 true: amortecido -- se g'*d >= 0 usa d = -g (salvaguarda)
%                 e alpha por pesquisa em linha
%       .modify   (só amortecido; false) se true e H_k não for definida
%                 positiva ([~, p] = chol(H), p > 0), usa H_k + mu*I com
%                 mu = 1e-3*max(1, norm(H_k,'fro')), 10 mu, 100 mu, ... até
%                 ser definida positiva
%       .ls       pesquisa em linha do amortecido: 'brent' (por omissão;
%                 tol 1e-10, a das figuras) ou 'fminbnd' -- ver LineSearch;
%                 ou 'armijo': alpha = 1, 1/2, 1/4, ... até
%                 f(x + alpha*d) <= f(x) + c1*alpha*g'*d, c1 = 1e-4
%                 (o valor aceite é reaproveitado como f no novo iterando)
%       .lstol    tolerância da pesquisa em linha ([] = a de LineSearch)
%       .verbose  (false) se true, imprime a tabela das iterações
%
%   Saídas
%     x      último iterando x_k
%     fx     f(x_k) (NaN se f = [])
%     info   estrutura com
%       .nfev     avaliações de f: f em cada iterando (x_0 incluído) + as
%                 avaliações das pesquisas em linha (amortecido)
%       .nfev_ls  avaliações de f só nas pesquisas em linha
%       .ngev     gradientes: nit + 1 (o último é o do teste de paragem)
%       .nhev     Hessianas: nit
%       .nit      número de iterações
%       .history  uma linha por k = 0, 1, ..., nit:
%                 [k  x_k'  f(x_k)  ||g_k||  alpha_{k-1}  lambda_min(H_{k-1})  dir_{k-1}]
%                 dir: 0 Newton, 1 H + mu*I, 2 -g (salvaguarda); NaN se não se
%                 aplica. lambda_min é só para leitura (não é usado).
%       .cols     nomes das colunas de .history
%       .flag     0 (||g|| <= tolg) ou 1 (atingiu kmax)
%       .message  mensagem de paragem
%
%   Newton procura g = 0, não mínimos: o puro pode ir para máximos ou selas
%   se H_k não for definida positiva (ver o exemplo de Himmelblau).
%   O amortecido precisa de LineSearch.m (code/matlab/common/). Os exemplos
%   tratam do caminho; fora deles, chamar antes
%   addpath('<repo>/code/matlab'); uc_setup
%
%   Otimização — deck 3.3.2 (Método de Newton).
%   Reproduz os exemplos dos slides: ver ex03_3_2_newton.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. a condição de Wolfe). Por omissão, o passo
%   amortecido é o minimizante da linha; com opts.ls = 'armijo', recuo
%   alpha = 1, 1/2, 1/4, ... até ao decréscimo suficiente de Armijo.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 7, opts = struct(); end
damped = campo(opts, 'damped', false);
modify = campo(opts, 'modify', false);
ls = campo(opts, 'ls', 'brent');
lstol = campo(opts, 'lstol', []);
verbose = campo(opts, 'verbose', false);
temf = ~isempty(f);
if damped && ~temf, error('NewtonND: o Newton amortecido precisa de f'); end

x = x0(:);
n = numel(x);
nfev = 0;  ngev = 0;  nhev = 0;  nfev_ls = 0;
if temf, fx = f(x); nfev = nfev + 1; else, fx = NaN; end   % f em x_0
hist = [0, x', fx, NaN, NaN, NaN, NaN];
flag = 1;
k = 0;
% --- núcleo (igual ao dos slides, com as contagens) ----------------------
while true
  if k >= kmax, flag = 1; break, end
  g = grad(x);  g = g(:);  ngev = ngev + 1;
  hist(k + 1, n + 3) = norm(g);
  if norm(g) <= tolg, flag = 0; break, end
  H = hess(x);  nhev = nhev + 1;
  lmin = min(eig((H + H')/2));                 % só para leitura
  tipo = 0;
  if damped && modify
    [~, p] = chol(H);                          % p = 0 se H for definida positiva
    if p > 0
      mu = 1e-3*max(1, norm(H, 'fro'));
      while p > 0
        H = H + mu*eye(n);  mu = 10*mu;
        [~, p] = chol(H);
      end
      tipo = 1;
    end
  end
  d = -(H \ g);                                % H*d = -g (sem inverter H)
  if damped
    if g'*d >= 0
      d = -g;  tipo = 2;                       % salvaguarda
    end
    if strcmpi(ls, 'armijo')                   % recuo até ao decréscimo suficiente
      gd = g'*d;  alpha = 1;  ft = f(x + d);  nls = 1;
      while ft > fx + 1e-4*alpha*gd && nls < 60
        alpha = alpha/2;  ft = f(x + alpha*d);  nls = nls + 1;
      end
    else
      phi = @(a) f(x + a*d);                   % pesquisa em linha
      [alpha, nls] = LineSearch(phi, ls, lstol);
    end
    nfev = nfev + nls;  nfev_ls = nfev_ls + nls;
  else
    alpha = 1;                                 % passo completo
  end
  x = x + alpha*d;
  if damped && strcmpi(ls, 'armijo')
    fx = ft;                                   % já avaliado no recuo
  elseif temf
    fx = f(x); nfev = nfev + 1;                % f no novo iterando
  end
  k = k + 1;
  hist(k + 1, :) = [k, x', fx, NaN, alpha, lmin, tipo];
end
% -------------------------------------------------------------------------

info.nfev = nfev;
info.nfev_ls = nfev_ls;
info.ngev = ngev;
info.nhev = nhev;
info.nit = k;
info.history = hist;
xc = cell(1, n);
for i = 1:n, xc{i} = sprintf('x%d', i); end
info.cols = [{'k'}, xc, {'f(x_k)', '||g_k||', 'alpha_k-1', 'lmin(H_k-1)', 'direcao_k-1'}];
info.flag = flag;
if flag == 0
  info.message = sprintf('||g|| <= tolg em %d iterações', k);
else
  info.message = sprintf('atingiu o número máximo de iterações (kmax = %d) sem ||g|| <= tolg', kmax);
  warning('NewtonND: %s', info.message);
end

if verbose
  nomes = {'Newton', 'Newton com H + mu I', '-g (salvaguarda)'};
  fprintf('%4s %24s %16s %11s %9s  %s\n', 'k', 'x_k', 'f(x_k)', '||g_k||', 'alpha', 'direção');
  for r = 1:size(hist, 1)
    if isnan(hist(r, n + 6)), s = ''; else, s = nomes{hist(r, n + 6) + 1}; end
    fprintf('%4d  (%s) %16.6e %11.2e %9.4f  %s\n', hist(r, 1), ...
            strjoin(arrayfun(@(v) sprintf('%9.4f', v), hist(r, 2:n + 1), 'UniformOutput', false), '; '), ...
            hist(r, n + 2:n + 4), s);
  end
end
end

function v = campo(s, nome, omissao)
% valor do campo nome da estrutura s, ou o valor por omissão
if isfield(s, nome) && ~isempty(s.(nome)), v = s.(nome); else, v = omissao; end
end
