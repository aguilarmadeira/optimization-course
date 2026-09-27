function [x, fx, info] = WeightedTchebycheff(F, X0, lb, ub, g, W, zid, verbose)
%WEIGHTEDTCHEBYCHEFF  Tchebycheff ponderado: min max_i w_i (f_i(x) - z_i^id).
%
%   [x, fx, info] = WeightedTchebycheff(F, X0, lb, ub, g, W, zid)
%   [x, fx, info] = WeightedTchebycheff(F, X0, lb, ub, g, W, zid, verbose)
%
%   O max não é diferenciável; como na nota do slide («Sem o max (4.2)»),
%   resolve-se a reformulação com mais uma variável theta:
%
%        min theta   s.a.   theta - w_i (f_i(x) - z_i^id) >= 0,  i = 1..m,
%                           g(x) >= 0,   lb <= x <= ub,
%
%   com fmincon ou, no GNU Octave, com sqp (variáveis z = [x theta]).
%
%   Entradas: como em WeightedSum, mais
%     zid      ponto ideal (vetor com m componentes)
%     X0       um ponto por linha; um só ponto => arranque a quente;
%              theta0 = max_i w_i (f_i(x0) - z_i^id)
%
%   Saídas
%     x, fx    um ponto por linha e F(x) de cada um
%     info     .nfev (chamadas a F: aqui F aparece nas restrições), .ngev,
%              .nhev (0), .nit, .history = [w  x  F(x)  theta  n_f(w)],
%              .cols, .exitflag, .solver, .flag, .message
%
%   Otimização — deck 6.3 (Métodos de agregação de objetivos).
%   Reproduz o exemplo dos slides: ver ex06_3_aggregation.m
%   Implementação didática — cada problema é resolvido localmente; não
%   inclui a versão aumentada (+ rho sum_i (f_i - z_i^id)).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 8, verbose = false; end
if exist('fmincon', 'file') == 2
  solver = 'fmincon';
elseif exist('sqp', 'file') == 2
  solver = 'sqp';
else
  error('WeightedTchebycheff: é preciso fmincon (MATLAB) ou sqp (GNU Octave).');
end

n = numel(lb);  lb = lb(:)';  ub = ub(:)';  zid = zid(:)';
X0 = reshape(X0', n, [])';
[nw, m] = size(W);
x = zeros(nw, n);  fx = zeros(nw, m);
hist = zeros(nw, m + n + m + 2);
info.exitflag = zeros(nw, 1);
info.nit = 0;
Fc = @(z) contar(F, z);
contar('reset');

x0 = X0(1, :);
for k = 1:nw
  w = W(k, :);
  % restrições em z = [x theta], na forma da UC (>= 0)
  if isempty(g)
    gz = @(z) z(n+1) - (w .* (Fc(z(1:n)) - zid))';
  else
    gz = @(z) [z(n+1) - (w .* (Fc(z(1:n)) - zid))'; reshape(g(z(1:n)), [], 1)];
  end
  antes = contar('get');
  melhor = Inf;
  for s = 1:size(X0, 1)
    if size(X0, 1) > 1, x0 = X0(s, :); end
    z0 = [x0, max(w .* (F(x0) - zid))];
    if strcmp(solver, 'fmincon')
      opt = optimset('Display', 'off');
      nonl = @(z) deal(-gz(z), []);
      [zs, ts, ef, out] = fmincon(@(z) z(n+1), z0, [], [], [], [], ...
                                  [lb -Inf], [ub Inf], nonl, opt);
      it = out.iterations;
    else
      ws = warning('off', 'all');            % ver a nota em WeightedSum
      try
        [zs, ts, ef, it] = sqp(z0(:), @(z) z(n+1), [], gz, [lb -Inf]', [ub Inf]');
        zs = zs(:)';
      catch
        zs = z0;  ts = Inf;  ef = -1;  it = 0;
      end
      warning(ws);
    end
    info.nit = info.nit + it;
    if ts < melhor
      melhor = ts;  zk = zs;  info.exitflag(k) = ef;
    end
  end
  x0 = zk(1:n);                            % a quente (só com um ponto inicial)
  x(k, :) = zk(1:n);
  fx(k, :) = F(zk(1:n));
  hist(k, :) = [w, zk(1:n), fx(k, :), zk(n+1), contar('get') - antes];
end

info.nfev = contar('get');
info.ngev = 0;
info.nhev = 0;
info.history = hist;
ws = cell(1, m);  xs = cell(1, n);  fs = cell(1, m);
for i = 1:m, ws{i} = sprintf('w%d', i);  fs{i} = sprintf('f%d', i); end
for i = 1:n, xs{i} = sprintf('x%d', i); end
info.cols = [ws, xs, fs, {'theta', 'nf'}];
info.solver = solver;
if strcmp(solver, 'fmincon')
  falhou = any(info.exitflag <= 0);
else
  falhou = any(~ismember(info.exitflag, [101 102 104]));
end
info.flag = double(falhou);
if falhou
  info.message = sprintf('%s: algum problema não terminou bem (ver info.exitflag)', solver);
else
  info.message = sprintf('%s: %d problemas escalares resolvidos', solver, nw);
end

if verbose
  fprintf('%s\n', sprintf(' %7s', info.cols{:}));
  for k = 1:nw
    fprintf(' %7.3f', hist(k, 1:m));
    fprintf(' %7.4f', hist(k, m+1:end-1));
    fprintf(' %7d\n', hist(k, end));
  end
end
end


function v = contar(fun, z)
%CONTAR  Conta as chamadas a F.  contar('reset') / contar('get').
persistent nc
if isempty(nc), nc = 0; end
if ischar(fun)
  if strcmp(fun, 'reset'), nc = 0; end
  v = nc;
  return
end
nc = nc + 1;
v = fun(z);
end
