function [x, fx, info] = WeightedSum(F, X0, lb, ub, g, W, verbose)
%WEIGHTEDSUM  Soma ponderada: min sum_i w_i f_i(x), um problema por vetor de pesos.
%
%   [x, fx, info] = WeightedSum(F, X0, lb, ub, g, W)
%   [x, fx, info] = WeightedSum(F, X0, lb, ub, g, W, verbose)
%
%   Para cada linha w de W resolve
%
%        min  w * F(x)'   s.a.   g(x) >= 0,   lb <= x <= ub
%
%   com fmincon (como no slide do 6.4: fmincon(@(x) w1*f1(x)+w2*f2(x), ...)
%   sem nonl) ou, no GNU Octave, com sqp.
%
%   Entradas
%     F        função (handle) que devolve o vetor linha [f1(x) ... fm(x)]
%     X0       pontos iniciais, um por linha (numel(lb) colunas):
%              - um só ponto: arranque a quente (a solução para um w é o
%                ponto inicial do w seguinte), como no slide;
%              - vários pontos: para cada w, parte de todos e fica com o de
%                menor soma ponderada (sem arranque a quente); útil em
%                problemas não convexos, onde o solver é local
%     lb, ub   limites das variáveis
%     g        restrições na forma da UC, g(x) >= 0 (handle), ou [] se não há
%     W        pesos, um vetor por linha (nw x m), w_i >= 0, sum_i w_i = 1
%     verbose  (opcional, false) se true, imprime a tabela dos pesos
%
%   Saídas
%     x        um ponto por linha (nw x n)
%     fx       F(x) de cada ponto (nw x m)
%     info     estrutura com
%       .nfev      chamadas a F (= avaliações da soma ponderada; incluem as
%                  dos gradientes por diferenças finitas do solver; com
%                  fmincon, é a soma de out.funcCount)
%       .ngev, .nhev   0 (os gradientes são numéricos, dentro do solver)
%       .nit       soma das iterações do solver
%       .history   uma linha por w: [w  x  F(x)  n_f(w)]
%       .cols      nomes das colunas de .history
%       .exitflag  código de saída do melhor arranque, por w
%       .solver    'fmincon' (MATLAB) ou 'sqp' (GNU Octave)
%       .flag      0 (todos terminaram bem) ou 1;  .message
%
%   Otimização — deck 6.3 (Métodos de agregação de objetivos).
%   Reproduz o exemplo dos slides: ver ex06_3_aggregation.m
%   Implementação didática — cada problema escalar é resolvido localmente; só
%   um minimizante global da soma ponderada tem garantia de ser eficiente.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 7, verbose = false; end
solver = escolher_solver();

n = numel(lb);  lb = lb(:)';  ub = ub(:)';
X0 = reshape(X0', n, [])';                 % um ponto por linha
[nw, m] = size(W);
x = zeros(nw, n);  fx = zeros(nw, m);
hist = zeros(nw, m + n + m + 1);
info.exitflag = zeros(nw, 1);
info.nit = 0;
Fc = @(z) contar(F, z);
contar('reset');

if ~isempty(g), g = @(z) reshape(g(z), [], 1); end   % restrições em coluna
x0 = X0(1, :);
for k = 1:nw
  w = W(k, :);
  fobj = @(z) w * Fc(z)';                  % soma ponderada
  antes = contar('get');
  melhor = Inf;
  for s = 1:size(X0, 1)
    if size(X0, 1) > 1, x0 = X0(s, :); end
    [xs, fs, ef, it] = resolver(solver, fobj, x0, lb, ub, g);
    info.nit = info.nit + it;
    if fs < melhor
      melhor = fs;  xk = xs;  info.exitflag(k) = ef;
    end
  end
  x0 = xk;                                 % a quente (só com um ponto inicial)
  x(k, :) = xk;
  fx(k, :) = F(xk);                        % valores (não contam)
  hist(k, :) = [w, xk, fx(k, :), contar('get') - antes];
end

info.nfev = contar('get');
info.ngev = 0;
info.nhev = 0;
info.history = hist;
ws = cell(1, m);  xs = cell(1, n);  fs = cell(1, m);
for i = 1:m, ws{i} = sprintf('w%d', i);  fs{i} = sprintf('f%d', i); end
for i = 1:n, xs{i} = sprintf('x%d', i); end
info.cols = [ws, xs, fs, {'nf'}];
info.solver = solver;
[info.flag, info.message] = estado(solver, info.exitflag, nw);

if verbose
  fprintf('%s\n', sprintf(' %7s', info.cols{:}));
  for k = 1:nw
    fprintf(' %7.2f', hist(k, 1:m));
    fprintf(' %7.4f', hist(k, m+1:end-1));
    fprintf(' %7d\n', hist(k, end));
  end
end
end


function solver = escolher_solver()
if exist('fmincon', 'file') == 2
  solver = 'fmincon';
elseif exist('sqp', 'file') == 2
  solver = 'sqp';
else
  error('é preciso fmincon (MATLAB) ou sqp (GNU Octave).');
end
end


function [x, fval, ef, it] = resolver(solver, fobj, x0, lb, ub, g)
%RESOLVER  min fobj(x) s.a. g(x) >= 0, lb <= x <= ub (fmincon ou sqp).
if strcmp(solver, 'fmincon')
  opt = optimset('Display', 'off');
  if isempty(g)
    nonl = [];
  else
    nonl = @(z) deal(-g(z), []);            % fmincon: c = -g <= 0
  end
  [x, fval, ef, out] = fmincon(fobj, x0, [], [], [], [], lb, ub, nonl, opt);
  it = out.iterations;
else
  % o sqp do Octave pode parar com erro quando o subproblema QP é ilimitado
  % (objetivo côncavo); esse arranque conta como falhado (ef = -1)
  ws = warning('off', 'all');
  try
    [x, fval, ef, it] = sqp(x0(:), fobj, [], g, lb(:), ub(:));
    x = x(:)';
  catch
    x = x0;  fval = Inf;  ef = -1;  it = 0;
  end
  warning(ws);
end
end


function [flag, msg] = estado(solver, ef, nw)
if strcmp(solver, 'fmincon')
  falhou = any(ef <= 0);
else
  falhou = any(~ismember(ef, [101 102 104]));   % paragem normal do sqp
end
flag = double(falhou);
if falhou
  msg = sprintf('%s: algum problema não terminou bem (ver info.exitflag)', solver);
else
  msg = sprintf('%s: %d problemas escalares resolvidos', solver, nw);
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
