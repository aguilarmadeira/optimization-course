function [x, fx, info] = EpsilonConstraint(f1, f2, x0, lb, ub, E, verbose, grafico)
%EPSILONCONSTRAINT  Método do epsilon-constrangimento (dois objetivos).
%
%   [x, fx, info] = EpsilonConstraint(f1, f2, x0, lb, ub, E)
%   [x, fx, info] = EpsilonConstraint(f1, f2, x0, lb, ub, E, verbose, grafico)
%
%   Para cada valor e de E resolve o subproblema do cap. 4
%
%        min f2(x)   s.a.   g(x) = e - f1(x) >= 0,   lb <= x <= ub
%
%   com arranque a quente (a solução para um e é o ponto inicial do seguinte)
%   e guarda o multiplicador u de f1 <= e (a taxa de troca, u = -df2*/de).
%
%   Entradas
%     f1, f2   funções (handles) de x, com valores escalares; f1 é o objetivo
%              limitado, f2 o objetivo minimizado (como no slide)
%     x0       ponto inicial do primeiro subproblema
%     lb, ub   limites das variáveis
%     E        valores de e a varrer (p. ex. 0.05:0.1:0.95); e tem as
%              unidades de f1 (não é uma tolerância numérica)
%     verbose  (opcional, false) se true, imprime a tabela do varrimento
%     grafico  (opcional, false) se true, desenha os pontos obtidos
%
%   Saídas
%     x        um ponto por linha (numel(E) x n)
%     fx       [f1(x) f2(x)] de cada ponto (numel(E) x 2)
%     info     estrutura com
%       .u         multiplicadores de f1 <= e (um por e)
%       .nfev      n_f do varrimento: soma de out.funcCount (fmincon, como no
%                  slide) ou das avaliações de f2 (sqp); inclui as dos
%                  gradientes por diferenças finitas do solver
%       .nf2, .nf1 chamadas a f2 e a f1 contadas diretamente
%       .ncalls    nf1 + nf2 (chamadas escalares)
%       .ngev, .nhev   0 (os gradientes são numéricos, dentro do solver)
%       .nit       soma das iterações do solver
%       .history   uma linha por e: [e  x  f1  f2  u  nf2(e)  nf1(e)]
%       .cols      nomes das colunas de .history
%       .exitflag  código de saída do solver em cada subproblema
%       .solver    'fmincon' (MATLAB) ou 'sqp' (GNU Octave)
%       .flag      0 (todos os subproblemas terminaram bem) ou 1
%       .message   mensagem de paragem
%
%   Solver: fmincon (Optimization Toolbox), como no slide, se existir;
%   senão (GNU Octave) sqp, com a mesma formulação g = e - f1 >= 0. Os pontos
%   e os multiplicadores coincidem; as contagens dependem do solver (as do
%   slide, 60/90/150, são do SLSQP do SciPy e só se verificam em Python).
%
%   Otimização — deck 6.4 (Método do epsilon-constrangimento).
%   Reproduz o exemplo dos slides: ver ex06_4_epsilon_constraint.m
%   Implementação didática — resolve cada subproblema localmente: num problema
%   não convexo as garantias do slide exigem subproblemas resolvidos
%   globalmente (p. ex. com vários pontos iniciais).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 7, verbose = false; end
if nargin < 8, grafico = false; end

if exist('fmincon', 'file') == 2
  solver = 'fmincon';
elseif exist('sqp', 'file') == 2
  solver = 'sqp';
else
  error('EpsilonConstraint: é preciso fmincon (MATLAB) ou sqp (GNU Octave).');
end

E = E(:)';  ne = numel(E);  n = numel(lb);
x0 = x0(:)';  lb = lb(:)';  ub = ub(:)';
x = zeros(ne, n);  fx = zeros(ne, 2);
hist = zeros(ne, n + 6);
info.u = zeros(ne, 1);
info.exitflag = zeros(ne, 1);
info.nit = 0;

% contadores diretos das chamadas a f1 e f2
f1c = @(z) contar(f1, z, 1);
f2c = @(z) contar(f2, z, 2);
contar('reset');

% --- núcleo (o do slide, com as contagens) --------------------------------
info.nfev = 0;
for k = 1:ne
  e = E(k);
  antes = contar('get');
  if strcmp(solver, 'fmincon')
    opt = optimset('Display', 'off');
    nonl = @(z) deal(f1c(z) - e, []);          % f1 <= e  (c = f1 - e <= 0)
    [xk, ~, ef, out, lambda] = fmincon(f2c, x0, ...
                  [], [], [], [], lb, ub, nonl, opt);
    u = lambda.ineqnonlin;                     % taxa de troca
    info.nfev = info.nfev + out.funcCount;
    it = out.iterations;
  else                                         % GNU Octave
    g = @(z) e - f1c(z);                       % forma da UC: g >= 0
    [xk, ~, ef, it, ~, lambda] = sqp(x0(:), f2c, [], g, lb(:), ub(:));
    u = lambda(1);                             % 1.º: a restrição g; depois os limites
    xk = xk(:)';
  end
  x0 = xk;                                     % a quente
  % ------------------------------------------------------------------------
  depois = contar('get');
  nfk = depois - antes;                        % [nf1(e) nf2(e)]
  if strcmp(solver, 'sqp')
    info.nfev = info.nfev + nfk(2);            % avaliações de f2
  end
  x(k, :) = xk;
  fx(k, :) = [f1(xk), f2(xk)];                 % valores (não contam)
  info.u(k) = u;
  info.exitflag(k) = ef;
  info.nit = info.nit + it;
  hist(k, :) = [e, xk, fx(k, :), u, nfk(2), nfk(1)];
end

nf = contar('get');
info.nf1 = nf(1);
info.nf2 = nf(2);
info.ncalls = nf(1) + nf(2);
info.ngev = 0;
info.nhev = 0;
info.history = hist;
xs = cell(1, n);
for i = 1:n, xs{i} = sprintf('x%d', i); end
info.cols = [{'e'}, xs, {'f1', 'f2', 'u', 'nf2', 'nf1'}];
info.solver = solver;
if strcmp(solver, 'fmincon')
  falhou = any(info.exitflag <= 0);
else
  falhou = any(~ismember(info.exitflag, [101 104]));   % 101/104: paragem normal
end
info.flag = double(falhou);
if falhou
  info.message = sprintf('%s: algum subproblema não terminou bem (ver info.exitflag)', solver);
else
  info.message = sprintf('%s: %d subproblemas resolvidos', solver, ne);
end

if verbose
  fprintf('%6s', 'e');
  for i = 1:n, fprintf(' %8s', xs{i}); end
  fprintf(' %8s %8s %8s %5s %5s\n', 'f1', 'f2', 'u', 'nf2', 'nf1');
  for k = 1:ne
    fprintf('%6.2f', hist(k, 1));
    fprintf(' %8.4f', hist(k, 2:n+1));
    fprintf(' %8.4f %8.4f %8.4f %5d %5d\n', hist(k, n+2:n+4), hist(k, n+5:n+6));
  end
end

if grafico
  plot(fx(:, 1), fx(:, 2), 'o');
  xlabel('f_1');  ylabel('f_2');
  title('\epsilon-constrangimento: um ponto por \epsilon');
end
end


function v = contar(fun, z, i)
%CONTAR  Conta as chamadas a f1 (i = 1) e f2 (i = 2).
%   contar('reset') põe os contadores a zero; contar('get') devolve [nf1 nf2].
persistent nc
if isempty(nc), nc = [0 0]; end
if ischar(fun)
  if strcmp(fun, 'reset'), nc = [0 0]; end
  v = nc;
  return
end
nc(i) = nc(i) + 1;
v = fun(z);
end
