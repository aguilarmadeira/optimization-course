function [x, fx, info] = LocalRandomSearch(f, x0, r0, m, gamma, tolx, kmax, ftarget, s, verbose)
%LOCALRANDOMSEARCH  Pesquisa aleatória localizada: m amostras em x + r[-1,1]^n.
%
%   [x, fx, info] = LocalRandomSearch(f, x0, r0)
%   [x, fx, info] = LocalRandomSearch(f, x0, r0, m, gamma, tolx, kmax, ftarget, s, verbose)
%
%   Em cada iteração gera m pontos x + xi, xi uniforme em [-r,r]^n; se o
%   melhor deles melhora f(x), aceita-o; senão, r <- gamma*r.
%   Para quando r <= tolx ou k >= kmax.
%
%   Entradas
%     f        função (handle) de x, vetor linha com n componentes
%     x0       ponto inicial (vetor linha)
%     r0       raio inicial
%     m        (opcional, 20) amostras por iteração
%     gamma    (opcional, 0.9) fator de redução do raio sem progresso
%     tolx     (opcional, 1e-8) tolerância no raio: para se r <= tolx
%     kmax     (opcional, 5000) número máximo de iterações
%     ftarget  (opcional, -Inf) alvo para info.nhit
%     s        (opcional, []) semente: se não for vazia, faz rng(s)
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        ponto atual (o melhor encontrado)
%     fx       f(x)
%     info     estrutura com
%       .nfev      avaliações de f: 1 + m*nit (inclui f(x0))
%       .nhit      primeira avaliação com f < ftarget, contando f(x0) como a
%                  1.ª e as m amostras de cada iteração pela ordem em que são
%                  avaliadas (Inf se nunca): o custo até atingir o alvo usado
%                  na comparação final do capítulo 3
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       iterações feitas (k)
%       .history   uma linha por k = 0, 1, ..., nit: [k  n_f  r  f(x)  x]
%       .cols      nomes das colunas de .history
%       .flag      0 (r <= tolx) ou 1 (atingiu kmax)
%       .message   mensagem de paragem
%
%   Método estocástico: o resultado é uma variável aleatória. Reporta-se a
%   mediana e os quartis de 20-30 corridas com sementes diferentes.
%
%   Otimização — deck 3.2.1 (Método de pesquisa aleatória, variante localizada).
%   Reproduz os exemplos dos slides: ver ex03_2_1_random_search.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não há limites para as variáveis).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 4 || isempty(m), m = 20; end
if nargin < 5 || isempty(gamma), gamma = 0.9; end
if nargin < 6 || isempty(tolx), tolx = 1e-8; end
if nargin < 7 || isempty(kmax), kmax = 5000; end
if nargin < 8 || isempty(ftarget), ftarget = -Inf; end
if nargin < 9, s = []; end
if nargin < 10, verbose = false; end
if ~isempty(s), rng(s); end      % semente

% --- núcleo (igual ao dos slides, com as contagens) ----------------------
x = x0(:)';  n = numel(x);  fx = f(x);  r = r0;  k = 0;
nfev = 1;
nhit = Inf;  if fx < ftarget, nhit = 1; end
hist = [k, nfev, r, fx, x];
while true
  X = repmat(x, m, 1) + r*(2*rand(m, n) - 1);   % m pontos em x + [-r,r]^n
  v = zeros(m, 1);
  for j = 1:m
    v(j) = f(X(j, :));
    if v(j) < ftarget && isinf(nhit), nhit = nfev + j; end
  end
  nfev = nfev + m;
  [vmin, i] = min(v);
  if vmin < fx                     % o melhor deles melhora f(x): aceitar
    x = X(i, :);  fx = vmin;
  else                             % sem progresso: encolher
    r = gamma*r;
  end
  k = k + 1;
  hist(end + 1, :) = [k, nfev, r, fx, x];
  if r <= tolx || k >= kmax, break; end
end
% -------------------------------------------------------------------------

cols = [{'k', 'n_f', 'r', 'f(x)'}, cell(1, n)];
for j = 1:n, cols{4 + j} = sprintf('x%d', j); end

info.nfev = nfev;
info.nhit = nhit;
info.ngev = 0;
info.nhev = 0;
info.nit = k;
info.history = hist;
info.cols = cols;
if r <= tolx
  info.flag = 0;  info.message = sprintf('r = %.3g <= tolx', r);
else
  info.flag = 1;  info.message = sprintf('atingiu kmax = %d iterações', kmax);
end

if verbose
  fprintf('%5s %7s %10s %12s', cols{1:4});  fprintf(' %10s', cols{5:end});  fprintf('\n');
  for i = 1:size(hist, 1)
    fprintf('%5d %7d %10.3e %12.4e', hist(i, 1:4));  fprintf(' %10.6f', hist(i, 5:end));
    fprintf('\n');
  end
end
end
