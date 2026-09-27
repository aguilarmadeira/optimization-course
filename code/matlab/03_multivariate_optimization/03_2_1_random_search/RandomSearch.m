function [x, fx, info] = RandomSearch(f, a, b, N, s, ftarget, verbose)
%RANDOMSEARCH  Pesquisa aleatória pura: min f(x) em D = [a,b], com N amostras.
%
%   [x, fx, info] = RandomSearch(f, a, b, N)
%   [x, fx, info] = RandomSearch(f, a, b, N, s, ftarget, verbose)
%
%   Entradas
%     f        função (handle) de x, vetor linha com n componentes
%     a, b     limites inferior e superior de D (vetores linha, a < b)
%     N        número de amostras uniformes em D (n_f = N)
%     s        (opcional, []) semente: se não for vazia, faz rng(s) antes
%              de gerar os pontos; se for vazia, usa o estado atual do gerador
%     ftarget  (opcional, -Inf) alvo para info.nhit
%     verbose  (opcional, false) se true, imprime as melhorias do melhor ponto
%
%   Saídas
%     x        melhor ponto encontrado (xbest)
%     fx       f(xbest)
%     info     estrutura com
%       .nfev      avaliações de f: N
%       .nhit      primeira avaliação com f < ftarget (Inf se nunca)
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       número de amostras geradas (= N)
%       .history   uma linha por melhoria do melhor ponto:
%                  [i  xbest  fbest]  (o melhor até agora nunca piora)
%       .cols      nomes das colunas de .history
%       .flag      0 (gerou as N amostras)
%       .message   mensagem de paragem
%
%   Método estocástico: o resultado é uma variável aleatória. Reporta-se a
%   mediana e os quartis de 20-30 corridas com sementes diferentes.
%
%   Otimização — deck 3.2.1 (Método de pesquisa aleatória).
%   Reproduz os exemplos dos slides: ver ex03_2_1_random_search.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não trata valores de f que sejam NaN).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 5, s = []; end
if nargin < 6 || isempty(ftarget), ftarget = -Inf; end
if nargin < 7, verbose = false; end
a = a(:)';  b = b(:)';  n = numel(a);
if any(~(a < b)), error('RandomSearch: é preciso a < b.'); end

nfev = 0;  nhit = Inf;  hist = zeros(0, n + 2);

% --- núcleo (igual ao dos slides, com as contagens) ----------------------
if ~isempty(s), rng(s); end      % semente
fbest = inf;
for i = 1:N
    x  = a + (b-a).*rand(1,n);   % em [a,b]
    fx = f(x);
    nfev = nfev + 1;
    if fx < ftarget && isinf(nhit), nhit = nfev; end
    if fx < fbest
        xbest = x;
        fbest = fx;
        hist(end + 1, :) = [i, xbest, fbest];
    end
end
% -------------------------------------------------------------------------

x = xbest;  fx = fbest;
cols = cell(1, n + 2);  cols{1} = 'i';  cols{end} = 'f_best';
for j = 1:n, cols{j + 1} = sprintf('x%d', j); end

info.nfev = nfev;
info.nhit = nhit;
info.ngev = 0;
info.nhev = 0;
info.nit = N;
info.history = hist;
info.cols = cols;
info.flag = 0;
info.message = sprintf('gerou as %d amostras pedidas', N);

if verbose
  fprintf('%6s', cols{1});  fprintf(' %10s', cols{2:end});  fprintf('\n');
  for i = 1:size(hist, 1)
    fprintf('%6d', hist(i, 1));  fprintf(' %10.4f', hist(i, 2:end));  fprintf('\n');
  end
end
end
