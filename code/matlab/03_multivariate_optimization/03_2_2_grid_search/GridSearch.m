function [x, fx, info] = GridSearch(f, a, b, m, verbose)
%GRIDSEARCH  Pesquisa em grelha: min f(x1,x2) em [a1,b1] x [a2,b2], m valores por variável.
%
%   [x, fx, info] = GridSearch(f, a, b, m)
%   [x, fx, info] = GridSearch(f, a, b, m, verbose)
%
%   Avalia f em todas as m^2 combinações de linspace(a1,b1,m) e
%   linspace(a2,b2,m) (a grelha do deck, construída com meshgrid) e devolve
%   o melhor nó. Como no deck, é para n = 2.
%
%   Entradas
%     f        função (handle) de x = [x1 x2]
%     a, b     limites inferior e superior: a = [a1 a2], b = [b1 b2]
%     m        número de valores por variável (m^2 avaliações)
%     verbose  (opcional, false) se true, imprime os 5 melhores nós
%
%   Saídas
%     x        melhor nó da grelha
%     fx       f(x)
%     info     estrutura com
%       .nfev      avaliações de f: m^2
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       0 (não há iterações: uma só passagem pela grelha)
%       .h         espaçamento da grelha, (b - a)/(m - 1)
%       .X1, .X2, .F   a grelha e os valores de f (matrizes m x m do meshgrid)
%       .history   todos os nós, do melhor para o pior: [x1  x2  f]
%                  (as primeiras linhas são os pontos iniciais naturais
%                  para um método local: «global grosseiro -> local fino»)
%       .cols      nomes das colunas de .history
%       .flag      0 (avaliou a grelha toda)
%       .message   mensagem de paragem
%
%   Otimização — deck 3.2.2 (Pesquisa em grelha).
%   Reproduz os exemplos dos slides: ver ex03_2_2_grid_search.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. só serve para n = 2, como no deck).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 5, verbose = false; end
if numel(a) ~= 2 || numel(b) ~= 2, error('GridSearch: esta versão é para n = 2.'); end
if any(~(a < b)) || m < 2, error('GridSearch: é preciso a < b e m >= 2.'); end
a1 = a(1);  a2 = a(2);  b1 = b(1);  b2 = b(2);

% --- núcleo (igual ao dos slides) ----------------------------------------
[X1, X2] = meshgrid(linspace(a1,b1,m), ...
                    linspace(a2,b2,m));
F = arrayfun(@(u,v) f([u v]), X1, X2);
[fmin, i] = min(F(:));  x = [X1(i) X2(i)];
% -------------------------------------------------------------------------

fx = fmin;
[Fs, j] = sort(F(:));

info.nfev = numel(F);           % m^2
info.ngev = 0;
info.nhev = 0;
info.nit = 0;
info.h = (b(:)' - a(:)')/(m - 1);
info.X1 = X1;  info.X2 = X2;  info.F = F;
info.history = [X1(j), X2(j), Fs];
info.cols = {'x1', 'x2', 'f'};
info.flag = 0;
info.message = sprintf('avaliou a grelha %d x %d (%d avaliações)', m, m, info.nfev);

if verbose
  fprintf('%4s %10s %10s %12s\n', 'i', 'x1', 'x2', 'f');
  for i = 1:min(5, numel(Fs))
    fprintf('%4d %10.4f %10.4f %12.4f\n', i, info.history(i, :));
  end
end
end
