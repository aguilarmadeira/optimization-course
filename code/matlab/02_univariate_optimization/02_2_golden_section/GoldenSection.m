function [x, fx, info] = GoldenSection(f, a, b, N, verbose)
%GOLDENSECTION  Método da secção áurea: min f(x) em [a,b], com N reduções.
%
%   [x, fx, info] = GoldenSection(f, a, b, N)
%   [x, fx, info] = GoldenSection(f, a, b, N, verbose)
%
%   Entradas
%     f        função (handle), unimodal em [a,b]
%     a, b     extremos do intervalo inicial, a < b
%     N        número de reduções do intervalo (L_N = tau^N * (b - a));
%              para garantir L_N <= tolx: N = ceil(log(tolx/(b-a))/log(tau))
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        ponto médio do intervalo final, |x - x*| <= L_N/2
%     fx       f(x)
%     info     estrutura com
%       .nfev      avaliações de f: N + 3 (N + 2 para localizar x*, mais f(x))
%       .nfev_loc  avaliações de f usadas para localizar x*: N + 2
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       número de reduções feitas (= N)
%       .history   tabela das iterações, uma linha por k = 0, 1, ..., N:
%                  [k  a  b  x1  x2  f(x1)  f(x2)  L_k]
%       .cols      nomes das colunas de .history
%       .flag      0 (fez as N reduções)
%       .message   mensagem de paragem
%
%   Otimização — deck 2.2 (Método da secção áurea).
%   Reproduz o exemplo dos slides: ver ex02_2_golden_section.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não verifica a unimodalidade de f).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 5, verbose = false; end
if ~(a < b), error('GoldenSection: é preciso a < b.'); end

nfev = 0;
hist = zeros(N + 1, 8);

% --- núcleo (igual ao dos slides, com as contagens) ----------------------
tau = (sqrt(5)-1)/2;  L = b - a;
x1 = b - tau*L;  f1 = f(x1);    % x1 < x2
x2 = a + tau*L;  f2 = f(x2);
nfev = nfev + 2;
hist(1, :) = [0, a, b, x1, x2, f1, f2, L];
for k = 1:N
  if f1 < f2                    % minimo em [a,x2]
    b = x2;  x2 = x1;  f2 = f1;
    L = b - a;  x1 = b - tau*L;  f1 = f(x1);
  else                          % minimo em [x1,b]
    a = x1;  x1 = x2;  f1 = f2;
    L = b - a;  x2 = a + tau*L;  f2 = f(x2);
  end
  nfev = nfev + 1;              % uma avaliação nova por redução
  hist(k + 1, :) = [k, a, b, x1, x2, f1, f2, L];
end
x = (a + b)/2;
% -------------------------------------------------------------------------

info.nfev_loc = nfev;           % N + 2: avaliações para localizar x*
fx = f(x);                      % avaliação no ponto médio final
nfev = nfev + 1;                % N + 3

info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = N;
info.history = hist;
info.cols = {'k', 'a', 'b', 'x1', 'x2', 'f(x1)', 'f(x2)', 'L_k'};
info.flag = 0;
info.message = sprintf('fez as %d reduções pedidas: L_N = %.4g', N, L);

if verbose
  fprintf('%3s %9s %9s %9s %9s %10s %10s %9s\n', 'k', 'a', 'b', 'x1', 'x2', ...
          'f(x1)', 'f(x2)', 'L_k');
  for i = 1:size(hist, 1)
    fprintf('%3d %9.4f %9.4f %9.4f %9.4f %10.4f %10.4f %9.4f\n', hist(i, :));
  end
end
end
