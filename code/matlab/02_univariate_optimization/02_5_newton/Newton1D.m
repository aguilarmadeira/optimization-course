function [x, fx, info] = Newton1D(f, df, ddf, x0, tolg, kmax, tolH, verbose)
%NEWTON1D  Método de Newton-Raphson em 1D: procura um zero de f'.
%
%   [x, fx, info] = Newton1D(f, df, ddf, x0, tolg, kmax)
%   [x, fx, info] = Newton1D(f, df, ddf, x0, tolg, kmax, tolH, verbose)
%
%   Entradas
%     f        função (handle) ou [] -- só para leitura (ver abaixo)
%     df, ddf  primeira e segunda derivadas f', f'' (handles)
%     x0       ponto inicial
%     tolg     tolerância no gradiente: para quando |f'(x_{k+1})| < tolg
%     kmax     número máximo de iterações (com aviso se for atingido)
%     tolH     (opcional, 1e-12) erro se |f''(x_k)| < tolH (curvatura ~ 0)
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        última iterada x_k
%     fx       f(x) (NaN se f = [])
%     info     estrutura com
%       .ngev      avaliações de f': 1 + nit (f'(x0) e uma por iteração)
%       .nhev      avaliações de f'': nit (uma por iteração)
%       .nfev      0: o método não usa f
%       .nit       número de iterações (passos de Newton)
%       .history   tabela das iterações, uma linha por k = 0, 1, ..., nit-1:
%                  [k  x_k  f(x_k)  f'(x_k)  f''(x_k)  x_{k+1}  f'(x_{k+1})]
%       .cols      nomes das colunas de .history
%       .flag      0 (|f'| < tolg) ou 1 (atingiu kmax)
%       .message   mensagem de paragem
%
%   Newton procura zeros de f', não mínimos: no fim, classificar o ponto com
%   f''(x) (> 0 mínimo, < 0 máximo). Essa avaliação extra não é feita aqui.
%   f só é usada para leitura (coluna f(x_k) e fx); essas avaliações não
%   entram em info.nfev (n_f = 0, como nos slides).
%
%   Otimização — deck 2.5 (Newton-Raphson).
%   Reproduz o exemplo dos slides: ver ex02_5_newton.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. passo amortecido ou salvaguarda por bisseção).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 7 || isempty(tolH), tolH = 1e-12; end
if nargin < 8, verbose = false; end
temf = ~isempty(f);

x = x0;
hist = zeros(0, 7);
flag = 1;
ngev = 0;  nhev = 0;
% --- núcleo (igual ao dos slides, com as contagens) ----------------------
g = df(x);
ngev = ngev + 1;
for k = 1:kmax
    H = ddf(x);
    nhev = nhev + 1;
    if abs(H) < tolH
        error('curvatura demasiado pequena')
    end
    xk = x;  gk = g;
    x = x - g/H;                 % passo de Newton
    g = df(x);                   % 1 avaliacao de f'
    ngev = ngev + 1;
    if temf, fk = f(xk); else, fk = NaN; end      % só para leitura
    hist(end + 1, :) = [k - 1, xk, fk, gk, H, x, g]; %#ok<AGROW>
    if abs(g) < tolg, flag = 0; break, end
end
% -------------------------------------------------------------------------

if temf, fx = f(x); else, fx = NaN; end          % só para leitura

info.ngev = ngev;
info.nhev = nhev;
info.nfev = 0;
info.nit = size(hist, 1);
info.history = hist;
info.cols = {'k', 'x_k', 'f(x_k)', 'df(x_k)', 'ddf(x_k)', 'x_k+1', 'df(x_k+1)'};
info.flag = flag;
if flag == 0
  info.message = sprintf('|f''(x)| < tolg em %d iterações', info.nit);
else
  info.message = sprintf('atingiu o número máximo de iterações (kmax = %d) sem |f''| < tolg', kmax);
  warning('Newton1D: %s', info.message);
end

if verbose
  fprintf('%3s %13s %13s %14s %12s %13s %14s\n', info.cols{:});
  for r = 1:size(hist, 1)
    fprintf('%3d %13.7f %13.6f %+14.6e %12.5f %13.7f %+14.6e\n', hist(r, :));
  end
end
end
