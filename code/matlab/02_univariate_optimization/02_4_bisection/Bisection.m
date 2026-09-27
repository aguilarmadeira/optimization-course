function [x, fx, info] = Bisection(f, df, a, b, N, verbose)
%BISECTION  Método da bisseção sobre f': min f(x) em [a,b], com N reduções.
%
%   [x, fx, info] = Bisection(f, df, a, b, N)
%   [x, fx, info] = Bisection(f, df, a, b, N, verbose)
%
%   Entradas
%     f        função (handle) ou [] -- só para leitura (ver abaixo)
%     df       derivada f' (handle); exige df(a) < 0 < df(b)
%     a, b     extremos do intervalo inicial, a < b
%     N        número de reduções (b_N - a_N = (b - a)/2^N);
%              para garantir b_N - a_N <= tolx: N = ceil(log2((b - a)/tolx))
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        ponto médio do intervalo final, |x - x*| <= (b_N - a_N)/2
%     fx       f(x) (NaN se f = [])
%     info     estrutura com
%       .ngev      avaliações de f': N + 2 (df(a), df(b) e uma por redução)
%       .nfev      0: o método só usa o sinal de f'
%       .nhev      0
%       .nit       número de reduções feitas (N, ou menos se f'(x_k) = 0)
%       .L         comprimento do intervalo final, b_N - a_N
%       .history   tabela das iterações, uma linha por k = 0, 1, ..., N-1:
%                  [k  a_k  b_k  x_k  f(x_k)  f'(x_k)]
%       .cols      nomes das colunas de .history
%       .flag      0 (fez as N reduções) ou 2 (f'(x_k) = 0: terminou)
%       .message   mensagem de paragem
%
%   f só é usada para leitura: na coluna f(x_k) da tabela e em fx. O método
%   não usa f, e por isso essas avaliações não entram em info.nfev (n_f = 0,
%   como nos slides). Com f = [] nada disso é calculado.
%
%   Otimização — deck 2.4 (Bisseção).
%   Reproduz o exemplo dos slides: ver ex02_4_bisection.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não verifica a unicidade do zero de f').
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 6, verbose = false; end
if ~(a < b), error('Bisection: é preciso a < b.'); end
temf = ~isempty(f);

% --- verificar o enquadramento: 2 avaliações de f' -----------------------
da = df(a);  db = df(b);
ngev = 2;
if ~(da<0 && db>0), error('enquadramento'), end

hist = zeros(0, 6);
flag = 0;
% --- núcleo (igual ao dos slides, com as contagens) ----------------------
for k = 1:N
    x = (a + b)/2;  d = df(x);
    ngev = ngev + 1;
    if temf, fk = f(x); else, fk = NaN; end      % só para leitura
    hist(end + 1, :) = [k - 1, a, b, x, fk, d]; %#ok<AGROW>
    if d == 0
        a = x;  b = x;  flag = 2;  break   % terminou
    elseif d < 0
        a = x;                  % minimo a direita
    else
        b = x;                  % minimo a esquerda
    end
end
x = (a + b)/2;
% -------------------------------------------------------------------------

if temf, fx = f(x); else, fx = NaN; end          % só para leitura

info.ngev = ngev;
info.nfev = 0;
info.nhev = 0;
info.nit = size(hist, 1);
info.L = b - a;
info.history = hist;
info.cols = {'k', 'a_k', 'b_k', 'x_k', 'f(x_k)', 'df(x_k)'};
info.flag = flag;
if flag == 2
  info.message = sprintf('f''(x_k) = 0 na iteração k = %d: terminou', info.nit - 1);
else
  info.message = sprintf('fez as %d reduções pedidas: b - a = %.2g', N, info.L);
end

if verbose
  fprintf('%3s %9s %9s %9s %10s %10s\n', 'k', 'a_k', 'b_k', 'x_k', 'f(x_k)', 'df(x_k)');
  for r = 1:size(hist, 1)
    fprintf('%3d %9.4f %9.4f %9.4f %10.4f %+10.4f\n', hist(r, :));
  end
end
end
