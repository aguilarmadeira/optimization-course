function [x, fx, info] = NewtonSafeguarded(f, df, ddf, a, b, x0, tolg, kmax, verbose)
%NEWTONSAFEGUARDED  Newton salvaguardado (bisseção + Newton): zero de f' em [a,b].
%
%   [x, fx, info] = NewtonSafeguarded(f, df, ddf, a, b, x0, tolg, kmax)
%   [x, fx, info] = NewtonSafeguarded(f, df, ddf, a, b, x0, tolg, kmax, verbose)
%
%   Mantém um intervalo [a,b] com f'(a) < 0 < f'(b). Em cada iteração:
%   (1) tenta o passo de Newton x_N = x_k - f'(x_k)/f''(x_k);
%   (2) se x_N sair de (a,b), ou se não der progresso suficiente
%       (|x_N - x_k| maior do que metade do passo de duas iterações antes,
%       o teste de rtsafe), faz um passo de bisseção, x_{k+1} = (a+b)/2;
%   (3) atualiza o intervalo com o sinal de f'(x_{k+1}): negativo -> a,
%       senão -> b.
%
%   Entradas
%     f        função (handle) ou [] -- só para leitura (fx)
%     df, ddf  primeira e segunda derivadas f', f'' (handles)
%     a, b     intervalo inicial, com f'(a) < 0 < f'(b)
%     x0       ponto inicial em [a,b]
%     tolg     tolerância no gradiente: para quando |f'(x_{k+1})| < tolg
%     kmax     número máximo de iterações (com aviso se for atingido)
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        última iterada x_k
%     fx       f(x) (NaN se f = [])
%     info     estrutura com
%       .ngev      avaliações de f': 3 + nit (f'(a), f'(b) para verificar o
%                  enquadramento, f'(x0) e uma por iteração)
%       .nhev      avaliações de f'': nit (uma por iteração)
%       .nfev      0: o método não usa f
%       .nit       número de iterações
%       .nbis      passos de bisseção
%       .a, .b     intervalo final
%       .history   uma linha por k = 0, 1, ..., nit-1:
%                  [k  x_k  x_N  tipo  x_{k+1}  a  b  f'(x_{k+1})]
%                  tipo = 0 (Newton) ou 1 (bisseção); [a, b] depois do passo
%       .cols      nomes das colunas de .history
%       .flag      0 (|f'| < tolg) ou 1 (atingiu kmax)
%       .message   mensagem de paragem
%
%   Otimização — deck 2.5 (Newton-Raphson), página «Newton salvaguardado».
%   Reproduz o exemplo dos slides: ver ex02_5_newton_safeguarded.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. tolerância no comprimento do intervalo).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 9, verbose = false; end
ga = df(a);  gb = df(b);
ngev = 2;
if ~(ga < 0 && gb > 0)
  error('NewtonSafeguarded: é preciso f''(a) < 0 < f''(b)');
end
if ~(a <= x0 && x0 <= b)
  error('NewtonSafeguarded: x0 tem de estar em [a,b]');
end

x = x0;
hist = zeros(0, 8);
flag = 1;
nhev = 0;  nbis = 0;
dxold = b - a;  dx = dxold;
% --- núcleo (igual ao dos slides, com as contagens) ----------------------
g = df(x);
ngev = ngev + 1;
for k = 1:kmax
  H = ddf(x);
  nhev = nhev + 1;
  if H ~= 0, xN = x - g/H; else, xN = Inf; end    % (1) passo de Newton
  if ~(a < xN && xN < b) || abs(xN - x) > abs(dxold)/2
    xn = (a + b)/2;                               % (2) bisseção
    tipo = 1;  nbis = nbis + 1;
  else
    xn = xN;
    tipo = 0;
  end
  dxold = dx;  dx = xn - x;
  gn = df(xn);                                    % 1 avaliação de f'
  ngev = ngev + 1;
  if gn < 0                                       % (3) novo intervalo
    a = xn;
  else
    b = xn;
  end
  hist(end + 1, :) = [k - 1, x, xN, tipo, xn, a, b, gn]; %#ok<AGROW>
  x = xn;  g = gn;
  if abs(g) < tolg, flag = 0; break, end
end
% -------------------------------------------------------------------------

if isempty(f), fx = NaN; else, fx = f(x); end      % só para leitura

info.ngev = ngev;
info.nhev = nhev;
info.nfev = 0;
info.nit = size(hist, 1);
info.nbis = nbis;
info.a = a;  info.b = b;
info.history = hist;
info.cols = {'k', 'x_k', 'x_Newton', 'tipo', 'x_k+1', 'a', 'b', 'df(x_k+1)'};
info.flag = flag;
if flag == 0
  info.message = sprintf('|f''(x)| < tolg em %d iterações (%d de bisseção)', info.nit, nbis);
else
  info.message = sprintf('atingiu o número máximo de iterações (kmax = %d) sem |f''| < tolg', kmax);
  warning('NewtonSafeguarded: %s', info.message);
end

if verbose
  tipos = {'   Newton', ' bisseção'};         % largura visível 9 (UTF-8)
  fprintf('%3s %12s %12s %9s %12s %24s %11s\n', 'k', 'x_k', 'x_Newton', 'tipo', 'x_k+1', '[a, b]', 'df(x_k+1)');
  for r = 1:size(hist, 1)
    fprintf('%3d %12.4e %12.4e %s %12.4e  [%10.3e; %10.3e] %+11.2e\n', hist(r, 1:3), ...
            tipos{hist(r, 4) + 1}, hist(r, 5:8));
  end
end
end
