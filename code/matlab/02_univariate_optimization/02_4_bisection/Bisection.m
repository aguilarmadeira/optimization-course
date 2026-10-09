function [x, fx, info] = Bisection(f, df, a, b, tolg, kmax, tolx, verbose)
%BISECTION  Método da bisseção sobre f': min f(x) em [a,b], como nas aulas (Deb).
%
%   [x, fx, info] = Bisection(f, df, a, b, tolg)
%   [x, fx, info] = Bisection(f, df, a, b, tolg, kmax, tolx, verbose)
%
%   Passo 1: escolher a < b com f'(a) < 0 < f'(b) e eps > 0; fazer x1 = a, x2 = b.
%   Passo 2: fazer z = (x1 + x2)/2 e calcular f'(z).
%   Passo 3: se |f'(z)| <= eps, terminar (x* ~ z); senão, se f'(z) < 0,
%            fazer x1 = z; se f'(z) > 0, fazer x2 = z; voltar ao passo 2.
%
%   Entradas
%     f        função (handle) ou [] -- só para leitura (ver abaixo)
%     df       derivada f' (handle); exige df(a) < 0 < df(b)
%     a, b     extremos do intervalo inicial, a < b
%     tolg     tolerância eps da paragem |f'(z)| <= eps
%     kmax     (opcional, 100) número máximo de iterações
%     tolx     (opcional, 0) variante: parar também quando x2 - x1 <= tolx e
%              devolver o ponto médio (|x - x*| <= tolx/2). Com tolg = 0 é a
%              bisseção com N = ceil(log2((b - a)/tolx)) reduções.
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        z da última iteração (na variante, o ponto médio final)
%     fx       f(x) (NaN se f = [])
%     info     estrutura com
%       .ngev      avaliações de f': nit + 2 (df(a), df(b) e uma por iteração)
%       .nfev      0: o método só usa f'
%       .nhev      0
%       .nit       número de iterações (avaliações de f'(z))
%       .L         x2 - x1 no fim
%       .history   tabela das iterações, uma linha por k = 1, ..., nit:
%                  [k  x1  x2  z  f(z)  f'(z)]  (x1, x2 no início da iteração)
%       .cols      nomes das colunas de .history
%       .flag      0 (|f'(z)| <= tolg), 1 (atingiu kmax) ou 2 (x2 - x1 <= tolx)
%       .message   mensagem de paragem
%
%   f só é usada para leitura: na coluna f(z) da tabela e em fx. O método
%   não usa f, e por isso essas avaliações não entram em info.nfev (n_f = 0,
%   como nos slides). Com f = [] nada disso é calculado.
%
%   Otimização — deck 2.4 (Bisseção).
%   Reproduz o exemplo dos slides: ver ex02_4_bisection.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não verifica a unicidade do zero de f').
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 6 || isempty(kmax), kmax = 100; end
if nargin < 7 || isempty(tolx), tolx = 0; end
if nargin < 8, verbose = false; end
if ~(a < b), error('Bisection: é preciso a < b.'); end
temf = ~isempty(f);

% --- verificar o enquadramento: 2 avaliações de f' -----------------------
da = df(a);  db = df(b);
ngev = 2;
if ~(da<0 && db>0), error('enquadramento'), end

hist = zeros(0, 6);
flag = 1;
x1 = a;  x2 = b;
z = (x1 + x2)/2;
% --- núcleo (igual ao dos slides, com as contagens) ----------------------
for k = 1:kmax
    z = (x1 + x2)/2;  d = df(z);                % passo 2
    ngev = ngev + 1;
    if temf, fk = f(z); else, fk = NaN; end      % só para leitura
    hist(end + 1, :) = [k, x1, x2, z, fk, d]; %#ok<AGROW>
    if abs(d) <= tolg                            % passo 3: terminou
        flag = 0;
        break
    end
    if d < 0
        x1 = z;                                  % minimo a direita
    else
        x2 = z;                                  % minimo a esquerda
    end
    if tolx > 0 && x2 - x1 <= tolx               % variante
        flag = 2;
        break
    end
end
% -------------------------------------------------------------------------

if flag == 2, x = (x1 + x2)/2; else, x = z; end
if temf, fx = f(x); else, fx = NaN; end          % só para leitura

info.ngev = ngev;
info.nfev = 0;
info.nhev = 0;
info.nit = size(hist, 1);
info.L = x2 - x1;
info.history = hist;
info.cols = {'k', 'x1', 'x2', 'z', 'f(z)', 'df(z)'};
info.flag = flag;
switch flag
  case 0
    info.message = sprintf('|f''(z)| <= tolg na iteração %d', info.nit);
  case 2
    info.message = sprintf('x2 - x1 = %.2g <= tolx ao fim de %d iterações (ponto médio)', info.L, info.nit);
  otherwise
    info.message = sprintf('atingiu o número máximo de iterações (kmax = %d)', kmax);
    warning('Bisection: %s', info.message);
end

if verbose
  fprintf('%3s %9s %9s %9s %10s %10s\n', 'k', 'x1', 'x2', 'z', 'f(z)', 'df(z)');
  for r = 1:size(hist, 1)
    fprintf('%3d %9.4f %9.4f %9.4f %10.4f %+10.4f\n', hist(r, :));
  end
end
end
