function [x, fx, info] = NewtonNumeric(f, x0, tolg, kmax, hrel, verbose)
%NEWTONNUMERIC  Newton 1D com f' e f'' por diferenças centradas (só usa f).
%
%   [x, fx, info] = NewtonNumeric(f, x0, tolg, kmax)
%   [x, fx, info] = NewtonNumeric(f, x0, tolg, kmax, hrel, verbose)
%
%   Em cada iteração: h = max(hrel*|x_k|, 1e-4); com f(x_k - h), f(x_k),
%   f(x_k + h) calcula df ~ f'(x_k) e ddf ~ f''(x_k) e dá o passo de Newton
%   x_{k+1} = x_k - df/ddf. Para quando |df| < tolg (df calculado em x_k) e
%   devolve x_{k+1}.
%
%   Entradas
%     f        função (handle)
%     x0       ponto inicial
%     tolg     tolerância em |f'_num(x_k)|
%     kmax     número máximo de iterações (com aviso se for atingido)
%     hrel     (opcional, 0.01) passo relativo: com 0.01 é a regra da UC,
%              h = max(0.01|x|, 1e-4)
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        último iterando x_{k+1}
%     fx       f(x)
%     info     estrutura com
%       .nfev      avaliações de f: 3 por iteração + 1 (f(x) no fim)
%       .ngev, .nhev   0 (as derivadas são numéricas: contam em nfev)
%       .nit       número de iterações (passos de Newton)
%       .history   tabela das iterações, uma linha por k = 0, 1, ..., nit-1:
%                  [k  x_k  h  df  ddf  x_{k+1}]
%       .cols      nomes das colunas de .history
%       .flag      0 (|df| < tolg) ou 1 (atingiu kmax)
%       .message   mensagem de paragem
%
%   Com h fixo o método estagna perto de x*: converge para o zero de f'_num,
%   que não é o zero de f' (erro de truncatura ~ h^2 |f'''|/6).
%
%   Otimização — deck 2.6 (Derivadas numéricas; Newton numérico).
%   Reproduz o exemplo dos slides: ver ex02_6_finite_differences.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. ddf <= 0, passo amortecido).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 5 || isempty(hrel), hrel = 0.01; end
if nargin < 6, verbose = false; end

x = x0;
nfev = 0;
hist = zeros(0, 6);
flag = 1;
for k = 1:kmax
  % --- núcleo (igual ao dos slides) ---------------------------------------
  h  = max(hrel*abs(x), 1e-4);
  fm = f(x - h);  f0 = f(x);  fp = f(x + h);
  df  = (fp - fm) / (2*h);
  ddf = (fp - 2*f0 + fm) / h^2;
  if ddf == 0, error('f''''_num = 0: aumentar h'), end
  xk = x;
  x   = x - df/ddf;              % passo de Newton
  % -----------------------------------------------------------------------
  nfev = nfev + 3;
  hist(end + 1, :) = [k - 1, xk, h, df, ddf, x]; %#ok<AGROW>
  if abs(df) < tolg, flag = 0; break, end
end
fx = f(x);
nfev = nfev + 1;

info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = size(hist, 1);
info.history = hist;
info.cols = {'k', 'x_k', 'h', 'df_num', 'ddf_num', 'x_k+1'};
info.flag = flag;
if flag == 0
  info.message = sprintf('|f''_num(x_k)| < tolg em k = %d; devolve x_%d', info.nit - 1, info.nit);
else
  info.message = sprintf('atingiu o número máximo de iterações (kmax = %d) sem |f''_num| < tolg', kmax);
  warning('NewtonNumeric: %s', info.message);
end

if verbose
  fprintf('%3s %11s %8s %14s %12s %11s\n', info.cols{:});
  for r = 1:size(hist, 1)
    fprintf('%3d %11.7f %8.4f %+14.6e %12.6f %11.7f\n', hist(r, :));
  end
end
end
