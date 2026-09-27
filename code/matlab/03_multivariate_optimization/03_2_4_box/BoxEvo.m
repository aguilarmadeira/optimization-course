function [x, fx, info] = BoxEvo(f, x0, Delta, tolx, kmax, verbose)
%BOXEVO  Método de Box (versão determinística simplificada de EVOP).
%
%   [x, fx, info] = BoxEvo(f, x0, Delta, tolx)
%   [x, fx, info] = BoxEvo(f, x0, Delta, tolx, kmax, verbose)
%
%   Em cada iteração avalia f nos 2^n vértices x + (1/2)(±Delta_1,...,±Delta_n)
%   da caixa centrada no ponto atual; se o melhor vértice melhora f(x), move-se
%   para ele (Delta fica igual); senão encolhe a caixa, Delta <- Delta/2.
%   Para quando norm(Delta) < tolx.
%
%   Entradas
%     f        função (handle) de um vetor linha x (1 x n)
%     x0       ponto inicial (vetor; é tratado como linha)
%     Delta    lados da caixa, um por coordenada (escalar = igual em todas)
%     tolx     tolerância: para quando norm(Delta) < tolx
%     kmax     (opcional, 10000) número máximo de iterações
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x        último centro da caixa (vetor linha)
%     fx       f(x) (já avaliado; não há avaliação extra no fim)
%     info     estrutura com
%       .nfev     avaliações de f: 1 + 2^n * nit (f(x0) conta; avaliam-se
%                 SEMPRE os 2^n vértices, sem reaproveitar pontos — é a
%                 convenção das contagens do deck 3.2.4)
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit      número de iterações (cada uma: 2^n avaliações e mover ou encolher)
%       .history  uma linha por iteração k = 0, ..., nit-1 (o centro x_k e a
%                 decisão tomada a partir dele, como na tabela do deck):
%                 [k  x_k(1:n)  Delta(1:n)  f(x_k)  decisão  f_novo  nfev]
%                 decisão = 1 (mover) ou 0 (encolher); nfev acumulado no fim da iteração
%       .cols     nomes das colunas de .history
%       .ftrace   valores de f pela ordem em que foram avaliados (1 x nfev);
%                 p. ex. find(info.ftrace < 1e-4, 1) = avaliações até f < 1e-4
%       .Delta    lados da caixa no fim
%       .flag     0 se norm(Delta) < tolx; 1 se atingiu kmax
%       .message  mensagem de paragem
%
%   Otimização — deck 3.2.4 (Método de Box — versão simplificada de EVOP).
%   Reproduz os exemplos dos slides: ver ex03_2_4_box.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não reaproveita vértices já avaliados e o
%   custo 2^n por iteração torna-o impraticável para n grande).
%   Uma caixa pequena não certifica estacionariedade a Delta finito.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 5 || isempty(kmax), kmax = 10000; end
if nargin < 6, verbose = false; end

x = x0(:).';                     % vetor linha, como no núcleo do slide
n = numel(x);
Delta = Delta(:).' .* ones(1, n);
if any(Delta <= 0), error('BoxEvo: os lados Delta têm de ser positivos.'); end

fx = f(x);  nfev = 1;            % f(x0) conta
ftrace = zeros(1, 1 + 2^n * min(kmax, 1000));
ftrace(1) = fx;
hist = zeros(0, 2*n + 5);

S = dec2bin(0:2^n-1) - '0';
S = 2*S - 1;                     % linhas de +-1 (em 2D: (-1,-1), (-1,1), (1,-1), (1,1))
k = 0;
while norm(Delta) >= tolx && k < kmax
  xk = x;  Dk = Delta;  fk = fx;
  % --- núcleo (igual ao dos slides, com as contagens) --------------------
  V  = x + 0.5*S.*Delta;         % vertices
  Fv = arrayfun(@(i) f(V(i,:)), 1:2^n);
  [fmin, i] = min(Fv);
  if fmin < fx
    x = V(i,:);  fx = fmin;      % mover
    dec = 1;
  else
    Delta = Delta/2;             % encolher
    dec = 0;
  end
  % -----------------------------------------------------------------------
  ftrace(nfev + (1:2^n)) = Fv;
  nfev = nfev + 2^n;             % avaliam-se sempre os 2^n vértices
  hist(k + 1, :) = [k, xk, Dk, fk, dec, fx, nfev];
  k = k + 1;
end
ftrace = ftrace(1:nfev);

cols = cell(1, 2*n + 5);
cols{1} = 'k';
for j = 1:n
  cols{1 + j} = sprintf('x%d', j);
  cols{1 + n + j} = sprintf('Delta%d', j);
end
cols(2*n + 2:end) = {'f(x_k)', 'decisão', 'f_novo', 'nfev'};

info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = k;
info.history = hist;
info.cols = cols;
info.ftrace = ftrace;
info.Delta = Delta;
if norm(Delta) < tolx
  info.flag = 0;
  info.message = sprintf('norm(Delta) = %.3g < tolx = %.3g ao fim de %d iterações', ...
                         norm(Delta), tolx, k);
else
  info.flag = 1;
  info.message = sprintf('atingiu kmax = %d iterações (norm(Delta) = %.3g)', kmax, norm(Delta));
end

if verbose
  fprintf('%5s', 'k');
  fprintf('%10s', cols{2:2*n+1});
  fprintf('%12s   decisão %12s %7s\n', 'f(x_k)', 'f_novo', 'n_f');
  nomes = {'encolhe', 'move'};
  for r = 1:size(hist, 1)
    fprintf('%5d', hist(r, 1));
    fprintf('%10.4f', hist(r, 2:2*n+1));
    fprintf('%12.4f %9s %12.4f %7d\n', hist(r, 2*n+2), nomes{hist(r, 2*n+3) + 1}, ...
            hist(r, 2*n+4), hist(r, 2*n+5));
  end
end
end
