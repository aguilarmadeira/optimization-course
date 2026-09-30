function [x, fx, info] = HookeJeeves(f, x0, a, P0, T, kmax, verbose, menor)
%HOOKEJEEVES  Método de Hooke–Jeeves (pesquisa em padrão), sem derivadas.
%
%   [x, fx, info] = HookeJeeves(f, x0, a, P0, T)
%   [x, fx, info] = HookeJeeves(f, x0, a, P0, T, kmax, verbose)
%   [x, fx, info] = HookeJeeves(f, x0, a, P0, T, kmax, verbose, menor)
%
%   Alterna o movimento exploratório (Exploratory.m, eixo a eixo: avalia +-P_j
%   e fica com o melhor, como em Deb, 2012) com o
%   movimento de padrão xt = xb + a (xe - xb), seguido de nova exploração em
%   torno de xt. Se a exploração em torno de xb falha, P <- P/2 (P só diminui:
%   não é reposto em P0 a cada sucesso). Para quando todos os P_j < T_j.
%
%   Entradas
%     f        função (handle) de um vetor linha x (1 x n)
%     x0       ponto inicial (vetor; é tratado como linha)
%     a        fator do movimento de padrão (tipicamente 2); [] = 2
%     P0       perturbações iniciais, uma por coordenada (escalar = igual em todas)
%     T        tolerâncias, uma por coordenada (escalar = igual em todas)
%     kmax     (opcional, 10000) número máximo de movimentos (explorações
%              a partir da base + movimentos de padrão)
%     verbose  (opcional, false) se true, imprime uma linha por movimento
%     menor    (opcional, @(u,v) u < v) o valor u é melhor do que v? (ver
%              Exploratory.m; deck 4.3.3: regra de admissibilidade, com f a
%              devolver [v(x) f(x)])
%
%   Saídas
%     x        ponto base final (vetor linha)
%     fx       f(x) (já avaliado; não há avaliação extra no fim)
%     info     estrutura com
%       .nfev     avaliações de f: 1 (f(x0)) + 2n por exploração
%                 + 1 por ponto tentativo xt (cada padrão gasta 1 + 2n)
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit      número de movimentos (linhas de .history)
%       .history  uma linha por movimento:
%                 [m  tipo  xb(1:n)  xt(1:n)  f(xt)  xe(1:n)  f(xe)  f_ref  P(1:n)  sucesso  nfev]
%                 tipo = 0: exploração em torno de xb (xt e f(xt) = NaN; f_ref = f(xb));
%                 tipo = 1: padrão a partir da base xb, com ponto tentativo xt e
%                 exploração em torno de xt até xe' (nas colunas xe, f(xe)); f_ref = f(xe)
%                 anterior; sucesso = 1 se f(xe) < f_ref; nfev acumulado
%       .cols     nomes das colunas de .history
%       .ftrace   valores de f pela ordem em que foram avaliados (1 x nfev);
%                 p. ex. find(info.ftrace < 1e-4, 1) = avaliações até f < 1e-4
%       .P        perturbações no fim
%       .flag     0 se todos os P_j < T_j; 1 se atingiu kmax
%       .message  mensagem de paragem
%
%   Otimização — deck 3.2.5 (Método de Hooke–Jeeves).
%   Reproduz os exemplos dos slides: ver ex03_2_5_hooke_jeeves.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não escala as variáveis: P0 e T devem ser
%   escolhidos relativamente à escala de cada variável). Parar a passo
%   finito não certifica estacionariedade.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if isempty(a), a = 2; end
if nargin < 6 || isempty(kmax), kmax = 10000; end
if nargin < 7 || isempty(verbose), verbose = false; end
if nargin < 8 || isempty(menor), menor = @(u, v) u < v; end   % a ordem habitual

xb = x0(:).';
n = numel(xb);
P = P0(:).' .* ones(1, n);
T = T(:).' .* ones(1, n);
if any(P <= 0) || any(T <= 0), error('HookeJeeves: P0 e T têm de ser positivos.'); end

fb = f(xb);  nfev = 1;           % f(x0) conta
ftrace = fb(end);
hist = zeros(0, 4*n + 7);
m = 0;
flag = 1;
xe = xb;  fe = fb;

while m < kmax                   % arranque / reinício
  [xe, fe, ie] = Exploratory(f, xb, fb, P, menor);   % exploração em torno de xb
  nfev = nfev + ie.nfev;  ftrace = [ftrace, ie.ftrace];
  m = m + 1;
  hist(m, :) = [m, 0, xb, NaN(1, n), NaN, xe, fe(end), fb(end), P, menor(fe, fb), nfev];
  if ~menor(fe, fb)              % falhou (fe >= fb)
    P = P/2;
    if all(P < T), flag = 0; break; end
  else
    while m < kmax               % movimentos de padrão
      xt = xb + a*(xe - xb);  ft = f(xt);            % padrão: ponto tentativo
      nfev = nfev + 1;  ftrace(end + 1) = ft(end);
      [xn, fn, ie] = Exploratory(f, xt, ft, P, menor);   % exploração em torno de xt
      nfev = nfev + ie.nfev;  ftrace = [ftrace, ie.ftrace];
      m = m + 1;
      hist(m, :) = [m, 1, xb, xt, ft(end), xn, fn(end), fe(end), P, menor(fn, fe), nfev];
      if ~menor(fn, fe)          % rejeitado (fn >= fe): recuar para xe e reiniciar
        xb = xe;  fb = fe;
        break
      else                       % aceite: continuar o padrão
        xb = xe;  fb = fe;
        xe = xn;  fe = fn;
      end
    end
  end
end
if menor(fe, fb)                 % saiu por kmax a meio de um padrão aceite
  xb = xe;  fb = fe;
end
x = xb;  fx = fb;

cols = cell(1, 4*n + 7);
cols(1:2) = {'m', 'tipo'};
for j = 1:n
  cols{2 + j} = sprintf('xb%d', j);
  cols{2 + n + j} = sprintf('xt%d', j);
  cols{3 + 2*n + j} = sprintf('xe%d', j);
  cols{5 + 3*n + j} = sprintf('P%d', j);
end
cols{3 + 2*n} = 'f(xt)';
cols(4 + 3*n:5 + 3*n) = {'f(xe)', 'f_ref'};
cols(6 + 4*n:7 + 4*n) = {'sucesso', 'nfev'};

info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = m;
info.history = hist;
info.cols = cols;
info.ftrace = ftrace;
info.P = P;
info.flag = flag;
if flag == 0
  info.message = sprintf('todos os P_j < T_j (max P_j = %.3g) ao fim de %d movimentos', max(P), m);
else
  info.message = sprintf('atingiu kmax = %d movimentos (max P_j = %.3g)', kmax, max(P));
end

if verbose
  fprintf('%4s %-9s', 'm', 'movimento');
  fprintf('%9s', cols{3:2+n});
  fprintf('%9s', cols{3+n:2+2*n});
  fprintf('%10s', 'f(xt)');
  fprintf('%9s', cols{4+2*n:3+3*n});
  fprintf('%10s', 'f(xe)');
  fprintf('%9s', cols{6+3*n:5+4*n});
  fprintf('  %-10s %5s\n', 'resultado', 'n_f');
  for r = 1:m
    h = hist(r, :);
    tipo = h(2);
    if tipo == 0, fprintf('%4d explor.  ', h(1));
    else,         fprintf('%4d padrão   ', h(1));
    end
    fprintf('%9.4f', h(3:2+n));
    if tipo == 0
      fprintf('%s', repmat(sprintf('%9s', '-'), 1, n));
      fprintf('%10s', '-');
    else
      fprintf('%9.4f', h(3+n:2+2*n));
      fprintf('%10.4f', h(3+2*n));
    end
    fprintf('%9.4f', h(4+2*n:3+3*n));
    fprintf('%10.4f', h(4+3*n));
    fprintf('%9.4f', h(6+3*n:5+4*n));
    if tipo == 0
      if h(6+4*n), res = 'melhora'; else, res = 'falha:P/2'; end
    else
      if h(6+4*n), res = 'aceite'; else, res = 'rejeitado'; end
    end
    fprintf('  %-10s %5d\n', res, h(7+4*n));
  end
end
end
