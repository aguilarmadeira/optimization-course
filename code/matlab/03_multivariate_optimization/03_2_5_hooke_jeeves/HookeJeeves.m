function [x, fx, info] = HookeJeeves(f, x0, a, P0, T, kmax, verbose, menor, variante, R)
%HOOKEJEEVES  Método de Hooke–Jeeves (pesquisa em padrão), sem derivadas.
%
%   [x, fx, info] = HookeJeeves(f, x0, a, P0, T)
%   [x, fx, info] = HookeJeeves(f, x0, a, P0, T, kmax, verbose)
%   [x, fx, info] = HookeJeeves(f, x0, a, P0, T, kmax, verbose, menor)
%   [x, fx, info] = HookeJeeves(f, x0, a, P0, T, kmax, verbose, menor, variante, R)
%
%   Alterna o movimento exploratório (Exploratory.m, eixo a eixo: avalia
%   +-Delta_j e fica com o melhor dos três pontos, como em Deb, 2012) com o
%   movimento de padrão X^P_{k+1} = 2 X_k - X_{k-1} (a = 2), seguido de nova
%   exploração em torno de X^P_{k+1}.
%   Variante 'deb' (por omissão; a das aulas, Deb, 2012, sec. 3.3.3):
%     passo 2  explorar em torno de X_k; se melhorar, X_{k+1} = resultado, passo 4;
%              senão, passo 3
%     passo 3  se ||Delta|| < eps, parar; senão Delta = Delta/R e voltar ao passo 2
%     passo 4  padrão X^P_{k+1} = 2 X_k - X_{k-1}
%     passo 5  explorar em torno de X^P_{k+1}; o resultado é X_{k+1}
%     passo 6  se f(X_{k+1}) < f(X_k), voltar ao passo 4; senão, passo 3
%   Aqui T é eps (um escalar) e R o fator de redução (2 por omissão).
%   Variante '1961' (Hooke e Jeeves, 1961): quando o padrão falha, recua-se
%   para X_k e explora-se de novo com o mesmo Delta; Delta só se divide por 2
%   quando a exploração em torno da base falha; pára quando todos os
%   Delta_j < T_j. Em ambas, Delta só diminui.
%
%   Entradas
%     f        função (handle) de um vetor linha x (1 x n)
%     x0       ponto inicial (vetor; é tratado como linha)
%     a        fator do movimento de padrão (tipicamente 2); [] = 2
%     P0       passos iniciais Delta_0, um por coordenada (escalar = igual em todas)
%     T        variante 'deb': eps (escalar), pára quando ||Delta|| < eps;
%              variante '1961': tolerâncias, uma por coordenada
%     kmax     (opcional, 10000) número máximo de movimentos (explorações
%              a partir da base + movimentos de padrão)
%     verbose  (opcional, false) se true, imprime uma linha por movimento
%     menor    (opcional, @(u,v) u < v) o valor u é melhor do que v? (ver
%              Exploratory.m; deck 4.3.3: regra de admissibilidade, com f a
%              devolver [v(x) f(x)])
%     variante (opcional, 'deb') 'deb' ou '1961'
%     R        (opcional, 2) fator de redução de Delta na variante 'deb'
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
%       .flag     0 se parou pelo critério de paragem; 1 se atingiu kmax
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
if nargin < 9 || isempty(variante), variante = 'deb'; end
if nargin < 10 || isempty(R), R = 2; end
if ~any(strcmp(variante, {'deb', '1961'})), error('HookeJeeves: variante tem de ser ''deb'' ou ''1961''.'); end
deb = strcmp(variante, 'deb');

xb = x0(:).';
n = numel(xb);
P = P0(:).' .* ones(1, n);
if deb
  if numel(T) ~= 1, error('HookeJeeves: na variante ''deb'', T é eps (um escalar).'); end
  epsD = T;  T = epsD*ones(1, n);
else
  T = T(:).' .* ones(1, n);
end
if any(P <= 0) || any(T <= 0), error('HookeJeeves: P0 e T têm de ser positivos.'); end

fb = f(xb);  nfev = 1;           % f(x0) conta
ftrace = fb(end);
hist = zeros(0, 4*n + 7);
m = 0;
flag = 1;
xe = xb;  fe = fb;

while m < kmax                   % passo 2: exploração em torno da base xb (= X_k)
  [xe, fe, ie] = Exploratory(f, xb, fb, P, menor);
  nfev = nfev + ie.nfev;  ftrace = [ftrace, ie.ftrace];
  m = m + 1;
  hist(m, :) = [m, 0, xb, NaN(1, n), NaN, xe, fe(end), fb(end), P, menor(fe, fb), nfev];
  if ~menor(fe, fb)              % falhou (fe >= fb): passo 3
    if deb
      if norm(P) < epsD, flag = 0; break; end
      P = P/R;
    else
      P = P/2;
      if all(P < T), flag = 0; break; end
    end
  else
    while m < kmax               % movimentos de padrão
      xt = xb + a*(xe - xb);  ft = f(xt);            % padrão: ponto tentativo
      nfev = nfev + 1;  ftrace(end + 1) = ft(end);
      [xn, fn, ie] = Exploratory(f, xt, ft, P, menor);   % exploração em torno de xt
      nfev = nfev + ie.nfev;  ftrace = [ftrace, ie.ftrace];
      m = m + 1;
      hist(m, :) = [m, 1, xb, xt, ft(end), xn, fn(end), fe(end), P, menor(fn, fe), nfev];
      if ~menor(fn, fe)          % rejeitado (fn >= fe): recuar para xe (= X_k)
        xb = xe;  fb = fe;
        if deb                   % passo 3 logo a seguir
          if norm(P) < epsD, flag = 0; else, P = P/R; end
        end
        break
      else                       % aceite: continuar o padrão
        xb = xe;  fb = fe;
        xe = xn;  fe = fn;
      end
    end
    if flag == 0, break; end
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
if flag == 0 && deb
  info.message = sprintf('||Delta|| = %.3g < eps = %.3g ao fim de %d movimentos', norm(P), epsD, m);
elseif flag == 0
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
      if h(6+4*n), res = 'melhora'; elseif deb, res = 'falha'; else, res = 'falha:P/2'; end
    else
      if h(6+4*n), res = 'aceite'; else, res = 'rejeitado'; end
    end
    fprintf('  %-10s %5d\n', res, h(7+4*n));
  end
end
end
