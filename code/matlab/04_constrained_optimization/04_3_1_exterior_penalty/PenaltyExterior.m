function [x, fx, info] = PenaltyExterior(f, g, h, x0, R0, c, tolviol, tolx, tmax, interno, verbose)
%PENALTYEXTERIOR  Método da função de penalização exterior (EPFM).
%
%   [x, fx, info] = PenaltyExterior(f, g, h, x0, R0, c, tolviol, tolx, tmax)
%   [x, fx, info] = PenaltyExterior(..., tmax, interno)
%   [x, fx, info] = PenaltyExterior(..., tmax, interno, verbose)
%
%   min f(x) s.a. g_j(x) >= 0 (j = 1..J), h_l(x) = 0 (l = 1..K): minimiza
%   P(x,R) = f(x) + R [ sum_j <g_j(x)>^2 + sum_l h_l(x)^2 ], <g> = min(g,0),
%   para R crescente (R <- c R, c > 1), com arranque a quente.
%
%   Entradas
%     f        função objetivo (handle) de x
%     g        handle x -> valor ou vetor das desigualdades g_j(x) >= 0
%              (@(x) [] se não houver)
%     h        handle x -> valor ou vetor das igualdades h_l(x) = 0
%              (@(x) [] se não houver)
%     x0       ponto inicial x^(0) (não precisa de ser admissível)
%     R0       R^(0) > 0
%     c        fator de aumento de R, c > 1
%     tolviol  para no ciclo t se viol(x^(t)) <= tolviol e
%     tolx       ||x^(t) - x^(t-1)|| / max(1, ||x^(t-1)||) <= tolx
%     tmax     número máximo de ciclos exteriores
%     interno  (opcional) minimizador sem restrições do cap. 3, handle
%              [xn, Pn, out] = interno(P, x), com out.nfev = chamadas a P.
%              Por omissão, o do slide: fminsearch(P, x, opt), com
%              opt = optimset('TolX', 1e-10, 'TolFun', 1e-10,
%              'MaxFunEvals', 2000, 'MaxIter', 2000) e out.nfev = funcCount.
%              Os exemplos passam @(P, x) NelderMeadComp(P, x, 1e-10, 1e-12, 2000)
%              (NelderMeadComp está em code/matlab/common/; ver uc_setup).
%     verbose  (opcional, false) se true, imprime a tabela dos ciclos
%
%   Saídas
%     x        x^(t) do último ciclo
%     fx       f(x)
%     info     estrutura com
%       .nfev      TODAS as chamadas a f: chamadas a P nas minimizações
%                  internas (cada uma avalia f uma vez) + a avaliação final
%       .u         -2R<g_j(x)> (estimativa de u*, KKT)
%       .lambda    -2R h_l(x) (estimativa do multiplicador beta* das aulas, com L = f - beta h)
%       .R         R usado no último ciclo
%       .ngev, .nhev   0 (o minimizador interno por omissão só usa f)
%       .nit       ciclos exteriores t feitos
%       .history   uma linha por t = 0, 1, ..., nit (t = 0: ponto inicial):
%                  [t R n_f(ciclo) n_f(acum.) viol x_1..x_n u_1..u_J lambda_1..lambda_K]
%                  (na linha t, R é o valor usado para obter x^(t); na
%                  linha 0, R, u e lambda são NaN)
%       .cols      nomes das colunas de .history
%       .flag      0 (critério de paragem) ou 1 (atingiu tmax)
%       .message   mensagem de paragem
%
%   Otimização — deck 4.3.1 (Método da função de penalização exterior).
%   Reproduz os exemplos dos slides: ver ex04_3_1_exterior_penalty.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não põe as restrições em escalas comparáveis
%   nem adapta c ou a tolerância interna ao longo dos ciclos).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 10 || isempty(interno), interno = @interno_fminsearch; end
if nargin < 11, verbose = false; end

col = @(v) v(:);                          % valores de g e h em coluna
x = x0(:)';  n = numel(x);
J = numel(g(x));  K = numel(h(x));
br = @(z) min(z, 0);                      % operador <.>
viol0 = max([abs(br(col(g(x)))); abs(col(h(x))); 0]);
hist = [0, NaN, 0, 0, viol0, x, NaN(1, J + K)];
info.nfev = 0;
R = R0;  Rult = R;  flag = 1;

% --- núcleo (o dos slides, com as contagens) -----------------------------
for t = 1:tmax
  P = @(x) f(x) + R*(sum(br(g(x)).^2) + sum(h(x).^2));
  [xn, ~, out] = interno(P, x);           % a partir de x
  xn = xn(:)';
  info.nfev = info.nfev + out.nfev;
  viol = max([abs(br(col(g(xn)))); abs(col(h(xn))); 0]);
  Rult = R;
  hist(end + 1, :) = [t, R, out.nfev, info.nfev, viol, xn, ...
                      -2*R*br(col(g(xn)))', -2*R*col(h(xn))'];
  if viol <= tolviol && ...
     norm(xn - x)/max(1, norm(x)) <= tolx
    x = xn;  flag = 0;  break
  end
  x = xn;  R = c*R;                       % aumentar R
end
% -------------------------------------------------------------------------

fx = f(x);  info.nfev = info.nfev + 1;    % avaliação final
info.u = -2*Rult*br(col(g(x)));           % u_R = -2R<g>
info.lambda = -2*Rult*col(h(x));          % beta_R = -2R h
info.R = Rult;
info.ngev = 0;
info.nhev = 0;
info.nit = size(hist, 1) - 1;
info.history = hist;
cols = {'t', 'R', 'n_f(ciclo)', 'n_f(acum)', 'viol'};
for i = 1:n, cols{end + 1} = sprintf('x%d', i); end
for j = 1:J, cols{end + 1} = sprintf('u%d', j); end
for l = 1:K, cols{end + 1} = sprintf('lambda%d', l); end
info.cols = cols;
info.flag = flag;
if flag == 0
  info.message = sprintf('viol <= tolviol e passo relativo <= tolx no ciclo %d (R = %g)', ...
                         info.nit, Rult);
else
  info.message = sprintf('atingiu tmax = %d ciclos (R = %g)', tmax, Rult);
end

if verbose
  fprintf('%3s %9s %10s %9s %10s', cols{1:5});  fprintf(' %10s', cols{6:end});  fprintf('\n');
  for i = 1:size(hist, 1)
    if isnan(hist(i, 2)), Rs = sprintf('%9s', '-'); else, Rs = sprintf('%9.1e', hist(i, 2)); end
    fprintf('%3d %s %10d %9d %10.2e', hist(i, 1), Rs, hist(i, 3), hist(i, 4), hist(i, 5));
    for v = hist(i, 6:end)
      if isnan(v), fprintf(' %10s', '-'); else, fprintf(' %10.6f', v); end
    end
    fprintf('\n');
  end
end
end

% -------------------------------------------------------------------------
function [x, Px, out] = interno_fminsearch(P, x0)
% minimizador interno por omissão (o do slide): fminsearch, out.nfev = funcCount
opt = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 2000, 'MaxIter', 2000);
[x, Px, ~, o] = fminsearch(P, x0, opt);
out.nfev = o.funcCount;
end
