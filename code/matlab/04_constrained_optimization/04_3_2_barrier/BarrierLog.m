function [x, fx, info] = BarrierLog(f, g, x0, R0, c, tolcomp, tolx, tmax, interno, verbose)
%BARRIERLOG  Método da barreira logarítmica (penalização interior).
%
%   [x, fx, info] = BarrierLog(f, g, x0, R0, c, tolcomp, tolx, tmax)
%   [x, fx, info] = BarrierLog(..., tmax, interno)
%   [x, fx, info] = BarrierLog(..., tmax, interno, verbose)
%
%   min f(x) s.a. g_j(x) >= 0 (j = 1..J), a partir de x0 com todos os
%   g_j(x0) > 0: minimiza P(x,R) = f(x) - R sum_j ln g_j(x) (Pbar; Inf se
%   algum g_j <= 0) para R decrescente (R <- c R, 0 < c < 1), com arranque
%   a quente.
%
%   Entradas
%     f        função objetivo (handle) de x
%     g        handle x -> valor ou vetor das desigualdades g_j(x) >= 0
%     x0       ponto inicial ESTRITAMENTE admissível (g_j(x0) > 0)
%     R0       R^(0) > 0
%     c        fator de redução de R, 0 < c < 1
%     tolcomp  para no ciclo t se J*R <= tolcomp (J*R = sum(u.*g),
%              complementaridade perturbada) e
%     tolx       ||x^(t) - x^(t-1)|| / max(1, ||x^(t-1)||) <= tolx
%     tmax     número máximo de ciclos exteriores
%     interno  (opcional) minimizador sem restrições do cap. 3, handle
%              [xn, Pn, out] = interno(P, x), com out.nfev = chamadas a P;
%              tem de aceitar P = Inf. Por omissão, o do slide:
%              fminsearch(P, x, opt), com opt = optimset('TolX', 1e-10,
%              'TolFun', 1e-10, 'MaxFunEvals', 2000, 'MaxIter', 2000) e
%              out.nfev = funcCount. Os exemplos passam
%              @(P, x) NelderMeadComp(P, x, 1e-10, 1e-12, 2000)
%              (NelderMeadComp está em code/matlab/common/; ver uc_setup).
%     verbose  (opcional, false) se true, imprime a tabela dos ciclos
%
%   Saídas
%     x        x^(t) do último ciclo (estritamente admissível)
%     fx       f(x)
%     info     estrutura com
%       .nfev      chamadas a P nos ciclos + 1 (f(x) final), como no slide:
%                  majorante de n_f
%       .nPev      chamadas a P nos ciclos (sem a avaliação final)
%       .nfora     chamadas a P com x fora de X (avaliam g, mas não f)
%       .nfev_efetivo  avaliações de f de facto: nPev - nfora + 1
%       .u         R./g_j(x) (estimativa de u*; u_j g_j = R)
%       .R         R usado no último ciclo
%       .ngev, .nhev   0 (o minimizador interno por omissão só usa P)
%       .nit       ciclos exteriores t feitos
%       .history   uma linha por t = 0, 1, ..., nit (t = 0: ponto inicial):
%                  [t R P(ciclo) fora(ciclo) P(acum.) x_1..x_n u_1..u_J]
%                  (na linha t, R é o valor usado para obter x^(t); na
%                  linha 0, R e u são NaN)
%       .cols      nomes das colunas de .history
%       .flag      0 (critério de paragem) ou 1 (atingiu tmax)
%       .message   mensagem de paragem
%
%   Otimização — deck 4.3.2 (Método da função de penalização interior).
%   Reproduz os exemplos dos slides: ver ex04_3_2_barrier.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não tem Fase I: exige um ponto estritamente
%   admissível; não trata igualdades).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 9 || isempty(interno), interno = @interno_fminsearch; end
if nargin < 10, verbose = false; end

x = x0(:)';  n = numel(x);
g0 = g(x);  J = numel(g0);
if any(g0 <= 0)
  error('BarrierLog: x0 tem de ser estritamente admissível (g_j(x0) > 0).');
end
hist = [0, NaN, 0, 0, 0, x, NaN(1, J)];
Pbar('reset');                              % contador das chamadas fora de X
info.nPev = 0;
R = R0;  Rult = R;  flag = 1;

% --- núcleo (o dos slides, com as contagens) -----------------------------
for t = 1:tmax
  P = @(x) Pbar(f, g, x, R);                % Inf fora de X
  fora0 = Pbar('fora');
  [xn, ~, out] = interno(P, x);             % a partir de x
  xn = xn(:)';
  info.nPev = info.nPev + out.nfev;
  u = R ./ g(xn);                           % u_j = R/g_j
  Rult = R;
  hist(end + 1, :) = [t, R, out.nfev, Pbar('fora') - fora0, info.nPev, xn, u(:)'];
  if numel(u)*R <= tolcomp && ...           % J*R = sum(u.*g)
     norm(xn - x)/max(1, norm(x)) <= tolx
    x = xn;  flag = 0;  break
  end
  x = xn;  R = c*R;                         % reduzir R
end
% -------------------------------------------------------------------------

fx = f(x);                                  % avaliação final
info.nfev = info.nPev + 1;                  % como no slide: chamadas a P + 1
info.nfora = Pbar('fora');
info.nfev_efetivo = info.nPev - info.nfora + 1;
gx = g(x);
info.u = Rult ./ gx(:);
info.R = Rult;
info.ngev = 0;
info.nhev = 0;
info.nit = size(hist, 1) - 1;
info.history = hist;
cols = {'t', 'R', 'P(ciclo)', 'fora(ciclo)', 'P(acum)'};
for i = 1:n, cols{end + 1} = sprintf('x%d', i); end
for j = 1:J, cols{end + 1} = sprintf('u%d', j); end
info.cols = cols;
info.flag = flag;
if flag == 0
  info.message = sprintf('J*R <= tolcomp e passo relativo <= tolx no ciclo %d (R = %g)', ...
                         info.nit, Rult);
else
  info.message = sprintf('atingiu tmax = %d ciclos (R = %g)', tmax, Rult);
end

if verbose
  fprintf('%3s %9s %9s %11s %8s', cols{1:5});  fprintf(' %10s', cols{6:end});  fprintf('\n');
  for i = 1:size(hist, 1)
    if isnan(hist(i, 2)), Rs = sprintf('%9s', '-'); else, Rs = sprintf('%9.1e', hist(i, 2)); end
    fprintf('%3d %s %9d %11d %8d', hist(i, 1), Rs, hist(i, 3), hist(i, 4), hist(i, 5));
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
