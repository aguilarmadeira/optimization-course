function [x, fx, info] = PowellQuadratic(f, x1, Delta, tolx, tolf, kmax, verbose)
%POWELLQUADRATIC  Interpolação quadrática sucessiva (método de Powell, 1D).
%
%   [x, fx, info] = PowellQuadratic(f, x1, Delta, tolx, tolf, kmax)
%   [x, fx, info] = PowellQuadratic(f, x1, Delta, tolx, tolf, kmax, verbose)
%
%   Entradas
%     f        função (handle)
%     x1       ponto inicial
%     Delta    passo inicial (x2 = x1 + Delta)
%     tolx     tolerância em x:  |X_min - xb| <= tolx
%     tolf     tolerância em f:  |F_min - f(xb)| <= tolf   (as duas juntas)
%     kmax     (opcional, 100) número máximo de iterações
%     verbose  (opcional, false) se true, imprime a tabela das iterações
%
%   Saídas
%     x, fx    o melhor dos quatro pontos da última iteração e o seu valor
%     info     estrutura com
%       .nfev      avaliações de f: 3 iniciais + 1 por iteração + 1 por reflexão
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       número de iterações (cada uma calcula um xb)
%       .nrefl     número de reflexões (salvaguarda para a2 <= 0)
%       .history   tabela das iterações, uma linha por iteração k = 0, 1, ...:
%                  [k  x1  x2  x3  f1  f2  f3  a1  a2  xb  f(xb)]
%       .cols      nomes das colunas de .history
%       .flag      0 (critério satisfeito) ou 1 (atingiu kmax)
%       .message   mensagem de paragem
%
%   Salvaguardas incluídas: se a2 <= 0 (parábola sem mínimo), reflete-se o
%   pior ponto no melhor, 2*X_min - x_pior, e repete-se; depois de cada
%   iteração guardam-se o melhor ponto e os dois que o enquadram (se não
%   existirem, os três melhores).
%
%   Otimização — deck 2.3 (Interpolação quadrática sucessiva).
%   Reproduz o exemplo dos slides: ver ex02_3_powell.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. limitar o passo, pontos repetidos).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 6 || isempty(kmax), kmax = 100; end
if nargin < 7, verbose = false; end

% --- arranque: 3 avaliações ----------------------------------------------
x2 = x1 + Delta;
f1 = f(x1);  f2 = f(x2);
if f1 > f2
  x3 = x1 + 2*Delta;           % f desce: o terceiro ponto vai para a frente
else
  x3 = x1 - Delta;             % f sobe: o terceiro ponto vai para trás
end
f3 = f(x3);
nfev = 3;
X = [x1, x2, x3];  F = [f1, f2, f3];

hist = zeros(0, 11);
k = 0;  nrefl = 0;  flag = 1;
while k < kmax
  % --- núcleo de uma iteração (igual ao dos slides) ----------------------
  [X, i] = sort(X);  F = F(i);         % x1 < x2 < x3
  [Fmin, imin] = min(F);  Xmin = X(imin);
  a1 = (F(2)-F(1)) / (X(2)-X(1));
  a2 = ((F(3)-F(1))/(X(3)-X(1)) - a1) / (X(3)-X(2));
  if a2 <= 0          % parabola sem minimo:
    % refletir o pior ponto no melhor e repetir
    [~, iw] = max(F);
    X(iw) = 2*Xmin - X(iw);  F(iw) = f(X(iw));
    nfev = nfev + 1;  nrefl = nrefl + 1;
    if nrefl > kmax, flag = 1; break; end
    continue
  else
    xb = (X(1)+X(2))/2 - a1/(2*a2);
    fb = f(xb);                        % 1 avaliacao
    nfev = nfev + 1;
  end
  % -----------------------------------------------------------------------
  hist(end + 1, :) = [k, X, F, a1, a2, xb, fb]; %#ok<AGROW>
  k = k + 1;

  % parar se xb e fb mudam pouco
  if abs(Fmin - fb) <= tolf && abs(Xmin - xb) <= tolx
    flag = 0;
    break
  end

  % senao, guardar o melhor ponto e os dois que o enquadram
  [Xa, i] = sort([X, xb]);  Fa = [F, fb];  Fa = Fa(i);
  [~, ib] = min(Fa);
  if ib > 1 && ib < 4
    X = Xa(ib-1:ib+1);  F = Fa(ib-1:ib+1);
  else                                 % sem enquadramento: os três melhores
    [~, j] = sort(Fa);
    X = Xa(j(1:3));  F = Fa(j(1:3));
  end
end

% o melhor dos pontos disponíveis (na paragem: os quatro da última iteração)
if k > 0
  Xa = [X, xb];  Fa = [F, fb];
else
  Xa = X;  Fa = F;
end
[fx, ib] = min(Fa);  x = Xa(ib);

info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = k;
info.nrefl = nrefl;
info.history = hist;
info.cols = {'k', 'x1', 'x2', 'x3', 'f1', 'f2', 'f3', 'a1', 'a2', 'xb', 'f(xb)'};
info.flag = flag;
if flag == 0
  info.message = sprintf('|F_min - f(xb)| <= tolf e |X_min - xb| <= tolx em %d iterações', k);
else
  info.message = sprintf('atingiu o número máximo de iterações (kmax = %d)', kmax);
  warning('PowellQuadratic: %s', info.message);
end

if verbose
  fprintf('%3s %9s %9s %9s %10s %10s %10s %9s %8s %9s %10s\n', info.cols{:});
  for r = 1:size(hist, 1)
    fprintf('%3d %9.4f %9.4f %9.4f %10.4f %10.4f %10.4f %9.3f %8.3f %9.4f %10.4f\n', hist(r, :));
  end
end
end
