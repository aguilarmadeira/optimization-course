function [x, fx, info] = SimulatedAnnealing(f, x0, T0, alpha, n, Ts, sigma, lb, ub, nfmax)
%SIMULATEDANNEALING  Simulated annealing como nas aulas (Deb, 2012, cap. 6): min f(x) em X = [lb, ub].
%
%   [x, fx, info] = SimulatedAnnealing(f, x0, T0, alpha, n, Ts, sigma)
%   [x, fx, info] = SimulatedAnnealing(f, x0, T0, alpha, n, Ts, sigma, lb, ub, nfmax)
%
%   Passo 1: escolher x0, a temperatura inicial T0 (alta), o número n de
%            pontos aceites a cada temperatura, alpha (0 < alpha < 1; nas
%            aulas entre 0.5 e 0.99) e a temperatura mínima Ts; fazer t = 0.
%   Passo 2: gerar um vizinho x' = x_t + sigma*zeta, zeta ~ N(0, I),
%            projetado em X.
%   Passo 3: Delta f = f(x') - f(x_t). Se Delta f < 0, aceitar
%            (x_{t+1} = x', t = t + 1); senão, gerar r ~ U(0, 1): se
%            r < exp(-Delta f/T), aceitar; senão, voltar ao passo 2.
%   Passo 4: se T < Ts, terminar; senão, se mod(t, n) = 0, fazer
%            T = alpha*T. Voltar ao passo 2.
%
%   O contador t só avança quando se aceita: a T baixa (quase tudo
%   rejeitado) cada patamar custa muitas avaliações; por isso há também um
%   limite nfmax de avaliações (info.flag = 1). Devolve-se o ponto final x_t
%   (como nas aulas) e também o melhor visitado (info.xbest), que na
%   prática convém guardar.
%
%   Entradas
%     f        função (handle) de um vetor x
%     x0       ponto inicial (escalar, linha ou coluna)
%     T0       temperatura inicial, T0 > 0 (na escala de f)
%     alpha    fator de redução da temperatura, 0 < alpha < 1
%     n        número de pontos aceites a cada temperatura
%     Ts       temperatura mínima (termina quando T < Ts)
%     sigma    desvio-padrão da perturbação gaussiana (vizinhança)
%     lb, ub   limites de X (escalares ou vetores com a forma de x0);
%              (opcionais, -Inf e Inf) o vizinho é projetado em X
%     nfmax    (opcional, 100000) limite de avaliações de f
%
%   Saídas
%     x        ponto final x_t
%     fx       f(x_t)
%     info     estrutura com
%       .xbest, .fbest   melhor ponto visitado e f nesse ponto
%       .nfev      avaliações de f (conta f(x0))
%       .nacc      pontos aceites (= t final); .nit = .nacc
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .T         temperatura no fim; .nT reduções de temperatura feitas
%       .history   uma linha por avaliação:
%                  [n_f  t  T  aceite  f(x')  f(x_t)  f_best  x'(1:n)  x_t(1:n)]
%                  (a linha 1 é o ponto inicial; T é a temperatura usada)
%       .cols      nomes das colunas de .history
%       .flag      0 (T < Ts) ou 1 (atingiu nfmax)
%       .message   mensagem de paragem
%
%   Usa randn e rand: para repetir uma corrida, fixar a semente antes da
%   chamada com rng(s).
%
%   Otimização — deck 5.2 (Simulated annealing).
%   Reproduz os exemplos dos slides: ver ex05_2_simulated_annealing.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não adapta sigma nem T ao longo da corrida,
%   não tem reaquecimentos).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 8 || isempty(lb), lb = -Inf; end
if nargin < 9 || isempty(ub), ub = Inf; end
if nargin < 10 || isempty(nfmax), nfmax = 100000; end
x = x0;  nv = numel(x);
cap = 1024;  hist = zeros(cap, 7 + 2*nv);  nh = 1;      % cresce por duplicação

% --- núcleo (igual ao dos slides, com as contagens) ----------------------
fx = f(x);  nfev = 1;                                    % passo 1
xbest = x;  fbest = fx;  T = T0;  t = 0;  nT = 0;  flag = 1;
hist(1, :) = [1, 0, T0, 1, fx, fx, fbest, x(:)', x(:)'];
while nfev < nfmax
    xn = min(max(x + sigma*randn(size(x)), lb), ub);    % passo 2: vizinho em X
    fn = f(xn);  df = fn - fx;  nfev = nfev + 1;
    Tk = T;
    if df < 0 || rand < exp(-df/T)                       % passo 3: Metropolis
        x = xn;  fx = fn;  t = t + 1;  aceite = 1;
        if fx < fbest, xbest = x; fbest = fx; end
    else
        aceite = 0;                                      % rejeitado: passo 2
    end
    nh = nh + 1;
    if nh > cap, hist = [hist; zeros(cap, 7 + 2*nv)]; cap = 2*cap; end %#ok<AGROW>
    hist(nh, :) = [nfev, t, Tk, aceite, fn, fx, fbest, xn(:)', x(:)'];
    if aceite                                            % passo 4
        if T < Ts, flag = 0; break, end
        if mod(t, n) == 0, T = alpha*T; nT = nT + 1; end
    end
end
% -------------------------------------------------------------------------

info.xbest = xbest;
info.fbest = fbest;
info.nfev = nfev;
info.nacc = t;
info.nit = t;
info.ngev = 0;
info.nhev = 0;
info.T = T;
info.nT = nT;
info.history = hist(1:nh, :);
xc = arrayfun(@(i) sprintf('x''_%d', i), 1:nv, 'UniformOutput', false);
xk = arrayfun(@(i) sprintf('x_t,%d', i), 1:nv, 'UniformOutput', false);
info.cols = [{'n_f', 't', 'T', 'aceite', 'f(x'')', 'f(x_t)', 'f_best'}, xc, xk];
info.flag = flag;
if flag == 0
  info.message = sprintf('T = %.3g < Ts: %d aceites, %d avaliações', T, t, nfev);
else
  info.message = sprintf('atingiu nfmax = %d avaliações (T = %.3g, %d aceites)', nfmax, T, t);
end
end
