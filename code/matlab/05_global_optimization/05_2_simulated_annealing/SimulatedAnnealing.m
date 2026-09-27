function [x, fx, info] = SimulatedAnnealing(f, x0, T0, c, sigma, Nmax, lb, ub)
%SIMULATEDANNEALING  Simulated annealing: min f(x) em X = [lb, ub].
%
%   [x, fx, info] = SimulatedAnnealing(f, x0, T0, c, sigma, Nmax, lb, ub)
%
%   Um único ponto x; em cada iteração propõe-se um vizinho
%   x' = x + sigma*zeta, zeta ~ N(0, I), projetado em X, e aceita-se pelo
%   critério de Metropolis: sempre se Delta f <= 0; se Delta f > 0, com
%   probabilidade exp(-Delta f/T). Depois T = c*T (arrefecimento geométrico).
%   Devolve-se o melhor ponto visitado.
%
%   Entradas
%     f        função (handle) de um vetor x
%     x0       ponto inicial (escalar, linha ou coluna)
%     T0       temperatura inicial, T0 > 0 (na escala de f)
%     c        fator de arrefecimento, 0 < c < 1
%     sigma    desvio-padrão da perturbação gaussiana (passo)
%     Nmax     número de iterações (uma avaliação de f por iteração)
%     lb, ub   limites de X (escalares ou vetores com a forma de x0);
%              (opcionais, -Inf e Inf) o vizinho é projetado em X
%
%   Saídas
%     x        melhor ponto visitado (x_best)
%     fx       f(x_best)
%     info     estrutura com
%       .nfev      avaliações de f: Nmax + 1 (conta f(x0))
%       .nacc      propostas aceites
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       iterações feitas (= Nmax)
%       .T         temperatura no fim: T0*c^Nmax
%       .xk, .fk   ponto corrente final e f nesse ponto
%       .history   uma linha por k = 0, 1, ..., Nmax:
%                  [k  T  aceite  f(x')  f(x_k)  f_best  x'(1:n)  x_k(1:n)]
%                  (T: a temperatura usada na iteração k, T0*c^(k-1);
%                  x' é a proposta e x_k o ponto corrente depois da
%                  iteração k; na linha 0, x' = x0 e T = T0)
%       .cols      nomes das colunas de .history
%       .flag      0 (fez as Nmax iterações)
%       .message   mensagem de paragem
%
%   Usa randn e rand: para repetir uma corrida, fixar a semente antes da
%   chamada com rng(s).
%
%   Otimização — deck 5.2 (Simulated annealing).
%   Reproduz os exemplos dos slides: ver ex05_2_simulated_annealing.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não adapta sigma nem T ao longo da corrida,
%   não tem reaquecimentos nem critério de paragem além de Nmax).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 7 || isempty(lb), lb = -Inf; end
if nargin < 8 || isempty(ub), ub = Inf; end

x = x0;
n = numel(x);
hist = zeros(Nmax + 1, 6 + 2*n);

% --- núcleo (igual ao dos slides, com as contagens) ----------------------
fx = f(x);  nfev = 1;
xbest = x;  fbest = fx;  T = T0;  nacc = 0;
hist(1, :) = [0, T0, 1, fx, fx, fbest, x(:)', x(:)'];
for k = 1:Nmax
    xn = x + sigma*randn(size(x));              % vizinho
    xn = min(max(xn, lb), ub);                  % projetar em X
    fn = f(xn);  df = fn - fx;
    nfev = nfev + 1;                            % uma avaliação por iteração
    if df <= 0 || rand < exp(-df/T)             % Metropolis
        x = xn;  fx = fn;  aceite = 1;  nacc = nacc + 1;
        if fx < fbest, xbest = x; fbest = fx; end
    else
        aceite = 0;
    end
    hist(k + 1, :) = [k, T, aceite, fn, fx, fbest, xn(:)', x(:)'];
    T = c*T;                                    % arrefecer
end
% -------------------------------------------------------------------------

info.nfev = nfev;
info.nacc = nacc;
info.ngev = 0;
info.nhev = 0;
info.nit = Nmax;
info.T = T;
info.xk = x;
info.fk = fx;
info.history = hist;
cols = {'k', 'T', 'aceite', 'f(x'')', 'f(x_k)', 'f_best'};
for i = 1:n, cols{end + 1} = sprintf('x''_%d', i); end
for i = 1:n, cols{end + 1} = sprintf('x_k,%d', i); end
info.cols = cols;
info.flag = 0;
info.message = sprintf('fez as %d iterações pedidas: %d aceites, T final = %.3g', ...
                       Nmax, nacc, T);

x = xbest;
fx = fbest;
end
