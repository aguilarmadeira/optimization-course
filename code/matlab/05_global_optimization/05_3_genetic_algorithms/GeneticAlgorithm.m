function [xb, fb, info] = GeneticAlgorithm(f, n, lb, ub, N, G, pc, pm, sigma, e, s)
%GENETICALGORITHM  Algoritmo genético com codificação real: min f(x) em X = [lb, ub].
%
%   [xb, fb, info] = GeneticAlgorithm(f, n, lb, ub, N, G, pc, pm, sigma, e)
%   [xb, fb, info] = GeneticAlgorithm(f, n, lb, ub, N, G, pc, pm, sigma, e, s)
%
%   Seleção por torneio de tamanho s, cruzamento BLX-0.5 (cruzarBLX, aos
%   pares, com probabilidade pc), mutação gaussiana por gene (probabilidade
%   pm, desvio sigma), projeção em X e elitismo (as e melhores soluções de
%   P_t substituem os e primeiros descendentes).
%
%   Entradas
%     f        função (handle) de um vetor linha x com n componentes
%     n        número de variáveis
%     lb, ub   limites de X (escalares ou vetores com n componentes)
%     N        tamanho da população
%     G        número de gerações
%     pc       probabilidade de cruzamento (por par)
%     pm       probabilidade de mutação (por gene)
%     sigma    desvio-padrão da mutação gaussiana
%     e        número de elites (0 <= e <= N)
%     s        (opcional, 2) tamanho do torneio; com s = 2 são as linhas
%              i1, i2 do slide
%
%   Saídas
%     xb       melhor solução da população final P_G
%     fb       f(xb)
%     info     estrutura com
%       .nfev      avaliações de f: N (G + 1) (população inicial + N por geração)
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       gerações feitas (= G)
%       .P, .F     população final e os seus valores de f
%       .history   uma linha por geração t = 0, 1, ..., G:
%                  [t  melhor f  média de f  x_melhor(1:n)]
%       .cols      nomes das colunas de .history
%       .flag      0 (fez as G gerações)
%       .message   mensagem de paragem
%
%   Como no slide, F guarda os valores de f (menor = melhor), não uma
%   aptidão a maximizar. Usa rand, randn e randi: para repetir uma corrida,
%   fixar a semente antes da chamada com rng(s).
%
%   Otimização — deck 5.3 (Algoritmos genéticos).
%   Reproduz os exemplos dos slides: ver ex05_3_genetic_algorithms.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. reavalia as elites em vez de guardar os
%   valores já calculados, não tem critério de paragem por estagnação nem
%   por diversidade, e não trata restrições além dos limites).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 11 || isempty(s), s = 2; end
lb = lb(:)';  ub = ub(:)';                  % escalares ou vetores linha

P = lb + (ub - lb).*rand(N, n);             % população inicial em X
F = arrayfun(@(i) f(P(i,:)), 1:N)';         % F: valores de f (menor = melhor)
nfev = N;
hist = zeros(G + 1, 3 + n);
[fb, ib] = min(F);  hist(1, :) = [0, fb, mean(F), P(ib, :)];

% --- núcleo (igual ao dos slides, com as contagens) ----------------------
for t = 1:G
  [~, o] = sort(F);  E = P(o(1:e), :);     % elites
  i1 = randi(N,N,1);                       % torneio (s candidatos por pai;
  for r = 2:s                              %  s = 2: i1, i2 como no slide)
    i2 = randi(N,N,1);  m = F(i2) < F(i1);
    i1(m) = i2(m);                         % fica o melhor
  end
  pais = P(i1,:);
  Q = cruzarBLX(pais, pc);                 % BLX-0.5
  M = rand(N,n) < pm;                      % mutacao
  Q = Q + M.*sigma.*randn(N,n);            % gaussiana
  Q = min(max(Q,lb),ub);  Q(1:e,:) = E;    % X; elites
  P = Q;  F = arrayfun(@(i) f(P(i,:)), 1:N)';
  nfev = nfev + N;                         % N avaliações por geração
  [fb, ib] = min(F);  hist(t + 1, :) = [t, fb, mean(F), P(ib, :)];
end
% -------------------------------------------------------------------------

xb = P(ib, :);
info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = G;
info.P = P;
info.F = F;
info.history = hist;
cols = {'t', 'melhor f', 'média f'};
for i = 1:n, cols{end + 1} = sprintf('x_%d', i); end
info.cols = cols;
info.flag = 0;
info.message = sprintf('fez as %d gerações pedidas: melhor f = %.4g', G, fb);
end
