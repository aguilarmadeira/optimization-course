function [X, F, info] = NSGA2(f, lb, ub, N, G, s, cv, par, verbose)
%NSGA2  NSGA-II: algoritmo genético multiobjetivo (rank + crowding, P U Q).
%
%   [X, F, info] = NSGA2(f, lb, ub, N, G)
%   [X, F, info] = NSGA2(f, lb, ub, N, G, s, cv, par, verbose)
%
%   Entradas
%     f        função (handle) de x, vetor linha com n componentes, que
%              devolve o vetor linha dos m objetivos [f1 ... fm]
%     lb, ub   limites da caixa [x^(L), x^(U)] (vetores linha)
%     N        tamanho da população (par)
%     G        número de gerações
%     s        (opcional, []) semente: se não for vazia, faz rng(s) antes de
%              gerar P_0; se for vazia, usa o estado atual do gerador
%     cv       (opcional, []) handle de x -> violação das restrições (>= 0;
%              0 = admissível), tratada pela regra de Deb; com cv = [], só a
%              caixa, garantida por projeção
%     par      (opcional) estrutura com os parâmetros (os que faltarem
%              tomam o valor por omissão):
%                .pc      probabilidade de cruzamento (0.9)
%                .eta_c   índice de distribuição do SBX (15)
%                .pm      probabilidade de mutação por variável (1/n)
%                .eta_m   índice de distribuição da mutação polinomial (20)
%                .ref     ponto de referência r do HV no histórico ([]: sem HV)
%     verbose  (opcional, false) se true, imprime o histórico de G/10 em
%              G/10 gerações
%
%   Saídas
%     X        soluções da frente F_1 de P_G (uma por linha)
%     F        objetivos dessas soluções
%     info     estrutura com
%       .nfev      avaliações do vetor f: N + G*N (a população inicial conta)
%       .ngev, .nhev   0 (o método só usa valores de f)
%       .nit       gerações (= G)
%       .P, .FP, .cvP  população final P_G, objetivos e violações
%       .rank, .crowd  rank e crowding em P_G
%       .history   uma linha por população P_t, t = 0..G:
%                  [t  n_f  |F_1|  inadm  HV]
%                  (|F_1| = não dominados admissíveis; HV da frente F_1
%                  com r = par.ref, NaN sem par.ref)
%       .cols      nomes das colunas de .history
%       .flag      0 (fez as G gerações)
%       .message   mensagem de paragem
%
%   Seleção para reprodução: torneio binário por <_n (rank; em empate, maior
%   crowding). Filhos: SBX (com prob. pc) e mutação polinomial, projetados
%   na caixa. Seleção ambiental (elitismo): R_t = P_t U Q_t, ordenado em
%   frentes; P_{t+1} enche-se frente a frente e a primeira que não cabe é
%   truncada pelos maiores d_p.
%
%   Método estocástico: o resultado é uma variável aleatória. Reporta-se a
%   mediana e os quartis de 20-30 corridas com sementes diferentes.
%
%   Otimização — deck 6.6 (NSGA-II), frame «O algoritmo».
%   Reproduz o exemplo dos slides: ver ex06_6_nsga2.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não trata valores de f que sejam NaN; o SBX
%   é a versão simples, com todas as variáveis cruzadas).
%   Com biblioteca: gamultiobj (Global Optimization Toolbox), ver
%   ex06_6_gamultiobj.m.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 6, s = []; end
if nargin < 7, cv = []; end
if nargin < 8, par = struct(); end
if nargin < 9, verbose = false; end
lb = lb(:)';  ub = ub(:)';  n = numel(lb);
if ~isfield(par, 'pc'),    par.pc = 0.9;     end
if ~isfield(par, 'eta_c'), par.eta_c = 15;   end
if ~isfield(par, 'pm'),    par.pm = 1/n;     end
if ~isfield(par, 'eta_m'), par.eta_m = 20;   end
if ~isfield(par, 'ref'),   par.ref = [];     end
if mod(N, 2) ~= 0, error('NSGA2: N tem de ser par.'); end
if ~isempty(s), rng(s); end

nfev = 0;
hist = zeros(G + 1, 5);

% --- núcleo (o do slide «O algoritmo», com as contagens) -----------------
P = repmat(lb, N, 1) + repmat(ub - lb, N, 1).*rand(N, n);   % P_0 na caixa
[FP, CP] = avaliar(f, cv, P);  nfev = nfev + N;             % avaliar f em P_0
for t = 0:G-1
  [fronts, rank] = NonDominatedSort(FP, cvarg(cv, CP));     % frentes de P_t
  d = zeros(N, 1);
  for k = 1:numel(fronts)
    d(fronts{k}) = CrowdingDistance(FP, fronts{k});         % d_p em cada frente
  end
  hist(t + 1, :) = registo(t, nfev, FP, CP, fronts, par.ref);

  Q = zeros(N, n);
  for k = 1:2:N                                             % N filhos
    i = randi(N, 1, 2);  j = randi(N, 1, 2);
    a = P(vence(i(1), i(2), rank, d), :);                   % torneio binário por <_n
    b = P(vence(j(1), j(2), rank, d), :);
    if rand < par.pc
      [c1, c2] = sbx(a, b, lb, ub, par.eta_c);
    else
      c1 = a;  c2 = b;
    end
    Q(k, :)     = mutar(c1, lb, ub, par.pm, par.eta_m);
    Q(k + 1, :) = mutar(c2, lb, ub, par.pm, par.eta_m);
  end
  [FQ, CQ] = avaliar(f, cv, Q);  nfev = nfev + N;           % avaliar Q_t

  R = [P; Q];  FR = [FP; FQ];  CR = [CP; CQ];
  fronts = NonDominatedSort(FR, cvarg(cv, CR));             % R_t = P_t U Q_t
  novos = [];
  for k = 1:numel(fronts)                                   % encher P_{t+1}
    fr = fronts{k};
    if numel(novos) + numel(fr) <= N
      novos = [novos, fr];                                  %#ok<AGROW>
    else                                                    % truncar pelos maiores d_p
      [~, o] = sort(-CrowdingDistance(FR, fr));
      novos = [novos, fr(o(1:N - numel(novos)))];           %#ok<AGROW>
      break
    end
  end
  P = R(novos, :);  FP = FR(novos, :);  CP = CR(novos);
end
% -------------------------------------------------------------------------

[fronts, rank] = NonDominatedSort(FP, cvarg(cv, CP));
d = zeros(N, 1);
for k = 1:numel(fronts)
  d(fronts{k}) = CrowdingDistance(FP, fronts{k});
end
hist(G + 1, :) = registo(G, nfev, FP, CP, fronts, par.ref);
F1 = fronts{1};
X = P(F1, :);
F = FP(F1, :);

info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = G;
info.P = P;  info.FP = FP;  info.cvP = CP;
info.rank = rank;  info.crowd = d;
info.history = hist;
info.cols = {'t', 'n_f', '|F_1|', 'inadm', 'HV'};
info.flag = 0;
info.message = sprintf('fez as %d gerações pedidas: %d não dominados em P_G', G, numel(F1));

if verbose
  passo = max(1, floor(G/10));
  fprintf('%5s %7s %7s %7s %9s\n', info.cols{:});
  for k = 1:size(hist, 1)
    if mod(hist(k, 1), passo) == 0 || hist(k, 1) == G
      fprintf('%5d %7d %7d %7d %9.4f\n', hist(k, :));
    end
  end
end
end

% ------------------------------------------------------------ auxiliares
function [FX, CX] = avaliar(f, cv, X)
% objetivos e violação das restrições, ponto a ponto
M = size(X, 1);
f1 = f(X(1, :));
FX = zeros(M, numel(f1));  FX(1, :) = f1;
for k = 2:M
  FX(k, :) = f(X(k, :));
end
CX = zeros(M, 1);
if ~isempty(cv)
  for k = 1:M
    CX(k) = cv(X(k, :));
  end
end
end

function c = cvarg(cv, C)
% sem restrições, NonDominatedSort usa a dominância pura
if isempty(cv), c = []; else, c = C; end
end

function w = vence(i, j, rank, d)
% torneio binário por <_n: menor rank; em empate, maior crowding
if rank(i) ~= rank(j)
  if rank(i) < rank(j), w = i; else, w = j; end
elseif d(i) >= d(j)
  w = i;
else
  w = j;
end
end

function [c1, c2] = sbx(a, b, lb, ub, eta_c)
% cruzamento binário simulado (Deb)
n = numel(a);
u = rand(1, n);
beta = zeros(1, n);
k = u <= 0.5;
beta(k) = (2*u(k)).^(1/(eta_c + 1));
beta(~k) = (1./(2 - 2*u(~k))).^(1/(eta_c + 1));
c1 = 0.5*((1 + beta).*a + (1 - beta).*b);
c2 = 0.5*((1 - beta).*a + (1 + beta).*b);
c1 = min(max(c1, lb), ub);
c2 = min(max(c2, lb), ub);
end

function y = mutar(x, lb, ub, pm, eta_m)
% mutação polinomial (Deb), cada variável com probabilidade pm
n = numel(x);
mk = rand(1, n) < pm;
u = rand(1, n);
delta = zeros(1, n);
k = u < 0.5;
delta(k) = (2*u(k)).^(1/(eta_m + 1)) - 1;
delta(~k) = 1 - (2 - 2*u(~k)).^(1/(eta_m + 1));
y = min(max(x + mk.*delta.*(ub - lb), lb), ub);
end

function row = registo(t, nfev, FX, CX, fronts, ref)
% linha do histórico: [t  n_f  |F_1|  inadm  HV]
F1 = fronts{1};
F1 = F1(CX(F1) <= 0);                        % não dominados admissíveis
if isempty(ref)
  hv = NaN;
elseif isempty(F1)
  hv = 0;
else
  hv = Hypervolume2D(FX(F1, :), ref);
end
row = [t, nfev, numel(F1), sum(CX > 0), hv];
end
