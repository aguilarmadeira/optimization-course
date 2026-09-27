% EX06_6_NSGA2  Reproduz os exemplos do deck 6.6 (NSGA-II) e o hipervolume do 6.2.
%
%   1. Oito soluções (frames «rank» e «crowding»): F_1 = {1,2,4,8},
%      F_2 = {3,5,6,7}; d_2 = 1.4, d_4 = 1.233; só cabem 3: sai o 4.
%   2. Hipervolume (deck 6.2): três pontos com r = (1.1; 1.3): 0.549;
%      convergência/diversidade com r = (1.1; 1.65): 0.654, 0.656, 1.095;
%      frentes exatas com r = (1.1; 1.65): convexa 1.648, não convexa 1.148.
%   3. NSGA-II no exemplo convexo da UC, N = 40, G = 50: geração 0 com
%      6 não dominados, 2040 avaliações, HV 1.637.
%   4. Tabela «DMS vs. NSGA-II» (exemplo não convexo): coluna do NSGA-II,
%      20 corridas, mediana e quartis do HV por orçamento.
%
%   Partes 1 e 2: determinísticas, comparação exata com os slides.
%   Partes 3 e 4: método estocástico. Os números dos slides foram gerados em
%   Python (numpy.random.default_rng); o gerador do MATLAB/Octave é outro,
%   por isso as MESMAS sementes dão populações DIFERENTES e os números não
%   coincidem um a um. A verificação aqui é ESTATÍSTICA:
%     - n_f = N + G*N = 2040 em todas as corridas (não depende do gerador);
%     - HV final: a mediana de 20 corridas (rng(0), ..., rng(19)) está a
%       menos de 0.005 de 1.637, e 1.637 está entre o mínimo e o máximo;
%     - geração 0: 6 não dominados está entre o mínimo e o máximo de 200
%       populações iniciais (rng(0), ..., rng(199));
%     - tabela: em cada orçamento, a mediana de 20 corridas (rng(100), ...,
%       rng(119)) está a menos de max(0.005, (Q3 - Q1)/2) da mediana do
%       slide (Q1, Q3 do slide: com 100 avaliações a dispersão é grande).
%   A reprodução exata está no exemplo Python (mesmas sementes, mesmo gerador).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 6.6.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
dentro = @(v, lo, hi) all(v(:) >= lo(:) & v(:) <= hi(:));
simnao = {'não', 'sim'};
lista = @(v) regexprep(sprintf('%d,', v), ',$', '');   % [1 2 4] -> '1,2,4'
% mediana e quartis [Q2 Q1 Q3] por interpolação linear (= numpy.percentile),
% sem depender da Statistics Toolbox
quartis = @(v) interp1(0:numel(v)-1, sort(v(:))', (numel(v) - 1)*[50 25 75]/100);
ok = true;

fprintf('Otimização — deck 6.6: NSGA-II (e hipervolume, deck 6.2)\n');
fprintf('(partes 3 e 4: gerador diferente do Python dos slides -> verificação estatística)\n');

% ------------------------------------------------ 1. oito soluções
Pex = [1 6; 2 3; 3 5; 4 1.5; 5 4; 6 2.5; 2.5 6.5; 7 1];
fprintf('\n1. Oito soluções (f1, f2): rank e crowding\n');
[fronts, rank] = NonDominatedSort(Pex);
F1 = fronts{1};
d = CrowdingDistance(Pex, F1);
dp = nan(8, 1);  dp(F1) = d;
fprintf('%3s %12s %5s %9s\n', 'p', '(f1; f2)', 'rank', 'd_p');
for i = 1:8
  if isnan(dp(i)), ds = '--'; else, ds = sprintf('%.3f', dp(i)); end
  fprintf('%3d %12s %5d %9s\n', i, sprintf('(%g; %g)', Pex(i, :)), rank(i), ds);
end
fprintf('F_1 = {%s}, F_2 = {%s}\n', lista(sort(F1)), lista(sort(fronts{2})));
[~, o] = sort(-d);                          % só cabem 3 dos 4
fica = F1(o(1:3));
sai = setdiff(F1, fica);
fprintf('Só cabem 3 dos 4 de F_1: ficam {%s}, sai {%s}\n', lista(sort(fica)), lista(sai));
c = [isequal(F1, [1 2 4 8]) && isequal(sort(fronts{2}), [3 5 6 7]), ...
     confere(dp(2), 1.4, 1), confere(dp(4), 1.233, 3), all(isinf(dp([1 8]))), ...
     isequal(sort(fica), [1 2 8]) && isequal(sai, 4)];
fprintf('  frentes: %s | d_2 = 1.4: %s | d_4 = 1.233: %s | extremos d = Inf: %s | sai o 4: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 2. hipervolume (6.2)
fprintf('\n2. Hipervolume (deck 6.2)\n');
P3 = [0.2 0.96; 0.5 0.75; 0.8 0.36];
hv3 = Hypervolume2D(P3, [1.1 1.3]);
REF = [1.1 1.65];
ff = @(s) [s(:), 1 - s(:).^2];              % frente não convexa f2 = 1 - f1^2
A = ff(linspace(0.40, 0.55, 10));           % perto, concentrados
B = ff(linspace(0, 1, 10)) + repmat([0.12 0.25], 10, 1);   % espalhados, longe
C = ff(linspace(0, 1, 10));                 % perto e espalhados
hvA = Hypervolume2D(A, REF);  hvB = Hypervolume2D(B, REF);  hvC = Hypervolume2D(C, REF);
t = linspace(0, 1, 2000)';
hv_cx = Hypervolume2D([t.^2, (1 - t).^2], REF);   % frente convexa, 2000 pontos
hv_nc = Hypervolume2D(ff(t), REF);
hv_cx_ex = 1.1*1.65 - 1/6;                  % exato: r1*r2 - int_0^1 (1 - sqrt(f1))^2 df1
hv_nc_ex = 1.1*1.65 - 2/3;                  % exato: r1*r2 - int_0^1 (1 - f1^2) df1
fprintf('três pontos, r = (1.1; 1.3): HV = %.3f\n', hv3);
fprintf('r = (1.1; 1.65): perto/concentrados %.3f, espalhados/longe %.3f, perto/espalhados %.3f\n', ...
        hvA, hvB, hvC);
fprintf('frente exata convexa: %.4f (2000 pontos: %.4f); não convexa: %.4f (2000 pontos: %.4f)\n', ...
        hv_cx_ex, hv_cx, hv_nc_ex, hv_nc);
c = [confere(hv3, 0.549, 3), confere(hvA, 0.654, 3), confere(hvB, 0.656, 3), confere(hvC, 1.095, 3), ...
     confere([hv_cx_ex hv_cx], 1.648, 3), confere([hv_nc_ex hv_nc], 1.148, 3)];
fprintf('  0.549: %s | 0.654: %s | 0.656: %s | 1.095: %s | 1.648: %s | 1.148: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 3. NSGA-II, exemplo convexo
f  = @(x) [x(1)^2 + x(2)^2, (x(1) - 1)^2 + x(2)^2];
cv = @(x) max(0, x(1)^2 + x(2)^2 - 1);     % x1^2 + x2^2 <= 1 (regra de Deb)
par.ref = REF;
fprintf('\n3. NSGA-II no exemplo convexo: N = 40, G = 50, r = (1.1; 1.65)\n');
fprintf('Uma corrida (rng(4); no Python, semente 4 dá o slide):\n');
[X, F, info] = NSGA2(f, [0 0], [1 1], 40, 50, 4, cv, par, true);
fprintf('%s; n_f = %d, HV = %.4f (exata %.4f)\n', info.message, info.nfev, info.history(end, 5), hv_cx_ex);
hv = zeros(20, 1);  nf = zeros(20, 1);
for s = 0:19
  [~, ~, in] = NSGA2(f, [0 0], [1 1], 40, 50, s, cv, par);
  hv(s + 1) = in.history(end, 5);  nf(s + 1) = in.nfev;
end
q = quartis(hv);
g0 = zeros(200, 1);
for s = 0:199                               % só a população inicial (G = 0)
  [~, ~, in] = NSGA2(f, [0 0], [1 1], 40, 0, s, cv);
  g0(s + 1) = in.history(1, 3);
end
fprintf('20 corridas: HV mediana %.4f [Q1 %.4f; Q3 %.4f], mín %.4f, máx %.4f (slide: 1.637, uma corrida)\n', ...
        q, min(hv), max(hv));
fprintf('200 populações iniciais: não dominados admissíveis entre %d e %d, mediana %g (slide: 6)\n', ...
        min(g0), max(g0), median(g0));
c = [all(nf == 2040), abs(q(1) - 1.637) <= 0.005, dentro(1.637, min(hv), max(hv)), ...
     dentro(6, min(g0), max(g0))];
fprintf('  n_f = 2040 (todas): %s | mediana a menos de 0.005 de 1.637: %s | 1.637 em [mín; máx]: %s | 6 em [mín; máx]: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 4. DMS vs. NSGA-II
fn = @(x) [x(1), 1 - x(1)^2 + x(2)];       % exemplo não convexo, x em [0,1]x[0,0.6]
fprintf('\n4. DMS vs. NSGA-II no exemplo não convexo (r = (1.1; 1.65)), coluna do NSGA-II\n');
BUDGET = 2000;  Gt = BUDGET/40 - 1;         % 49 gerações: 40 + 49*40 = 2000
HV = zeros(20, Gt + 1);
for s = 100:119
  [~, ~, in] = NSGA2(fn, [0 0], [1 0.6], 40, Gt, s, [], par);
  HV(s - 99, :) = in.history(:, 5)';
end
ev = in.history(:, 2);                      % n_f de cada geração: 40, 80, ..., 2000
orc = [100 200 400 800 2000];
dms = [1.064 1.094 1.123 1.139 1.144];      % coluna do DMS (slide; ver 06_5_dms)
Tn = [0.986 0.968 1.006; 1.083 1.066 1.096; 1.125 1.117 1.129; ...
      1.132 1.131 1.133; 1.132 1.132 1.133];  % slide: mediana [Q1; Q3]
fprintf('%8s %8s %10s %18s %24s\n', 'n_f', 'DMS', 'NSGA med.', '[Q1; Q3]', 'slide: med. [Q1; Q3]');
c = false(1, numel(orc));
for e = 1:numel(orc)
  k = find(ev <= orc(e), 1, 'last');        % última geração completa que não excede
  q = quartis(HV(:, k));
  fprintf('%8d %8.3f %10.3f   [%.3f; %.3f] %12.3f [%.3f; %.3f]\n', orc(e), dms(e), q, Tn(e, :));
  c(e) = abs(q(1) - Tn(e, 1)) <= max(0.005, (Tn(e, 3) - Tn(e, 2))/2);
end
fprintf('  mediana a menos de max(0.005, (Q3-Q1)/2) da do slide, n_f = 100, 200, 400, 800, 2000: %s | %s | %s | %s | %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides (verificação estatística nas partes 3 e 4): %s\n', simnao{ok + 1});
