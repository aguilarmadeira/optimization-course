% EX03_2_1_RANDOM_SEARCH  Reproduz os exemplos do deck 3.2.1 (pesquisa aleatória).
%
%   1. f = x1^2 + x2^2 em [-5,5]^2, uma corrida (rng(1)): melhor f com
%      N = 100 e N = 1000 amostras (slide 1: 1.406 e 0.053).
%   2. Análise experimental: 30 corridas por N, mediana e quartis (slide 3).
%   3. Localizada no Rosenbrock de (-1.5, 2), m = 20, r0 = 1, gamma = 0.9:
%      uma corrida, 150 iterações, n_f = 3001 (slide 5);
%      30 corridas (sementes 0-29): n_f até f < 1e-4 = 975 [867; 1052].
%   4. Pura vs. localizada, orçamento 3000, 30 corridas (slide 6).
%
%   Método estocástico. Os números dos slides foram gerados em Python
%   (numpy.random.default_rng); o gerador do MATLAB/Octave é outro, por isso
%   as MESMAS sementes dão amostras DIFERENTES e os números não coincidem
%   um a um. A verificação aqui é ESTATÍSTICA (30 corridas, como no slide):
%     - pesquisa pura: a mediana das 30 corridas está a menos de um fator 3
%       da mediana do slide (a mediana de 30 corridas varia ~25 %);
%     - localizada, custo até f < 1e-4 (975): mediana a menos de 15 %;
%     - localizada, melhor f com n_f = 3001 (~1e-12): mesma ordem de
%       grandeza (fator 10) da mediana do slide;
%     - uma corrida isolada do slide (0.053; 5.9e-13) está entre o mínimo e
%       o máximo das 30 corridas feitas aqui.
%   A corrida do slide com N = 100 (f = 1.406) não entra na verificação: com
%   100 amostras, P(f_best > 1.406) ~ 1 %, é uma corrida invulgarmente má.
%   A reprodução exata está no exemplo Python (mesmas sementes, mesmo gerador).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 3.2.1.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
dentro = @(v, lo, hi) all(v(:) >= lo(:) & v(:) <= hi(:));
fator = @(v, s, k) v >= s/k && v <= s*k;    % v a menos de um fator k de s?
simnao = {'não', 'sim'};
% mediana e quartis [Q1 Q3] por interpolação linear (= numpy.percentile),
% sem depender da Statistics Toolbox
prctile30 = @(v) interp1(0:numel(v)-1, sort(v(:))', (numel(v) - 1)*[50 25 75]/100);
ok = true;

fprintf('Otimização — deck 3.2.1: pesquisa aleatória pura e localizada\n');
fprintf('(MATLAB/Octave: gerador diferente do Python dos slides -> verificação estatística)\n');

sq  = @(x) x(1)^2 + x(2)^2;
ros = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;

% ------------------------------------------------------------ 1. uma corrida
a = [-5 -5];  b = [5 5];
fprintf('\n1. f = x1^2 + x2^2 em [-5,5]^2, uma corrida (rng(1)), N = 1000\n');
[x, fx, info] = RandomSearch(sq, a, b, 1000, 1, [], true);
H = info.history;
f100 = H(find(H(:, 1) <= 100, 1, 'last'), end);   % melhor f com as primeiras 100 amostras
fprintf('N = 100: melhor f = %.3f;  N = 1000: melhor f = %.3f em (%.4f, %.4f);  n_f = %d\n', ...
        f100, fx, x(1), x(2), info.nfev);
fprintf('  (slide, Python: 1.406 e 0.053; é uma corrida — ver a verificação em 2.)\n');
c = info.nfev == 1000;
fprintf('  n_f = N: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 2. 30 corridas
fprintf('\n2. Análise experimental: 30 corridas por N (rng(0), ..., rng(29))\n');
Ns = [10 100 1000 10000];
% slide: ||x_best|| mediana [Q1; Q3]  ->  f = ||x||^2
Tx = [1.5 1.2 2.1;  0.54 0.38 0.76;  0.14 0.09 0.17;  0.053 0.027 0.070];
Tf = [2.3 0.29 0.018 0.0028];
fprintf('%6s %10s %10s %10s %10s %9s %22s\n', 'N', 'f (med)', '||x|| med', 'Q1', 'Q3', 'slide f', 'slide ||x|| [Q1; Q3]');
c = false(1, numel(Ns));  fb = zeros(30, numel(Ns));
for j = 1:numel(Ns)
  for s = 0:29
    [~, fb(s + 1, j)] = RandomSearch(sq, a, b, Ns(j), s);
  end
  q = prctile30(fb(:, j));  nx = sqrt(q);
  fprintf('%6d %10.4g %10.3g %10.3g %10.3g %9.2g %8.2g [%.2g; %.2g]\n', Ns(j), q(1), nx, Tf(j), Tx(j, :));
  c(j) = fator(q(1), Tf(j), 3);                    % mediana a menos de um fator 3 da do slide
end
fprintf('  mediana de f a menos de um fator 3 da do slide, N = 10, 100, 1000, 10000: %s | %s | %s | %s\n', simnao{c + 1});
ok = ok && all(c);
c = dentro(0.053, min(fb(:, 3)), max(fb(:, 3)));
fprintf('  corrida do slide com N = 1000 (f = 0.053) dentro de [mín; máx] = [%.4f; %.4f] das 30: %s\n', ...
        min(fb(:, 3)), max(fb(:, 3)), simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 3. localizada
x0 = [-1.5 2];
fprintf('\n3. Localizada no Rosenbrock de (-1.5, 2), m = 20, r0 = 1, gamma = 0.9\n');
[x, fx, info] = LocalRandomSearch(ros, x0, 1, 20, 0.9, 1e-8, 150, [], 3);
fprintf('uma corrida (rng(3)): %d iterações, n_f = %d, x = (%.6f, %.6f), f = %.2e\n', ...
        info.nit, info.nfev, x(1), x(2), fx);
fprintf('  primeiras e últimas linhas de info.history [k, n_f, r, f(x), x1, x2]:\n');
fprintf('  %5d %6d %10.3e %12.4e %10.6f %10.6f\n', info.history([1:4, end-1:end], :)');
c = [info.nit == 150, info.nfev == 3001, norm(x - [1 1]) < 1e-3];
fprintf('  150 it.: %s | n_f = 3001: %s | chega a (1,1): %s\n', simnao{c + 1});
ok = ok && all(c);

nhit = zeros(1, 30);
for s = 0:29
  [~, ~, info] = LocalRandomSearch(ros, x0, 1, 20, 0.9, 1e-8, 5000, 1e-4, s);
  nhit(s + 1) = info.nhit;
end
q = prctile30(nhit(isfinite(nhit)));
fprintf('30 corridas (rng(0)-rng(29)): n_f até f < 1e-4 = %.0f [%.0f; %.0f] (mediana [Q1; Q3]); falhas: %d\n', ...
        q, sum(isinf(nhit)));
fprintf('  (slide, Python: 975 [867; 1052])\n');
c = [abs(q(1) - 975) <= 0.15*975, all(isfinite(nhit))];
fprintf('  mediana a menos de 15%% de 975: %s | 30 de 30 atingem f < 1e-4: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 4. pura vs. localizada
fprintf('\n4. Pura vs. localizada no Rosenbrock, orçamento 3000, 30 corridas (rng(0)-rng(29))\n');
fp = zeros(1, 30);  fl = zeros(1, 30);
for s = 0:29
  [~, fp(s + 1)] = RandomSearch(ros, [-2 -1], [2 3], 3000, s);
  [~, fl(s + 1)] = LocalRandomSearch(ros, x0, 1, 20, 0.9, 1e-8, 150, [], s);
end
qp = prctile30(fp);  ql = prctile30(fl);
fprintf('pura em [-2,2]x[-1,3], n_f = 3000:  melhor f = %.3f [%.3f; %.3f]   (slide: 0.014 [0.008; 0.022])\n', qp);
fprintf('localizada, 150 it., n_f = 3001:    melhor f = %.1e [%.1e; %.1e]   (slide: 2e-12 [5e-13; 8e-12])\n', ql);
c = [fator(qp(1), 0.014, 3), fator(ql(1), 2e-12, 10), dentro(5.9e-13, min(fl), max(fl))];
fprintf('  mediana: pura a menos de um fator 3 do slide: %s | localizada a menos de um fator 10: %s\n', simnao{c(1:2) + 1});
fprintf('  corrida do slide (f = 5.9e-13) dentro de [mín; máx] das 30: %s\n', simnao{c(3) + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides (verificação estatística): %s\n', simnao{ok + 1});

